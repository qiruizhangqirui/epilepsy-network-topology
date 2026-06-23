% =========================================================================
% Script: S5_hubmess_pipeline.m
% -------------------------------------------------------------------------
% Description:
%   Pipeline for "Hubness" map analysis.
%
%   1. ROI extraction from hubness NIfTI maps.
%   2. Normative modeling (PCNtoolkit) for W-scores.
%   3. Max-T permutation testing and statistical mapping.
%
% Note:
%   - Requires 'permuztest': https://github.com/mickcrosse/PERMUTOOLS
%   - Requires 'ENIGMA' toolbox (plot_cortical, parcel_to_surface): https://github.com/MICA-MNI/ENIGMA
%   - Requires PCNtoolkit setup in Python.
%
% Inputs:
%   - outputs/All_derivatives_struct.mat
%   - data/Subjects.xlsx
%   - data/SPARK_Derivatives/
%
% Outputs:
%   - outputs/hubness_statistic/
%
% Author: Qirui Zhang,  Farber Institute for Neuroscience, Department of Neurology, Thomas Jefferson University
% Email: qirui.zhang@jefferson.edu;fmrizhangqr@126.com
% Date: 12/24/2025
% =========================================================================

clear; clc;

%% =========================
% Part 1: Configuration & Inputs
% =========================
this_file = mfilename('fullpath');

if ~isempty(this_file)
    % Case 1: script/function is being executed from a .m file
    project_root = fileparts(fileparts(fileparts(this_file)));
else
    % Case 2: executed from Command Window or interactive context
    project_root = fileparts(fileparts(pwd));
end

data_dir     = fullfile(project_root, 'data');
result_dir   = fullfile(project_root, 'outputs');

if ~exist(result_dir, 'dir'); mkdir(result_dir); end
addpath(genpath(fullfile(project_root, 'src', 'utils')));

% Add src/utils to Python path for pcn_wscore.py
if count(py.sys.path, '') == 0
    insert(py.sys.path, int32(0), '');
end
P = py.sys.path;
utils_py_path = fullfile(project_root, 'src', 'utils');
if count(P, utils_py_path) == 0
    insert(P, int32(0), utils_py_path);
end

pathDerivatives = fullfile(result_dir, 'All_derivatives_struct.mat');
pathSubjects    = fullfile(data_dir, 'Subjects.xlsx');

out_dir = fullfile(result_dir, 'hubness_statistic');
if ~exist(out_dir,'dir'); mkdir(out_dir); end

pathOutMeasures = fullfile(result_dir, 'Hubness_Measures.mat');
pathOutWscore   = fullfile(result_dir, 'Hubness_Wscore.mat');
pathOutStats    = fullfile(result_dir, 'Hubness_stats.mat');
SPARK_dir       = fullfile(data_dir, 'SPARK_Derivatives');

% Volume atlases for ROI extraction
AtlasNii = struct();
AtlasNii.AS200K17_200ROI = fullfile(data_dir, 'NCT_atlases', 'atlases', 'FSLMNI2mm', 'YeoLab', 'AS200Y17.nii.gz');
AtlasNii.AS200K34        = fullfile(data_dir, 'AS17network_34ROIs_sorted.nii.gz');
AtlasNii.Subcortical     = fullfile(data_dir, 'Subcortical_FLS2mm.nii');

% 34->200 mapping
assign34_mat = fullfile(data_dir, 'assignment_34.mat');

%% =========================
% Part 2: Load Data
% =========================
D = load(pathDerivatives);
All_derivatives_struct = D.All_derivatives_struct;

% Use any subject as a container for template network names
any_sub = 'JLH0001';
if ~isfield(All_derivatives_struct, any_sub)
    fns = fieldnames(All_derivatives_struct);
    any_sub = fns{1};
end
Atlas_network_name = All_derivatives_struct.(any_sub);

subjectTable = readtable(pathSubjects, 'VariableNamingRule','preserve');

subjectIDs = cellstr(string(subjectTable.("ID_new")));
site       = cellstr(string(subjectTable.("site")));
group_raw  = cellstr(string(subjectTable.("Group")));

sex = double(subjectTable.("sex_b(M1F0)"));
age = double(subjectTable.("age"));
hm  = double(subjectTable.("mean_FD"));

Nsub = numel(subjectIDs);

%% =========================
% Part 3: Group Logic (Refine TLE)
% =========================
group_ref = group_raw;
isTLE = strcmp(group_raw,'TLE');

for i = 1:numel(group_raw)
    if isTLE(i)
        if subjectTable.("lateralization_L")(i) == 1
            group_ref{i} = 'TLE Left';
        elseif subjectTable.("lateralization_R")(i) == 1
            group_ref{i} = 'TLE Right';
        elseif subjectTable.("lateralization_UC")(i) == 1
            group_ref{i} = 'TLE Unclear';
        else
            group_ref{i} = 'TLE unknown';
        end
    end
end

%% =========================
% Part 4: Define Site Groups
% =========================
siteNames = {'JLH','TJU'};

ColnamesBySite = struct();
ColnamesBySite.JLH = {'Focal Epilepsy','TLE','EXE','GGE','SeLECTS','AE','TLE Left','TLE Right','HP'};
ColnamesBySite.TJU = {'Focal Epilepsy','TLE','EXE','TLE Left','TLE Right','HP'};

groupTableBySite = struct();
for iSite = 1:numel(siteNames)
    sitename = siteNames{iSite};
    idx_site = strcmp(site, sitename);

    columns = ColnamesBySite.(sitename);
    Ns = sum(idx_site);
    Ng = numel(columns);

    M = zeros(Ns, Ng);
    g_raw = group_raw(idx_site);
    g_ref = group_ref(idx_site);

    for i = 1:Ns
        gr = g_raw{i};
        gf = g_ref{i};

        if strcmp(gr, 'HP'); M(i, strcmp(columns,'HP')) = 1; continue; end

        if strcmp(gr, 'TLE')
            if any(strcmp(columns,'TLE'));            M(i, strcmp(columns,'TLE')) = 1; end
            if any(strcmp(columns,'Focal Epilepsy')); M(i, strcmp(columns,'Focal Epilepsy')) = 1; end
            if strcmp(gf,'TLE Left') && any(strcmp(columns,'TLE Left'));     M(i, strcmp(columns,'TLE Left')) = 1; end
            if strcmp(gf,'TLE Right') && any(strcmp(columns,'TLE Right'));   M(i, strcmp(columns,'TLE Right')) = 1; end

        elseif strcmp(gr, 'EXE')
            if any(strcmp(columns,'EXE'));            M(i, strcmp(columns,'EXE')) = 1; end
            if any(strcmp(columns,'Focal Epilepsy')); M(i, strcmp(columns,'Focal Epilepsy')) = 1; end

        elseif strcmp(gr, 'GGE') && any(strcmp(columns,'GGE'));         M(i, strcmp(columns,'GGE')) = 1;
        elseif strcmp(gr, 'SeLECTS') && any(strcmp(columns,'SeLECTS')); M(i, strcmp(columns,'SeLECTS')) = 1;
        elseif strcmp(gr, 'AE') && any(strcmp(columns,'AE'));           M(i, strcmp(columns,'AE')) = 1;
        end
    end
    groupTableBySite.(sitename) = array2table(M, 'VariableNames', columns);
end

% Load 34-network assignment first
network_assignment = [];
if exist(assign34_mat,'file') == 2
    Smap = load(assign34_mat);
    if isfield(Smap,'network_assignment'); network_assignment = Smap.network_assignment; end
end

dataTypes = {'RawAtoms','DenoisedAtoms'};
ROI_defs = struct();
ROI_defs.AS200K17_200ROI = 1:200;

if ~isempty(network_assignment)
    ROI_defs.AS200K34 = network_assignment.network_order_name;
else
    ROI_defs.AS200K34 = 1:34; % Fallback
end

ROI_defs.Subcortical = { ...
    'Left_Accumbens'; ...
    'Left_Amygdala'; ...
    'Left_Caudate'; ...
    'Left_Hippocampus'; ...
    'Left_Pallidum'; ...
    'Left_Putamen'; ...
    'Left_Thalamus'; ...
    'Right_Accumbens'; ...
    'Right_Amygdala'; ...
    'Right_Caudate'; ...
    'Right_Hippocampus'; ...
    'Right_Pallidum'; ...
    'Right_Putamen'; ...
    'Right_Thalamus'};

roi_types = fieldnames(ROI_defs);

%% =========================================================
% PART 5: ROI extraction
% =========================================================
fprintf('\n=== PART 5: ROI extraction (hubness maps) ===\n');

MeasuresHubness = struct();
for iType = 1:numel(dataTypes)
    dataType = dataTypes{iType};
    for r = 1:numel(roi_types)
        roi = roi_types{r};
        MeasuresHubness.(dataType).(roi) = nan(Nsub, numel(ROI_defs.(roi)));
    end
end

AtomN = struct();
AtomN.RawAtoms = nan(Nsub,1);

for i = 1:Nsub
    subjID = subjectIDs{i};

    if isfield(All_derivatives_struct, subjID)
        Ssub = All_derivatives_struct.(subjID);
        ref_atlas = 'AS200Y17';
        AtomN.RawAtoms(i) = numel(Ssub.(ref_atlas).Atoms);
    end

    hub_root = SPARK_dir;
    fn_raw   = fullfile(hub_root, subjID, sprintf('%s_khubness.nii.gz', subjID));
    fn_clean = fullfile(hub_root, subjID, sprintf('%s_khubness_clean.nii.gz', subjID));

    if exist(fn_raw,'file') == 2
        MeasuresHubness.RawAtoms.AS200K17_200ROI(i,:) = mean_by_atlas(fn_raw, AtlasNii.AS200K17_200ROI, ROI_defs.AS200K17_200ROI);
        MeasuresHubness.RawAtoms.AS200K34(i,:)        = mean_by_atlas(fn_raw, AtlasNii.AS200K34,        ROI_defs.AS200K34);
        MeasuresHubness.RawAtoms.Subcortical(i,:)     = mean_by_atlas(fn_raw, AtlasNii.Subcortical,    ROI_defs.Subcortical);
    else
        warning('Missing raw hubness: %s', fn_raw);
    end

    if exist(fn_clean,'file') == 2
        MeasuresHubness.DenoisedAtoms.AS200K17_200ROI(i,:) = mean_by_atlas(fn_clean, AtlasNii.AS200K17_200ROI, ROI_defs.AS200K17_200ROI);
        MeasuresHubness.DenoisedAtoms.AS200K34(i,:)        = mean_by_atlas(fn_clean, AtlasNii.AS200K34,        ROI_defs.AS200K34);
        MeasuresHubness.DenoisedAtoms.Subcortical(i,:)     = mean_by_atlas(fn_clean, AtlasNii.Subcortical,    ROI_defs.Subcortical);
    else
        warning('Missing clean hubness: %s', fn_clean);
    end

    if mod(i,200)==0; fprintf('  Extracted %d / %d\n', i, Nsub); end
end

save(pathOutMeasures, 'MeasuresHubness', 'AtomN', '-v7.3');
fprintf('Saved ROI extraction: %s\n', pathOutMeasures);

%% =========================================================
% PART 6: Normative modelling (W-scores)
% =========================================================
fprintf('\n=== PART 6: Normative modelling (W-scores) ===\n');

% WscoreHubness = struct();

for iSite = 1:numel(siteNames)
    sitename = siteNames{iSite};
    fprintf('\n--- Site: %s ---\n', sitename);

    idx_site = strcmp(site, sitename);
    idx_use  = idx_site;

    idx_HP = idx_use & strcmp(group_ref,'HP');
    idx_HP_site = idx_HP(idx_use);
    fprintf('  HP baseline n = %d\n', sum(idx_HP_site));

    X_cov_all = [sex, age, hm, nan(Nsub,1)];

    for iType = 1:numel(dataTypes)
        dataType = dataTypes{iType};
        fprintf('  DataType: %s\n', dataType);

        X_cov_all(:,4) = double(AtomN.RawAtoms);
        X_site = X_cov_all(idx_use,:);

        for r = 1:numel(roi_types)
            roi = roi_types{r};
            Y_all  = MeasuresHubness.(dataType).(roi);
            Y_site = Y_all(idx_use,:);

            [W, Resid] = compute_wscore_pcn(X_site, Y_site, idx_HP_site);

            WscoreHubness.(sitename).(dataType).(roi).W     = W;
            WscoreHubness.(sitename).(dataType).(roi).Resid = Resid;
            WscoreHubness.(sitename).(dataType).(roi).FeatureNames = ROI_defs.(roi)(:);

            % Calculate HP raw stats
            WscoreHubness.(sitename).(dataType).(roi).HP_mean = mean(Y_site(idx_HP_site, :), 1, 'omitnan');
            WscoreHubness.(sitename).(dataType).(roi).HP_std  = std(Y_site(idx_HP_site, :), 0, 1, 'omitnan');
        end
    end
end

save(pathOutWscore, 'WscoreHubness', '-v7.3');
fprintf('Saved W-scores: %s\n', pathOutWscore);

%% =========================================================
% PART 7: Statistics + outputs
% =========================================================
fprintf('\n=== PART 7: Statistics + outputs ===\n');

alpha   = 0.05;
nperm   = 10000;
z_range = [-0.6 0.6];
color_range_hp = [1 4];
cmap_name = 'RdBu_r';
cmap_name_hp = 'Reds';

network_assignment = [];
if exist(assign34_mat,'file') == 2
    Smap = load(assign34_mat);
    if isfield(Smap,'network_assignment'); network_assignment = Smap.network_assignment; end
end

StatsHubness = struct();

for iSite = 1:numel(siteNames)
    sitename = siteNames{iSite};
    groupTable = groupTableBySite.(sitename);
    Group_names = groupTable.Properties.VariableNames;

    for iType = 1:numel(dataTypes)
        dataType = dataTypes{iType};

        for r = 1:numel(roi_types)
            roi = roi_types{r};
            W = WscoreHubness.(sitename).(dataType).(roi).W;
            P = size(W,2);

            out_dir_roi = fullfile(out_dir);
            if ~exist(out_dir_roi,'dir'); mkdir(out_dir_roi); end

            for g = 1:numel(Group_names)
                gname = Group_names{g};
                if strcmpi(gname,'HP'); continue; end

                idx_g = logical(groupTable.(gname));
                if sum(idx_g) < 3; continue; end

                G = W(idx_g,:);
                mu = mean(G, 1, 'omitnan');

                sigma = std(G, 0, 1);
                sigma(sigma == 0) = eps;

                [Z_perm, p_perm] = permuztest(G, 0, sigma, ...
                    'nperm', nperm, 'alpha', alpha, 'tail', 'both', 'correct', true);

                h_perm = p_perm < alpha;
                mu_tmax = mu; mu_tmax(~h_perm) = 0;

                % Save Stats
                gfield = matlab.lang.makeValidName(strrep(gname,' ','_'));
                StatsHubness.(sitename).(dataType).(roi).(gfield).n = sum(idx_g);
                StatsHubness.(sitename).(dataType).(roi).(gfield).mu = mu;
                StatsHubness.(sitename).(dataType).(roi).(gfield).perm.Z = Z_perm;
                StatsHubness.(sitename).(dataType).(roi).(gfield).perm.p = p_perm;
                StatsHubness.(sitename).(dataType).(roi).(gfield).perm.h = h_perm;

                % Output Maps
                suffix = sprintf('%s_%s_%s_%s', sitename, dataType, roi, strrep(gname,' ','_'));

                switch roi
                    case 'AS200K17_200ROI'
                        parcel_to_volume(mu,      AtlasNii.AS200K17_200ROI, fullfile(out_dir_roi, [suffix '_meanW.nii']),      [1:size(ROI_defs.AS200K17_200ROI,2)]);
                        parcel_to_volume(mu_tmax, AtlasNii.AS200K17_200ROI, fullfile(out_dir_roi, [suffix '_meanW_Tmax.nii']), [1:size(ROI_defs.AS200K17_200ROI,2)]);

                        % Surface
                        pv = parcel_to_surface(mu, 'schaefer_200x17_conte69');
                        f = figure('Color','w','Position',[100 100 960 720]);
                        plot_cortical(pv, 'surface_name','conte69', 'color_range', z_range, 'cmap', cmap_name);
                        print(f, '-dtiff', '-r300', fullfile(out_dir_roi, [suffix '_meanW_conte69.tiff']));
                        close(f);

                        pv = parcel_to_surface(mu_tmax, 'schaefer_200x17_conte69');
                        f = figure('Color','w','Position',[100 100 960 720]);
                        plot_cortical(pv, 'surface_name','conte69', 'color_range', z_range, 'cmap', cmap_name);
                        print(f, '-dtiff', '-r300', fullfile(out_dir_roi, [suffix '_meanW_Tmax_conte69.tiff']));
                        close(f);

                    case 'Subcortical'
                        parcel_to_volume(mu,      AtlasNii.Subcortical, fullfile(out_dir_roi, [suffix '_meanW.nii']),      [1:size(ROI_defs.Subcortical,1)]);
                        parcel_to_volume(mu_tmax, AtlasNii.Subcortical, fullfile(out_dir_roi, [suffix '_meanW_Tmax.nii']), [1:size(ROI_defs.Subcortical,1)]);

                        parcel_to_Subsurface(mu,      fullfile(out_dir_roi, [suffix '_meanW.tiff']),      z_range, cmap_name);
                        parcel_to_Subsurface(mu_tmax, fullfile(out_dir_roi, [suffix '_meanW_Tmax.tiff']), z_range, cmap_name);

                    case 'AS200K34'
                        if ~isempty(network_assignment)
                            map = network_assignment.mapping;
                            mu200      = mu(map);
                            mu200_tmax = mu_tmax(map);

                            parcel_to_volume(mu200,      AtlasNii.AS200K17_200ROI, fullfile(out_dir_roi, [suffix '_meanW_on200.nii']), [1:size(ROI_defs.AS200K17_200ROI,2)]);
                            parcel_to_volume(mu200_tmax, AtlasNii.AS200K17_200ROI, fullfile(out_dir_roi, [suffix '_meanW_Tmax_on200.nii']), [1:size(ROI_defs.AS200K17_200ROI,2)]);

                            pv = parcel_to_surface(mu200, 'schaefer_200x17_conte69');
                            f = figure('Color','w','Position',[100 100 960 720]);
                            plot_cortical(pv, 'surface_name','conte69', 'color_range', z_range, 'cmap', cmap_name);
                            print(f, '-dtiff', '-r300', fullfile(out_dir_roi, [suffix '_meanW_on200_conte69.tiff']));
                            close(f);

                            pv = parcel_to_surface(mu200_tmax, 'schaefer_200x17_conte69');
                            f = figure('Color','w','Position',[100 100 960 720]);
                            plot_cortical(pv, 'surface_name','conte69', 'color_range', z_range, 'cmap', cmap_name);
                            print(f, '-dtiff', '-r300', fullfile(out_dir_roi, [suffix '_meanW_Tmax_on200_conte69.tiff']));
                            close(f);
                        end
                end

                % --- HP raw mean plot ---
                if isfield(WscoreHubness.(sitename).(dataType).(roi), 'HP_mean')
                    mu_hp = WscoreHubness.(sitename).(dataType).(roi).HP_mean;
                    suffix_hp = sprintf('%s_%s_%s_HP', sitename, dataType, roi);

                    switch roi
                        case 'AS200K17_200ROI'
                            pv = parcel_to_surface(mu_hp, 'schaefer_200x17_conte69');
                            f = figure('Color','w','Position',[100 100 960 720]);
                            plot_cortical(pv, 'surface_name','conte69', 'color_range', color_range_hp, 'cmap', cmap_name_hp);
                            print(f, '-dtiff', '-r300', fullfile(out_dir_roi, [suffix_hp '_meanW_conte69.tiff']));
                            close(f);

                        case 'Subcortical'
                            parcel_to_Subsurface(mu_hp, fullfile(out_dir_roi, [suffix_hp '_meanW.tiff']), color_range_hp, cmap_name_hp);

                        case 'AS200K34'
                            if ~isempty(network_assignment)
                                map = network_assignment.mapping;
                                mu200_hp = mu_hp(map);
                                pv = parcel_to_surface(mu200_hp, 'schaefer_200x17_conte69');
                                f = figure('Color','w','Position',[100 100 960 720]);
                                plot_cortical(pv, 'surface_name','conte69', 'color_range', color_range_hp, 'cmap', cmap_name_hp);
                                print(f, '-dtiff', '-r300', fullfile(out_dir_roi, [suffix_hp '_meanW_on200_conte69.tiff']));
                                close(f);
                            end
                    end
                end
            end

            % --- Group Boxplots for specific ROIs ---
            if ismember(roi, {'AS200K34', 'Subcortical'})
                out_dir_box = fullfile(out_dir_roi, 'Group_Boxplots');
                if ~exist(out_dir_box, 'dir'); mkdir(out_dir_box); end

                fprintf('    Generating boxplots for %s...\n', roi);

                for iR = 1:size(W, 2)
                    w_vec = W(:, iR);
                    % Use FeatureNames if available, otherwise generic
                    if isfield(WscoreHubness.(sitename).(dataType).(roi), 'FeatureNames')
                        fNames = WscoreHubness.(sitename).(dataType).(roi).FeatureNames;
                        if iR <= numel(fNames)
                            featName = char(fNames{iR}); % Ensure char
                        else
                            featName = sprintf('ROI%03d', iR);
                        end
                    else
                        featName = sprintf('ROI%03d', iR);
                    end

                    % Sanitize filename
                    featNameSafe = matlab.lang.makeValidName(featName);

                    out_name = fullfile(out_dir_box, sprintf('%s_%s_%s_%s_box.tiff', sitename, dataType, roi, featNameSafe));
                    label_str = sprintf('%s | %s | %s | %s', sitename, dataType, roi, featName);

                    % Reuse plotting function from S4 (plot_group_wscore)
                    % Function signature: [w_median, p_val, labs, w_std] = plot_group_wscore(data, groupTable, title_str, out_path)
                    plot_group_wscore(w_vec, groupTable, label_str, out_name);
                    close all;
                end
            end
        end
    end
end

save(pathOutStats, 'StatsHubness', '-v7.3');
fprintf('\nDONE.\n  Raw ROI: %s\n  Wscore:   %s\n  Stats:    %s\n  Figures:  %s\n', ...
    pathOutMeasures, pathOutWscore, pathOutStats, out_dir);


%% =========================================================
% PART 9: Correspondence Analysis (TJU RawAtoms only)
% =========================================================
fprintf('\n=== PART 9: Correspondence Analysis (TJU RawAtoms) ===\n');

% Only run if TJU RawAtoms
if strcmp(sitename, 'TJU') && strcmp(dataType, 'RawAtoms')

    out_dir_corr = fullfile(out_dir, 'Correspondence_Analysis');
    if ~exist(out_dir_corr, 'dir'); mkdir(out_dir_corr); end

    % Ensure required variables are defined (e.g. if running this section standalone)
    if ~exist('alpha', 'var'); alpha = 0.05; end
    if ~exist('nperm', 'var'); nperm = 10000; end
    if ~exist('z_range_r', 'var'); z_range_r = [-0.5 0.5]; end
    if ~exist('z_range_t', 'var'); z_range_t = [-5 5]; end

    csv_file = fullfile(project_root, 'outputs', 'PySuStaln', 'Correspondence', 'K02', 'assignment_K2.csv');

    if exist(csv_file, 'file')
        T_K2 = readtable(csv_file);

        % We need to match subjects in 'pat_indices' (TJU non-HP) with T_K2.SubID
        % Note: pat_indices uses the global 'subjectTable' order
        idx_site = strcmp(site, 'TJU');
        idx_hp = strcmp(group_raw, 'HP');
        idx_pat = idx_site & ~idx_hp;
        pat_indices = find(idx_pat);
        Npat = numel(pat_indices);
        idx_pat_in_W = ~strcmp(group_raw(idx_site), 'HP');

        % Initialize vectors for analysis
        % Subtype: 0 or 1. Stage: 0-N
        SubtypeVec = nan(Npat, 1);
        StageVec   = nan(Npat, 1);

        valid_sub_mask = false(Npat, 1);

        for p = 1:Npat
            % Global index
            g_idx = pat_indices(p);
            subID = subjectIDs{g_idx};

            % Find in K2 table
            row_k2 = strcmp(T_K2.SubID, subID);
            if any(row_k2)
                SubtypeVec(p) = T_K2.ml_subtype(row_k2);
                StageVec(p)   = T_K2.ml_stage(row_k2);
                valid_sub_mask(p) = true;
            end
        end

        fprintf('  Matched %d subjects from assignment_K2.csv\n', sum(valid_sub_mask));

        % --- Analysis Loop ---
        corr_rois = {'AS200K17_200ROI','AS200K34', 'Subcortical'};
        corr_rois = {'AS200K34'};
        for r = 1:numel(corr_rois)
            roi = corr_rois{r};
            if ~isfield(WscoreHubness.(sitename).(dataType), roi); continue; end

            % Extract W for patients using the site-specific index
            W_pat = WscoreHubness.(sitename).(dataType).(roi).W(idx_pat_in_W, :);

            % 1. Correspondence Stage Correlation
            valid = valid_sub_mask & ~isnan(StageVec) & ~all(isnan(W_pat),2);
            if sum(valid) > 5
                [r_corr, p_corr] = permucorr(W_pat(valid,:), StageVec(valid), ...
                    'nperm', nperm, 'alpha', alpha, 'tail', 'both', 'type','spearman','correct', true);

                suffix = sprintf('%s_%s_%s_Corr_Stage', sitename, dataType, roi);

                % Inline Plotting
                stat_val = r_corr; p_val = p_corr; z_lim = z_range_r;
                stat_sig = stat_val; stat_sig(p_val >= alpha) = 0;

                if strcmp(roi, 'Subcortical')
                    parcel_to_Subsurface(stat_val, fullfile(out_dir_corr, [suffix '_raw.tiff']), z_lim, 'RdBu_r');
                    parcel_to_Subsurface(stat_sig, fullfile(out_dir_corr, [suffix '_sig.tiff']), z_lim, 'RdBu_r');
                elseif strcmp(roi, 'AS200K17_200ROI')
                    % Direct mapping for 200 ROI
                    mu200 = stat_val;
                    pv = parcel_to_surface(mu200, 'schaefer_200x17_conte69');
                    f = figure('Color','w','Position',[100 100 960 720]);
                    plot_cortical(pv, 'surface_name','conte69', 'color_range', z_lim, 'cmap', 'RdBu_r');
                    print(f, '-dtiff', '-r300', fullfile(out_dir_corr, [suffix '_raw_conte69.tiff']));
                    close(f);
                    
                    mu200_sig = stat_sig;
                    pv = parcel_to_surface(mu200_sig, 'schaefer_200x17_conte69');
                    f = figure('Color','w','Position',[100 100 960 720]);
                    plot_cortical(pv, 'surface_name','conte69', 'color_range', z_lim, 'cmap', 'RdBu_r');
                    print(f, '-dtiff', '-r300', fullfile(out_dir_corr, [suffix '_sig_conte69.tiff']));
                    close(f);

                elseif strcmp(roi, 'AS200K34') && ~isempty(network_assignment)
                    map = network_assignment.mapping;

                    mu200 = stat_val(map);
                    pv = parcel_to_surface(mu200, 'schaefer_200x17_conte69');
                    f = figure('Color','w','Position',[100 100 960 720]);
                    plot_cortical(pv, 'surface_name','conte69', 'color_range', z_lim, 'cmap', 'RdBu_r');
                    print(f, '-dtiff', '-r300', fullfile(out_dir_corr, [suffix '_raw_conte69.tiff']));
                    close(f);

                    mu200_sig = stat_sig(map);
                    pv = parcel_to_surface(mu200_sig, 'schaefer_200x17_conte69');
                    f = figure('Color','w','Position',[100 100 960 720]);
                    plot_cortical(pv, 'surface_name','conte69', 'color_range', z_lim, 'cmap', 'RdBu_r');
                    print(f, '-dtiff', '-r300', fullfile(out_dir_corr, [suffix '_sig_conte69.tiff']));
                    close(f);
                    
                    % BAR PLOTS (S11 Order, Compound L/R)
                    local_plot_bar(stat_val(1:17),  p_val(1:17), ...
                                       stat_val(18:34), p_val(18:34), ...
                                       [suffix '_bar_L_R'], out_dir_corr);
                end
            end

            % 2. Correspondence Subtype T-test
            % Exclude Stage 0
            valid_subtype = valid_sub_mask & ~isnan(SubtypeVec) & (StageVec ~= 0) & ~all(isnan(W_pat),2);

            % Assuming Subtype is 0 and 1
            if sum(valid_subtype) > 5
                g1 = W_pat(valid_subtype & SubtypeVec == 0, :); % Subtype 0
                g2 = W_pat(valid_subtype & SubtypeVec == 1, :); % Subtype 1

                if size(g1,1) > 2 && size(g2,1) > 2
                    [t_val, p_val] = permuttest2(g1, g2, ...
                        'nperm', nperm, 'alpha', alpha, 'tail', 'both', 'correct', true);

                    suffix = sprintf('%s_%s_%s_Ttest_Subtype0vs1', sitename, dataType, roi);

                    % Inline Plotting
                    stat_val = t_val; p_val = p_val; z_lim = z_range_t;
                    stat_sig = stat_val; stat_sig(p_val >= alpha) = 0;

                    if strcmp(roi, 'Subcortical')
                        parcel_to_Subsurface(stat_val, fullfile(out_dir_corr, [suffix '_raw.tiff']), z_lim, 'RdBu_r');
                        parcel_to_Subsurface(stat_sig, fullfile(out_dir_corr, [suffix '_sig.tiff']), z_lim, 'RdBu_r');
                    elseif strcmp(roi, 'AS200K17_200ROI')
                        % Direct mapping for 200 ROI
                        mu200 = stat_val;
                        pv = parcel_to_surface(mu200, 'schaefer_200x17_conte69');
                        f = figure('Color','w','Position',[100 100 960 720]);
                        plot_cortical(pv, 'surface_name','conte69', 'color_range', z_lim, 'cmap', 'RdBu_r');
                        print(f, '-dtiff', '-r300', fullfile(out_dir_corr, [suffix '_raw_conte69.tiff']));
                        close(f);
                        
                        mu200_sig = stat_sig;
                        pv = parcel_to_surface(mu200_sig, 'schaefer_200x17_conte69');
                        f = figure('Color','w','Position',[100 100 960 720]);
                        plot_cortical(pv, 'surface_name','conte69', 'color_range', z_lim, 'cmap', 'RdBu_r');
                        print(f, '-dtiff', '-r300', fullfile(out_dir_corr, [suffix '_sig_conte69.tiff']));
                        close(f);

                    elseif strcmp(roi, 'AS200K34') && ~isempty(network_assignment)
                        map = network_assignment.mapping;

                        mu200 = stat_val(map);
                        pv = parcel_to_surface(mu200, 'schaefer_200x17_conte69');
                        f = figure('Color','w','Position',[100 100 960 720]);
                        plot_cortical(pv, 'surface_name','conte69', 'color_range', z_lim, 'cmap', 'RdBu_r');
                        print(f, '-dtiff', '-r300', fullfile(out_dir_corr, [suffix '_raw_conte69.tiff']));
                        close(f);

                        mu200_sig = stat_sig(map);
                        pv = parcel_to_surface(mu200_sig, 'schaefer_200x17_conte69');
                        f = figure('Color','w','Position',[100 100 960 720]);
                        plot_cortical(pv, 'surface_name','conte69', 'color_range', z_lim, 'cmap', 'RdBu_r');
                        print(f, '-dtiff', '-r300', fullfile(out_dir_corr, [suffix '_sig_conte69.tiff']));
                        close(f);
                        
                        % BAR PLOTS (S11 Order, Compound L/R)
                        local_plot_bar(stat_val(1:17),  p_val(1:17), ...
                                           stat_val(18:34), p_val(18:34), ...
                                           [suffix '_bar_L_R'], out_dir_corr);
                    end
            end
        end
        end

        % 3. Correspondence Stage Correlation (Subtype 0 only)
        valid_s0 = valid_sub_mask & (SubtypeVec == 0) & ~isnan(StageVec) & ~all(isnan(W_pat),2);
        if sum(valid_s0) > 5
            [r_corr, p_corr] = permucorr(W_pat(valid_s0,:), StageVec(valid_s0), ...
                'nperm', nperm, 'alpha', alpha, 'tail', 'both', 'type', 'spearman', 'correct', true);

            suffix = sprintf('%s_%s_%s_Corr_Stage_Subtype0', sitename, dataType, roi);

            % Inline Plotting
            stat_val = r_corr; p_val = p_corr; z_lim = z_range_r;
            stat_sig = stat_val; stat_sig(p_val >= alpha) = 0;

            if strcmp(roi, 'Subcortical')
                parcel_to_Subsurface(stat_val, fullfile(out_dir_corr, [suffix '_raw.tiff']), z_lim, 'RdBu_r');
                parcel_to_Subsurface(stat_sig, fullfile(out_dir_corr, [suffix '_sig.tiff']), z_lim, 'RdBu_r');
            elseif strcmp(roi, 'AS200K17_200ROI')
                % Direct mapping for 200 ROI
                mu200 = stat_val;
                pv = parcel_to_surface(mu200, 'schaefer_200x17_conte69');
                f = figure('Color','w','Position',[100 100 960 720]);
                plot_cortical(pv, 'surface_name','conte69', 'color_range', z_lim, 'cmap', 'RdBu_r');
                print(f, '-dtiff', '-r300', fullfile(out_dir_corr, [suffix '_raw_conte69.tiff']));
                close(f);
                
                mu200_sig = stat_sig;
                pv = parcel_to_surface(mu200_sig, 'schaefer_200x17_conte69');
                f = figure('Color','w','Position',[100 100 960 720]);
                plot_cortical(pv, 'surface_name','conte69', 'color_range', z_lim, 'cmap', 'RdBu_r');
                print(f, '-dtiff', '-r300', fullfile(out_dir_corr, [suffix '_sig_conte69.tiff']));
                close(f);

            elseif strcmp(roi, 'AS200K34') && ~isempty(network_assignment)
                map = network_assignment.mapping;

                mu200 = stat_val(map);
                pv = parcel_to_surface(mu200, 'schaefer_200x17_conte69');
                f = figure('Color','w','Position',[100 100 960 720]);
                plot_cortical(pv, 'surface_name','conte69', 'color_range', z_lim, 'cmap', 'RdBu_r');
                print(f, '-dtiff', '-r300', fullfile(out_dir_corr, [suffix '_raw_conte69.tiff']));
                close(f);

                mu200_sig = stat_sig(map);
                pv = parcel_to_surface(mu200_sig, 'schaefer_200x17_conte69');
                f = figure('Color','w','Position',[100 100 960 720]);
                plot_cortical(pv, 'surface_name','conte69', 'color_range', z_lim, 'cmap', 'RdBu_r');
                print(f, '-dtiff', '-r300', fullfile(out_dir_corr, [suffix '_sig_conte69.tiff']));
                close(f);
                
                % BAR PLOTS (S11 Order, Compound L/R)
                local_plot_bar(stat_val(1:17),  p_val(1:17), ...
                                   stat_val(18:34), p_val(18:34), ...
                                   [suffix '_bar_L_R'], out_dir_corr);
            end
        end
    end
else
    warning('assignment_K2.csv not found: %s', csv_file);
end



