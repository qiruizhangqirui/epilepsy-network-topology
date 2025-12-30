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

dataTypes = {'RawAtoms','DenoisedAtoms'};
ROI_defs = struct();
ROI_defs.AS200K17_200ROI = 1:200;
ROI_defs.AS200K34        = 1:34;
ROI_defs.Subcortical     = 1:14;
roi_types = fieldnames(ROI_defs);

%% =========================================================
% PART 1: ROI extraction
% =========================================================
fprintf('\n=== PART 1: ROI extraction (hubness maps) ===\n');

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
% PART 2: Normative modelling (W-scores)
% =========================================================
fprintf('\n=== PART 2: Normative modelling (W-scores) ===\n');

WscoreHubness = struct();

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
        end
    end
end

save(pathOutWscore, 'WscoreHubness', '-v7.3');
fprintf('Saved W-scores: %s\n', pathOutWscore);

%% =========================================================
% PART 3: Statistics + outputs
% =========================================================
fprintf('\n=== PART 3: Statistics + outputs ===\n');

alpha   = 0.05;
nperm   = 10000;
z_range = [-0.6 0.6];
cmap_name = 'RdBu_r';

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
                        parcel_to_volume(mu,      AtlasNii.AS200K17_200ROI, fullfile(out_dir_roi, [suffix '_meanW.nii']),      ROI_defs.AS200K17_200ROI);
                        parcel_to_volume(mu_tmax, AtlasNii.AS200K17_200ROI, fullfile(out_dir_roi, [suffix '_meanW_Tmax.nii']), ROI_defs.AS200K17_200ROI);

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
                        parcel_to_volume(mu,      AtlasNii.Subcortical, fullfile(out_dir_roi, [suffix '_meanW.nii']),      ROI_defs.Subcortical);
                        parcel_to_volume(mu_tmax, AtlasNii.Subcortical, fullfile(out_dir_roi, [suffix '_meanW_Tmax.nii']), ROI_defs.Subcortical);

                        parcel_to_Subsurface(mu,      fullfile(out_dir_roi, [suffix '_meanW.tiff']),      z_range, cmap_name);
                        parcel_to_Subsurface(mu_tmax, fullfile(out_dir_roi, [suffix '_meanW_Tmax.tiff']), z_range, cmap_name);

                    case 'AS200K34'
                        if ~isempty(network_assignment)
                            map = network_assignment.mapping;
                            mu200      = mu(map);
                            mu200_tmax = mu_tmax(map);

                            parcel_to_volume(mu200,      AtlasNii.AS200K17_200ROI, fullfile(out_dir_roi, [suffix '_meanW_on200.nii']), ROI_defs.AS200K17_200ROI);
                            parcel_to_volume(mu200_tmax, AtlasNii.AS200K17_200ROI, fullfile(out_dir_roi, [suffix '_meanW_Tmax_on200.nii']), ROI_defs.AS200K17_200ROI);

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
            end
        end
    end
end

save(pathOutStats, 'StatsHubness', '-v7.3');
fprintf('\nDONE.\n  Raw ROI: %s\n  Wscore:   %s\n  Stats:    %s\n  Figures:  %s\n', ...
    pathOutMeasures, pathOutWscore, pathOutStats, out_dir);
