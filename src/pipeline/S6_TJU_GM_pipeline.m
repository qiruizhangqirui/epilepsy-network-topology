% =========================================================================
% Script: S6_TJU_GM_pipeline.m
% -------------------------------------------------------------------------
% Description:
%   Pipeline for Gray Matter (GM) analysis on TJU dataset using CAT12 derivatives.
%
%   1. ROI extraction from CAT12 ROI files.
%   2. Normative modeling (PCNtoolkit) for W-scores.
%      - EXCLUDES subjects with QCfail > 0 during modeling.
%      - Subjects with QCfail > 0 are filled with NaN in outputs.
%   3. Max-T permutation testing and statistical mapping.
%
% Note:
%   - Requires 'permuztest'.
%   - Requires PCNtoolkit setup in Python.
%
% Inputs:
%   - data/CAT_Derivatives/CAT12.9_2577/
%   - data/Subjects.xlsx
%   - data/assignment_34.mat
%   - data/tian2aseg.mat
%
% Outputs:
%   - outputs/TJU_GM_raw.mat
%   - outputs/TJU_GM_Wscore.mat
%   - outputs/TJU_GM_stats.mat
%   - outputs/TJU_GM_statistic/
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


pathSubjects = fullfile(data_dir, 'Subjects.xlsx');
PathCAT12    = fullfile(data_dir, 'CAT_Derivatives', 'CAT12.9_2577');
pathAssign34 = fullfile(data_dir, 'assignment_34.mat');
pathTian2Aseg= fullfile(data_dir, 'tian2aseg.mat');

out_dir = fullfile(result_dir, 'TJU_GM_statistic');
if ~exist(out_dir,'dir'); mkdir(out_dir); end

pathOutMeasures = fullfile(result_dir, 'TJU_GM_Measures.mat');
pathOutWscore   = fullfile(result_dir, 'TJU_GM_Wscore.mat');
pathOutStats    = fullfile(result_dir, 'TJU_GM_stats.mat');

% Atlases for plotting (Using paths similar to S5, adjusted if necessary)
% Note: Assuming standard atlas location relative to data_dir
AtlasNii = struct();
AtlasNii.AS200K17_200ROI = fullfile(data_dir, 'NCT_atlases', 'atlases', 'FSLMNI2mm', 'YeoLab', 'AS200Y17.nii.gz');
AtlasNii.AS200K34        = fullfile(data_dir, 'AS17network_34ROIs_sorted.nii.gz');
AtlasNii.Subcortical     = fullfile(data_dir, 'Subcortical_FLS2mm.nii');

%% =========================
% Part 2: Load Metadata
% =========================
subjectTable = readtable(pathSubjects, 'VariableNamingRule','preserve');

subjectIDs = cellstr(string(subjectTable.("ID_new")));
site       = cellstr(string(subjectTable.("site")));
group_raw  = cellstr(string(subjectTable.("Group")));
QCfail = double(subjectTable.("T1w_failedQC"));

sex = double(subjectTable.("sex_b(M1F0)"));
age = double(subjectTable.("age"));
hm  = double(subjectTable.("mean_FD"));

% --- Filter: Keep only TJU subjects ---
target_site = 'TJU';
idx_keep = strcmp(site, target_site);

subjectIDs = subjectIDs(idx_keep);
site       = site(idx_keep);
group_raw  = group_raw(idx_keep);
QCfail     = QCfail(idx_keep);
sex        = sex(idx_keep);
age        = age(idx_keep);
hm         = hm(idx_keep);

fprintf('Filtered to TJU site only. N=%d\n', numel(subjectIDs));

TIV = nan(numel(subjectIDs), 1); % Will be loaded from CAT12

% ---- refine TLE lateralization ----
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

% Site Logic: TJU Only (Already filtered)
idx_site_all = true(numel(subjectIDs), 1);

%% =========================
% Part 3: Load Data (CAT12)
% =========================
fprintf('\n=== PART 1: Data Loading (CAT12) ===\n');

load(pathAssign34, 'network_assignment');
load(pathTian2Aseg, 'tian2aseg');

MeasuresGM = struct();
% dataTypes removed

% Pre-allocate
Nsub = numel(subjectIDs);
ROI_data_200ROI_All   = nan(Nsub, 200);
ROI_data_34ROI_All    = nan(Nsub, 34);
ROI_data_20ROI_All    = nan(Nsub, 20);
ROI_data_Subcortical_All = nan(Nsub, 14);

for i = 1:Nsub
    % Run for all subjects, filter later or use NaN for missing
    subjID = subjectIDs{i};

    GM_mat = fullfile(PathCAT12, 'label', ['catROI_', subjID, '.mat']);
    CAT_mat = fullfile(PathCAT12, 'report', ['cat_', subjID, '.mat']);

    if exist(GM_mat, 'file') == 2 && exist(CAT_mat, 'file') == 2
        try
            GM_data = load(GM_mat);
            % 200 ROI
            ROI_data_200ROI_All(i,:) = GM_data.S.Schaefer2018_200Parcels_17Networks_order.data.Vgm;

            % Subcortical
            Subcortical_data = GM_data.S.Tian_Subcortex_S4_7T.data.Vgm;
            % Subcortical_name = GM_data.S.Tian_Subcortex_S4_7T.names;
            ROI_data_Subcortical_All(i,:) = accumarray(tian2aseg, Subcortical_data, [14 1], @sum, 0);

            % 34 ROI
            ROI_data_34ROI_All(i,:) = accumarray(network_assignment.mapping, ROI_data_200ROI_All(i,:)', [34 1], @sum, 0);

            % 20 ROI (if mapping exists)
            if isfield(network_assignment, 'mapping_20')
                ROI_data_20ROI_All(i,:) = accumarray(network_assignment.mapping_20, ROI_data_200ROI_All(i,:)', [20 1], @sum, 0);
            end

            % TIV / IQR
            CAT_data = load(CAT_mat);
            TIV(i,1) = CAT_data.S.subjectmeasures.vol_TIV;
            % IQR(i,1) = CAT_data.S.qualityratings.IQR;
        catch ME
            warning('Error loading data for %s: %s', subjID, ME.message);
        end
    else
        % File missing
        % warning('Missing CAT12 files for %s', subjID);
    end

    if mod(i, 100) == 0; fprintf('  Loaded %d / %d\n', i, Nsub); end
end

% Store in struct (Flattened)
MeasuresGM.AS200K17_200ROI = ROI_data_200ROI_All;
MeasuresGM.AS200K34        = ROI_data_34ROI_All;
MeasuresGM.AS200K20        = ROI_data_20ROI_All;
MeasuresGM.Subcortical     = ROI_data_Subcortical_All;

ROI_defs = struct();
ROI_defs.AS200K17_200ROI = 1:200;
ROI_defs.AS200K34        = network_assignment.network_order_name;
ROI_defs.AS200K20        = network_assignment.network_order_name20;
ROI_defs.Subcortical     = 1:14;
roi_types = fieldnames(ROI_defs);

save(pathOutMeasures, 'MeasuresGM', 'TIV', '-v7.3');
fprintf('Saved raw data: %s\n', pathOutMeasures);


%% =========================
% Part 4: Normative Modeling (W-scores)
% =========================
fprintf('\n=== PART 2: Normative modelling (W-scores) ===\n');

WscoreGM = struct();
sitename = target_site;

% Site subset
idx_site = idx_site_all;

% 1. Identify valid subjects for modeling: Site match AND QC passed
idx_valid = idx_site & (QCfail == 0);

% 2. Identify HP baseline among valid subjects
idx_HP_valid = idx_valid & strcmp(group_ref, 'HP');

fprintf('Site: %s\n', sitename);
fprintf('  Total N: %d\n', sum(idx_site));
fprintf('  Valid N (QC pass): %d\n', sum(idx_valid));
fprintf('  HP Baseline N: %d\n', sum(idx_HP_valid));

% Prepare Covariates (Nsub x 4)
% [Sex, Age, MeanFD, TIV]
X_cov_all = [sex, age, hm, TIV];

for r = 1:numel(roi_types)
    roi = roi_types{r};
    Y_all = MeasuresGM.(roi);
    % --- Logic for exclusion ---
    % We only use 'idx_valid' subjects to train/predict.
    % We populate the full Nsub array with NaN for excluded subjects.

    % Subset to valid subjects
    X_in = X_cov_all(idx_valid, :);
    Y_in = Y_all(idx_valid, :);

    % Logical index of HP within the valid subset
    % idx_HP_valid corresponds to the full array.
    % We need boolean vector matching rows of X_in/Y_in.
    idx_HP_in = strcmp(group_ref(idx_valid), 'HP');

    % Compute W-score
    [W_valid, Resid_valid] = compute_wscore_pcn(X_in, Y_in, idx_HP_in);

    % Map back to full size (Nsub x Features)
    nFeats = size(Y_all, 2);
    W_full = nan(Nsub, nFeats);
    R_full = nan(Nsub, nFeats);

    W_full(idx_valid, :) = W_valid;
    R_full(idx_valid, :) = Resid_valid;

    % Store - using Sitename.ROI structure
    WscoreGM.(sitename).(roi).W     = W_full;
    WscoreGM.(sitename).(roi).Resid = R_full;
    WscoreGM.(sitename).(roi).FeatureNames = ROI_defs.(roi)(:);

    WscoreGM.(sitename).(roi).idx_valid = idx_valid; % Store mask for reference
end

save(pathOutWscore, 'WscoreGM', '-v7.3');
fprintf('Saved W-scores: %s\n', pathOutWscore);


%% =========================
% Part 5: Statistics
% =========================
fprintf('\n=== PART 3: Statistics + Outputs ===\n');

alpha   = 0.05;
nperm   = 10000;
z_range = [-0.6 0.6];
cmap_name = 'RdBu_r';

% Define TJU Groups manually if not carried over from S5 logic
ColnamesTJU = {'Focal Epilepsy','TLE','EXE','TLE Left','TLE Right'};
% Build Group Table for TJU
M = zeros(sum(idx_site), numel(ColnamesTJU));
g_raw_site = group_raw(idx_site);
g_ref_site = group_ref(idx_site);

for i = 1:numel(g_raw_site)
    gr = g_raw_site{i};
    gf = g_ref_site{i};

    if strcmp(gr, 'TLE')
        if any(strcmp(ColnamesTJU,'TLE')); M(i, strcmp(ColnamesTJU,'TLE')) = 1; end
        if any(strcmp(ColnamesTJU,'Focal Epilepsy')); M(i, strcmp(ColnamesTJU,'Focal Epilepsy')) = 1; end
        if strcmp(gf,'TLE Left') && any(strcmp(ColnamesTJU,'TLE Left')); M(i, strcmp(ColnamesTJU,'TLE Left')) = 1; end
        if strcmp(gf,'TLE Right') && any(strcmp(ColnamesTJU,'TLE Right')); M(i, strcmp(ColnamesTJU,'TLE Right')) = 1; end
    elseif strcmp(gr, 'EXE')
        if any(strcmp(ColnamesTJU,'EXE')); M(i, strcmp(ColnamesTJU,'EXE')) = 1; end
        if any(strcmp(ColnamesTJU,'Focal Epilepsy')); M(i, strcmp(ColnamesTJU,'Focal Epilepsy')) = 1; end
    end
end
T_Group = array2table(M, 'VariableNames', ColnamesTJU);

StatsGM = struct();

for r = 1:numel(roi_types)
    roi = roi_types{r};
    if ~isfield(WscoreGM.(sitename), roi); continue; end

    W_full = WscoreGM.(sitename).(roi).W;

    out_dir_roi = fullfile(out_dir);

    for g = 1:numel(ColnamesTJU)
        gname = ColnamesTJU{g};

        % Get group indices relative to the SITE subset
        % idx_site selects the site rows from the total subject list.
        % T_Group corresponds to these rows.
        idx_g_in_site = logical(T_Group.(gname));

        % Map back to full Nsub array
        idx_g_full = false(Nsub, 1);
        site_indices = find(idx_site);
        idx_g_full(site_indices(idx_g_in_site)) = true;

        % Intersect with Valid QC subjects
        % (W_full already has NaNs for invalid subjects, but strictly we filter)
        idx_use_stats = idx_g_full & (QCfail == 0);

        if sum(idx_use_stats) < 3; continue; end

        G = W_full(idx_use_stats, :);

        % Double check for NaNs (should be none if logic holds)
        if any(isnan(G(:)))
            warning('NaNs found in valid group data. Removing rows with NaNs.');
            G = G(all(~isnan(G),2), :);
        end
        if size(G,1) < 3; continue; end

        mu = mean(G, 1);
        sigma = std(G, 0, 1);
        sigma(sigma == 0) = eps;

        [Z_perm, p_perm] = permuztest(G, 0, sigma, ...
            'nperm', nperm, 'alpha', alpha, 'tail', 'both', 'correct', true);

        h_perm = p_perm < alpha;
        mu_tmax = mu;
        mu_tmax(~h_perm) = 0;

        % Save Stats
        gfield = matlab.lang.makeValidName(strrep(gname,' ','_'));
        StatsGM.(sitename).(roi).(gfield).n = size(G,1);
        StatsGM.(sitename).(roi).(gfield).mu = mu;
        StatsGM.(sitename).(roi).(gfield).perm.Z = Z_perm;
        StatsGM.(sitename).(roi).(gfield).perm.p = p_perm;
        StatsGM.(sitename).(roi).(gfield).perm.h = h_perm;

        % Output Maps
        suffix = sprintf('%s_%s_%s', sitename, roi, strrep(gname,' ','_'));

        % Plotting logic dependent on ROI type
        switch roi
            case 'AS200K17_200ROI'
                % Volume
                if exist(AtlasNii.AS200K17_200ROI, 'file')
                    parcel_to_volume(mu,      AtlasNii.AS200K17_200ROI, fullfile(out_dir_roi, [suffix '_meanW.nii']),      ROI_defs.(roi));
                    parcel_to_volume(mu_tmax, AtlasNii.AS200K17_200ROI, fullfile(out_dir_roi, [suffix '_meanW_Tmax.nii']), ROI_defs.(roi));
                end
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
                if exist(AtlasNii.Subcortical, 'file')
                    parcel_to_volume(mu,      AtlasNii.Subcortical, fullfile(out_dir_roi, [suffix '_meanW.nii']),      ROI_defs.(roi));
                    parcel_to_volume(mu_tmax, AtlasNii.Subcortical, fullfile(out_dir_roi, [suffix '_meanW_Tmax.nii']), ROI_defs.(roi));
                end
                parcel_to_Subsurface(mu,      fullfile(out_dir_roi, [suffix '_meanW.tiff']),      z_range, cmap_name);
                parcel_to_Subsurface(mu_tmax, fullfile(out_dir_roi, [suffix '_meanW_Tmax.tiff']), z_range, cmap_name);

            case 'AS200K34'
                if exist(AtlasNii.AS200K17_200ROI, 'file') && exist(pathAssign34,'file')
                    % Map 34 -> 200 for plotting
                    map = network_assignment.mapping;
                    mu200      = mu(map);
                    mu200_tmax = mu_tmax(map);

                    parcel_to_volume(mu200,      AtlasNii.AS200K17_200ROI, fullfile(out_dir_roi, [suffix '_meanW_on200.nii']), 1:200);
                    parcel_to_volume(mu200_tmax, AtlasNii.AS200K17_200ROI, fullfile(out_dir_roi, [suffix '_meanW_Tmax_on200.nii']), 1:200);

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
            case 'AS200K20'
                if exist(AtlasNii.AS200K17_200ROI, 'file') && exist(pathAssign34,'file')
                    % Map 34 -> 200 for plotting
                    map = network_assignment.mapping_20;
                    mu200      = mu(map);
                    mu200_tmax = mu_tmax(map);

                    parcel_to_volume(mu200,      AtlasNii.AS200K17_200ROI, fullfile(out_dir_roi, [suffix '_meanW_on200.nii']), 1:200);
                    parcel_to_volume(mu200_tmax, AtlasNii.AS200K17_200ROI, fullfile(out_dir_roi, [suffix '_meanW_Tmax_on200.nii']), 1:200);

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

save(pathOutStats, 'StatsGM', '-v7.3');
fprintf('\nDONE.\n  Raw Data: %s\n  Wscore:   %s\n  Stats:    %s\n', ...
    pathOutMeasures, pathOutWscore, pathOutStats);
