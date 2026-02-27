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
%   - Requires 'permuztest': https://github.com/mickcrosse/PERMUTOOLS
%   - Requires 'ENIGMA' toolbox (plot_cortical, parcel_to_surface): https://github.com/MICA-MNI/ENIGMA
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


