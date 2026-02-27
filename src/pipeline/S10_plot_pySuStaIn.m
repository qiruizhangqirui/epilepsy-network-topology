% =========================================================================
% Script: S10_plot_pySuStaIn.m
% -------------------------------------------------------------------------
% Description:
%   Visualizes SuStaIn (Subtype and Stage Inference) results for epilepsy
%   cohorts across different disease progression models (K01-K04).
%
%   - Loads subject data and network W-scores from Excel.
%   - Groups features into clinical, cognitive, GMV, Correspondence, and Hubness.
%   - Loads SuStaIn assignment results (subtypes and stages).
%   - Creates stage range summaries and count matrices.
%   - Generates cortical and subcortical surface visualizations for:
%       1. Gray Matter Volume (GMV)
%       2. Hubness (network centrality)
%       3. Correspondence (network similarity)
%
% Inputs:
%   - data/Network_Wscore.xlsx (subject features)
%   - outputs/PySuStaln/GMV/K01-K04/assignment_K*.csv (SuStaIn results)
%   - Atlas files for cortical (AS200Y17) and subcortical regions
%   - Atlas/assignment_34.mat (network mapping)
%
% Outputs:
%   - outputs/PySuStaln/GMV/plot/K0X/*.tiff (Surface visualizations)
%   - outputs/PySuStaln/GMV/plot/K0X/*.nii (Volume files)
%   - outputs/PySuStaln/GMV/plot/K0X/*.csv (Summary tables)
%
% Author: Qirui Zhang, Farber Institute for Neuroscience, Department of Neurology, Thomas Jefferson University
% Email: qirui.zhang@jefferson.edu; fmrizhangqr@126.com
% Date: 12/24/2025
% =========================================================================

clear; clc;

%% =========================
% Part 1: Configuration & Paths
% =========================
this_file = mfilename('fullpath');

if ~isempty(this_file)
    % Case 1: script/function is being executed from a .m file
    project_root = fileparts(fileparts(fileparts(this_file)));
else
    % Case 2: executed from Command Window or interactive context
    project_root = fileparts(fileparts(pwd));
end

data_dir    = fullfile(project_root, 'data');
result_dir  = fullfile(project_root, 'outputs');
%SuStaIn_dir = fullfile(result_dir, 'PySuStaln', 'GMV'); % if plot GMV
SuStaIn_dir = fullfile(result_dir, 'PySuStaln', 'Correspondence'); % if plot Correspondence
out_dir     = fullfile(SuStaIn_dir, 'plot');

if ~exist(out_dir, 'dir'); mkdir(out_dir); end

addpath(genpath(fullfile(project_root, 'src', 'utils')));

% Atlas paths (using project root for portability)

Atlas200 = fullfile(data_dir, 'NCT_atlases', 'atlases', 'FSLMNI2mm', 'YeoLab', 'AS200Y17.nii.gz');
Atlas34  = fullfile(data_dir, 'AS17network_34ROIs_sorted.nii.gz');
Atlas_Subcortical = fullfile(data_dir, 'Subcortical_FLS2mm.nii');
assignment_file   = fullfile(data_dir, 'assignment_34.mat');

% Load network assignment mapping (34 ROIs to 200 parcels)
load(assignment_file, 'network_assignment');

% Define stage ranges for grouping
stage_ranges = {[0 0], [1 5], [6 10],[11 20], [21 50]};

%% =========================
% Part 2: Load Data and Define Groups
% =========================
% Read subject data from Excel
filePath = fullfile(data_dir, 'Network_Wscore.xlsx');
opts = detectImportOptions(filePath, 'TextType', 'string');
opts.VariableNamingRule = 'preserve';
T = readtable(filePath, opts);

% --- Subject Filtering ---
% Exclude subjects with T1notUse = 1
if ismember("T1notUse", T.Properties.VariableNames)
    n_before = height(T);
    T = T(T.T1notUse ~= 1, :);
    n_after = height(T);
    if n_before ~= n_after
        fprintf('[Filter] Excluded %d subjects with T1notUse=1. Remaining: %d\n', n_before - n_after, n_after);
    end
else
    fprintf('[Filter] Warning: Column ''T1notUse'' not found. Skipping exclusion.\n');
end

% % Exclude subjects with missing IQ (no cognitive data)
% if ismember("IQ", T.Properties.VariableNames)
%     n_before = height(T);
%     T = T(~isnan(T.IQ) & ~ismissing(T.IQ), :);
%     n_after = height(T);
%     if n_before ~= n_after
%         fprintf('[Filter] Excluded %d subjects with missing IQ. Remaining: %d\n', n_before - n_after, n_after);
%     end
% else
%     fprintf('[Filter] Warning: Column ''IQ'' not found. Skipping exclusion.\n');
% end

SubID    = T.("SubID");
varNames = string(T.Properties.VariableNames);

% Helper functions for data extraction
getByPrefix = @(prefix) varNames(startsWith(varNames, prefix));
getCols     = @(names) T(:, cellstr(names));
toMatrix    = @(tbl) table2array(varfun(@double, tbl));

% Define clinical variables
clinical_vars = [
    "Age", "HeadMotion", "EpilepsyType", "Lateralization", "Pathology", "FTBTC", ...
    "AgeOnset", "SeizureDuration"
    ];

% Define cognitive variables
cognitive_vars = [
    "IQ", "TMTA", "TMTB", "CVLT TL", "CVLT LDFR", "BNT", "Letter Fluency", "Semantic fluency", ...
    "Vocabulary", "Similarities", "Matrix Reasoning", "Digit Span", "Block Design", "Coding", ...
    "WCST PR", "WCST CC", "Logical Memory 1", "Logical Memory 2", "ROCF Copy", "Pegboard R", "Pegboard L"
    ];

% Extract feature groups by prefix
gmv_vars            = getByPrefix("GMV_");
Correspondence_vars = getByPrefix("Correspondence_");
hubness_vars        = getByPrefix("Hubness_");

% Build feature groups structure
Groups = struct();

% Clinical group (preserve original format including categorical)
clinical_exist = clinical_vars(ismember(clinical_vars, varNames));
Groups.clinical.data  = getCols(clinical_exist);
Groups.clinical.names = clinical_exist(:);

% Cognitive group
cog_exist = cognitive_vars(ismember(cognitive_vars, varNames));
Groups.cognitive.data  = getCols(cog_exist);
Groups.cognitive.names = cog_exist(:);

% GMV group
gmv_tbl = getCols(gmv_vars);
Groups.GMV.X     = toMatrix(gmv_tbl);
Groups.GMV.names = gmv_vars(:);

% Correspondence group
correspondence_tbl = getCols(Correspondence_vars);
Groups.Correspondence.X     = toMatrix(correspondence_tbl);
Groups.Correspondence.names = Correspondence_vars(:);

% Hubness group
hub_tbl = getCols(hubness_vars);
Groups.Hubness.X     = toMatrix(hub_tbl);
Groups.Hubness.names = hubness_vars(:);

%% =========================
% Part 3: Load SuStaIn Results and Create Stage Range Summaries
% =========================
Kfolders = {'K01', 'K02', 'K03', 'K04'};
All_K_Data = struct();

fprintf('\n=== Loading SuStaIn Assignment Results ===\n');

for k = 1:numel(Kfolders)
    Kdir = fullfile(SuStaIn_dir, Kfolders{k});
    files = dir(fullfile(Kdir, 'assignment_K*.csv'));

    allData = table();
    for i = 1:numel(files)
        fpath = fullfile(files(i).folder, files(i).name);
        T_sustain = readtable(fpath);

        % Record source file for tracking
        T_sustain.SourceFile = repmat(string(files(i).name), height(T_sustain), 1);

        allData = [allData; T_sustain];
    end

    % Group by subtype and stage
    G = findgroups(allData.ml_subtype, allData.ml_stage);

    summaryTbl = table;
    summaryTbl.Subtype = splitapply(@(x) x(1), allData.ml_subtype, G);
    summaryTbl.Stage   = splitapply(@(x) x(1), allData.ml_stage, G);

    % Collect SubID and subject_index lists (convert to 1-based indexing)
    summaryTbl.SubID_list = splitapply(@(x) {x}, allData.SubID, G);
    summaryTbl.SubjectIndex_list = splitapply(@(x) {x+1}, allData.subject_index, G);

    % Count subjects per group
    summaryTbl.Count = splitapply(@numel, allData.SubID, G);

    % Store in structure
    All_K_Data.(Kfolders{k}) = summaryTbl;

    fprintf('  Loaded %s: %d subjects, %d groups\n', Kfolders{k}, height(allData), height(summaryTbl));
end

%% Create Stage Range Summaries
fprintf('\n=== Creating Stage Range Summaries ===\n');


for k = 1:numel(Kfolders)
    outdir_perK = fullfile(out_dir, Kfolders{k});
    if ~exist(outdir_perK, 'dir'); mkdir(outdir_perK); end

    % Get summary table for this K
    Data = All_K_Data.(Kfolders{k});

    % Handle subtypes starting from 0
    Subtype_min = min(Data.Subtype);
    Subtype_max = max(Data.Subtype);

    result_cell = {};
    row = 1;

    for s = Subtype_min:Subtype_max
        for r = 1:numel(stage_ranges)
            range = stage_ranges{r};

            % Select subjects in this subtype & stage range
            idx = (Data.Subtype == s) & (Data.Stage >= range(1)) & (Data.Stage <= range(2));

            if any(idx)
                % Combine subject indices from all matching groups
                subj_list = vertcat(Data.SubjectIndex_list{idx});
            else
                subj_list = [];
            end

            % Store in cell array for table creation
            result_cell{row, 1} = s;
            result_cell{row, 2} = sprintf('%dto%d', range(1), range(2));
            result_cell{row, 3} = numel(subj_list);
            result_cell{row, 4} = subj_list(:)';
            row = row + 1;
        end
    end

    % Convert to table
    ResultTbl = cell2table(result_cell, ...
        'VariableNames', {'Subtype', 'StageRange', 'Count', 'SubjectIndex_list'});
    All_K_Data.([Kfolders{k}, '_Stage_range']) = ResultTbl;

    % Create wide-format count matrix (rows=Subtype, cols=Stage ranges)
    uniq_sub = (Subtype_min:Subtype_max).';
    colnames = cellfun(@(r) sprintf('%dto%d', r(1), r(2)), stage_ranges, 'UniformOutput', false);
    CountMat = array2table(nan(numel(uniq_sub), numel(stage_ranges)), ...
        'VariableNames', colnames);
    CountMat.Subtype = uniq_sub;
    CountMat = movevars(CountMat, 'Subtype', 'Before', 1);

    % Populate count matrix
    for i = 1:height(ResultTbl)
        s  = ResultTbl.Subtype(i);
        rg = ResultTbl.StageRange{i};
        c  = ResultTbl.Count(i);
        CountMat{CountMat.Subtype == s, rg} = c;
    end

    % Save tables
    outfile_long = fullfile(outdir_perK, sprintf('Subtype_StageRange_summary_%s.csv', Kfolders{k}));
    outfile_wide = fullfile(outdir_perK, sprintf('Subtype_StageRange_countWide_%s.csv', Kfolders{k}));
    writetable(ResultTbl, outfile_long);
    writetable(CountMat, outfile_wide);

    fprintf('  %s completed: %s | %s\n', Kfolders{k}, outfile_long, outfile_wide);
end

%% =========================
% Part 4: Generate Cortical and Subcortical Visualizations
% =========================
fprintf('\n=== Generating Surface Visualizations ===\n');

% Visualization parameters
color_range = [-1, 0];
cmap_name   = 'Blues_r';

for k = 1:numel(Kfolders)
    Data = All_K_Data.([Kfolders{k}, '_Stage_range']);
    outdir_perK = fullfile(out_dir, Kfolders{k});

    fprintf('  Processing %s...\n', Kfolders{k});

    for i = 0:max(Data.Subtype)
        range_types = unique(Data.StageRange);

        for j = 1:length(range_types)
            try
                % Get subject indices for this subtype and stage range
                idx_match = Data.Subtype == i & strcmp(Data.StageRange, range_types{j});
                Sub_index = Data.SubjectIndex_list{idx_match};

                if isempty(Sub_index); continue; end

                fprintf('    Subtype %d, Stage %s: %d subjects\n', i, range_types{j}, numel(Sub_index));

                % -------------------------------------------------
                % GMV Visualization
                % -------------------------------------------------
                GMV = Groups.GMV.X(Sub_index, :);
                Mean_GMV = mean(GMV, 1);  % 1-34: cortical, 47-60: subcortical

                % Cortical GMV (first 34 regions)
                Mean_GMV_C = Mean_GMV(1:34);
                mean_GMV_C_200 = Mean_GMV_C(network_assignment.mapping);

                parcel_to_volume(mean_GMV_C_200, Atlas200, ...
                    fullfile(outdir_perK, sprintf('GMV_Subtype%d_Stage%s_Cortical.nii', i, range_types{j})), 1:200);

                pv = parcel_to_surface(mean_GMV_C_200, 'schaefer_200x17_conte69');
                f = figure('Color', 'w', 'Position', [100 100 960 720]);
                plot_cortical(pv, 'surface_name', 'conte69', 'color_range', color_range, 'cmap', cmap_name);
                print(f, '-dtiff', '-r300', fullfile(outdir_perK, sprintf('GMV_Subtype%d_Stage%s_Cortical.tiff', i, range_types{j})));
                close(f);

                % Subcortical GMV (regions 47-60)
                Mean_GMV_S = Mean_GMV(47:60);
                parcel_to_volume(Mean_GMV_S, Atlas_Subcortical, ...
                    fullfile(outdir_perK, sprintf('GMV_Subtype%d_Stage%s_Subcortical.nii', i, range_types{j})), 1:14);
                parcel_to_Subsurface(Mean_GMV_S, ...
                    fullfile(outdir_perK, sprintf('GMV_Subtype%d_Stage%s_Subcortical.tiff', i, range_types{j})), ...
                    color_range, cmap_name);


                % -------------------------------------------------
                % Correspondence Visualization
                % -------------------------------------------------
                Correspondence = Groups.Correspondence.X(Sub_index, :);
                Mean_Correspondence = mean(Correspondence, 1);

                % Cortical Correspondence (17 networks, duplicated for bilateral)
                Mean_Correspondence_C = Mean_Correspondence(1:17);
                Mean_Correspondence_C = [Mean_Correspondence_C, Mean_Correspondence_C];
                mean_Correspondence_C_200 = Mean_Correspondence_C(network_assignment.mapping);

                parcel_to_volume(mean_Correspondence_C_200, Atlas200, ...
                    fullfile(outdir_perK, sprintf('Correspondence_Subtype%d_Stage%s_Cortical.nii', i, range_types{j})), 1:200);

                pv = parcel_to_surface(mean_Correspondence_C_200, 'schaefer_200x17_conte69');
                f = figure('Color', 'w', 'Position', [100 100 960 720]);
                plot_cortical(pv, 'surface_name', 'conte69', 'color_range', color_range, 'cmap', cmap_name);
                print(f, '-dtiff', '-r300', fullfile(outdir_perK, sprintf('Correspondence_Subtype%d_Stage%s_Cortical.tiff', i, range_types{j})));
                close(f);

                % Subcortical Correspondence (7 regions, duplicated for bilateral)
                Mean_Correspondence_S = Mean_Correspondence(18:24);
                Mean_Correspondence_S = [Mean_Correspondence_S, Mean_Correspondence_S];

                parcel_to_volume(Mean_Correspondence_S, Atlas_Subcortical, ...
                    fullfile(outdir_perK, sprintf('Correspondence_Subtype%d_Stage%s_Subcortical.nii', i, range_types{j})), 1:14);
                parcel_to_Subsurface(Mean_Correspondence_S, ...
                    fullfile(outdir_perK, sprintf('Correspondence_Subtype%d_Stage%s_Subcortical.tiff', i, range_types{j})), ...
                    color_range, cmap_name);

            catch ME
                warning('Failed to process Subtype %d, Stage %s: %s', i, range_types{j}, ME.message);
            end
        end
    end
end

fprintf('\nDONE. All visualizations saved to:\n  %s\n', out_dir);
