% =========================================================================
% Script: S11_pySuStaIn_statisitic.m
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
SuStaIn_GMV_dir = fullfile(result_dir, 'PySuStaln', 'GMV','K03');
SuStaIn_Correspondence_dir = fullfile(result_dir, 'PySuStaln', 'Correspondence','K02');
out_dir     = fullfile(result_dir, 'PySuStaln_statistic');

if ~exist(out_dir, 'dir'); mkdir(out_dir); end

addpath(genpath(fullfile(project_root, 'src', 'utils')));


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
    "Age", "HeadMotion", "EpilepsyType", "Lateralization", "Pathology", "FBTCS", ...
    "AgeOnset", "SeizureDuration"
    ];

% Define cognitive variables
cognitive_vars = [
    "IQ", "TMTA", "TMTB", "CVLT TL", "CVLT LDFR", "BNT", "Letter Fluency", "Semantic fluency", ...
    "Vocabulary", "Similarities", "Matrix Reasoning", "Digit Span", "Block Design", "Coding", ...
    "WCST PR", "WCST CC", "Logical Memory 1", "Logical Memory 2", "ROCF Copy", "Pegboard R", "Pegboard L"
    ];

% Define correspondence variables (names only)
correspondence_vars = [
    "Correspondence_VisualA"
    "Correspondence_VisualB"
    "Correspondence_SomatomotorA"
    "Correspondence_SomatomotorB"
    "Correspondence_DorsAttnB"
    "Correspondence_Sal_VenAttnA"
    "Correspondence_DorsAttnA"
    "Correspondence_TempPar"
    "Correspondence_Sal_VenAttnB"
    "Correspondence_ControlA"
    "Correspondence_LimbicA"
    "Correspondence_ControlC"
    "Correspondence_LimbicB"
    "Correspondence_DefaultC"
    "Correspondence_ControlB"
    "Correspondence_DefaultB"
    "Correspondence_DefaultA"];


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

corr_exist = correspondence_vars(ismember(correspondence_vars, varNames));

Groups.correspondence.data  = getCols(corr_exist);
Groups.correspondence.names = corr_exist(:);



%% =========================
% Part 3: Load SuStaIn Results and Create Stage Range Summaries
% =========================
filePath = fullfile(SuStaIn_GMV_dir, 'assignment_K3.csv');
opts = detectImportOptions(filePath, 'TextType', 'string');
opts.VariableNamingRule = 'preserve';
T_GMV = readtable(filePath, opts);
GMV_subtype = T_GMV.ml_subtype;
GMV_stage = T_GMV.ml_stage;

filePath = fullfile(SuStaIn_Correspondence_dir, 'assignment_K2.csv');
opts = detectImportOptions(filePath, 'TextType', 'string');
opts.VariableNamingRule = 'preserve';
T_Correspondence = readtable(filePath, opts);
Correspondence_subtype = T_Correspondence.ml_subtype;
Correspondence_stage = T_Correspondence.ml_stage;


mdl = fitlm(GMV_stage, Correspondence_stage);
R2  = mdl.Rsquared.Ordinary;

VIF_X1 = 1 / (1 - R2);

%% =========================
% Part 5: Correspondence Thresholding & Stacked Bar (Split by Subtype)
% =========================

% 1. Get Data
Z_corr = table2array(Groups.correspondence.data)*-1;
Z_names = strrep(Groups.correspondence.names, 'Correspondence_', '');
Z_names = erase(Z_names, "_")

% 2. Thresholding
Z_discrete = zeros(size(Z_corr));
Z_discrete(Z_corr < 1) = 0;
Z_discrete(Z_corr >= 1 & Z_corr < 2) = 1;
Z_discrete(Z_corr >= 2 & Z_corr < 3) = 2;
Z_discrete(Z_corr >= 3) = 3;

% Define Codes and Labels
codes_corr = [3, 2, 1, 0];
labels_corr = ["Z<-3", "-3<Z<-2", "-2<Z<-1", "Z>-1"];

% 3. Subtype Logic
% Subtype 1: Correspondence_subtype == 0
% Subtype 2: Correspondence_subtype == 1
% (Assuming Correspondence_subtype is loaded from Part 3)

subtypes = [0, 1];
subtype_names = ["Subtype 1 (transmodal-predominant)", "Subtype 2 (unimodal-predominant)"];

% Define Colors
% Subtype 1 (Red)
colors_s1 = [
    0.6, 0.0, 0.0;   % Z=3: Deep Red
    1.0, 0.4, 0.4;   % Z=2: Light Red
    1.0, 0.7, 0.7;   % Z=1: Lighter Red
    0.98, 0.90, 0.90 % Z=0: Very pale red
    ];

% Subtype 2 (Blue)
% Deep Blue -> Light Blue -> Lighter Blue -> Pale Blue
colors_s2 = [
    0.0, 0.2, 0.6;   % Z=3: Deep Blue (Navy)
    0.4, 0.6, 1.0;   % Z=2: Light Blue
    0.7, 0.8, 1.0;   % Z=1: Lighter Blue
    0.90, 0.94, 0.98 % Z=0: Very pale blue
    ];

colors_cell = {colors_s1, colors_s2};


for i = 1:numel(subtypes)
    st_code = subtypes(i);
    st_name = subtype_names(i);
    st_colors = colors_cell{i};

    % Filter subjects
    % Note: Correspondence_subtype comes from SuStaIn assignment.
    % Ensure Correspondence_subtype is valid (not NaN)
    idx_st = (Correspondence_subtype == st_code);

    if sum(idx_st) == 0
        fprintf('Warning: No subjects found for %s (Code %d)\n', st_name, st_code);
        continue;
    end

    Z_sub = Z_discrete(idx_st, :);

    % Reshape for plotting
    nSub_st = size(Z_sub, 1);

    X_flat = reshape(Z_sub, [], 1);
    Group_flat = repelem(string(Z_names), nSub_st, 1);
    Group_categorical = categorical(Group_flat, Z_names, 'Ordinal', true);

    % Plot
    % Width: 1200, YLim: 25%
    out_name = sprintf("StackedBar_Correspondence_%s.png", strrep(st_name, " ", ""));
    out_path = fullfile(out_dir, out_name);

    local_stacked_bar(Group_categorical, X_flat, st_name, ...
        out_path, codes_corr, labels_corr, st_colors, 1200, 60);

    % Save network values to CSV
    try
        M_sub = zeros(numel(Z_names), numel(codes_corr));
        for r_idx = 1:numel(Z_names)
            xi = Z_sub(:, r_idx);
            for c_idx = 1:numel(codes_corr)
                M_sub(r_idx, c_idx) = sum(xi == codes_corr(c_idx));
            end
        end
        P_sub = M_sub ./ max(sum(M_sub, 2), 1) * 100;
        
        numRows = numel(Z_names);
        numCols = 1 + 2 * numel(codes_corr) + 1;
        colNames = cell(1, numCols);
        colNames{1} = 'Network';
        for c_idx = 1:numel(codes_corr)
            labelStr = char(labels_corr{c_idx});
            colNames{2*c_idx} = ['Count_', labelStr];
            colNames{2*c_idx+1} = ['Percent_', labelStr];
        end
        colNames{end} = 'Total_Count';
        
        dataCell = cell(numRows, numCols);
        for r_idx = 1:numRows
            dataCell{r_idx, 1} = char(Z_names(r_idx));
            rowTotal = sum(M_sub(r_idx, :));
            for c_idx = 1:numel(codes_corr)
                dataCell{r_idx, 2*c_idx} = M_sub(r_idx, c_idx);
                dataCell{r_idx, 2*c_idx+1} = P_sub(r_idx, c_idx);
            end
            dataCell{r_idx, end} = rowTotal;
        end
        
        csv_path = fullfile(out_dir, sprintf("StackedBar_Correspondence_%s.csv", strrep(st_name, " ", "")));
        writecell([colNames; dataCell], csv_path);
        fprintf('  Saved stacked bar table to: %s\n', csv_path);
    catch ME
        warning('Failed to export stacked bar table: %s', ME.message);
    end

    fprintf('Generated stacked bar for %s (%d subjects)\n', st_name, nSub_st);
end


%% =========================
% Part 4: PCA of Cognition
% =========================

idx_keep = all(isfinite(table2array(Groups.cognitive.data)), 2);
idx_keep_list = find(idx_keep);

X = table2array(Groups.cognitive.data);
X = X(idx_keep, :);

[Xz, mu, sigma] = zscore(X);
sigma(sigma==0) = 1;

[coeff, score, latent, ~, explained, mu_pca] = pca(Xz, 'Algorithm','svd');

Y_keep = [Correspondence_stage,GMV_stage];
Y_keep = Y_keep(idx_keep_list,:)
[r,p] = corr(Correspondence_stage,GMV_stage,'Type','Spearman')

[r,p] = corr(score(:,1), Y_keep, 'Type','Spearman', 'Rows','complete')


Correspondence_subtype_keep = Correspondence_subtype(idx_keep_list)

%% =========================
% Part 6: Clinical Association
% =========================

%% -------------------------
% 6.1 Define SuStaIn groups
% -------------------------
% Correspondence: Normal / Subtype1 / Subtype2
CorrGroup = strings(height(T),1);
CorrGroup(Correspondence_stage == 0) = "Normal";
CorrGroup(Correspondence_stage > 0 & Correspondence_subtype == 0) = "Subtype1";
CorrGroup(Correspondence_stage > 0 & Correspondence_subtype == 1) = "Subtype2";
CorrGroup = categorical(CorrGroup, ["Normal","Subtype1","Subtype2"], 'Ordinal', true);

% GMV: Normal / Subtype1 / Subtype2 / Subtype3
GMVGroup = strings(height(T),1);
GMVGroup(GMV_stage == 0) = "Normal";
GMVGroup(GMV_stage > 0 & GMV_subtype == 0) = "Subtype1";
GMVGroup(GMV_stage > 0 & GMV_subtype == 1) = "Subtype2";
GMVGroup(GMV_stage > 0 & GMV_subtype == 2) = "Subtype3";
GMVGroup = categorical(GMVGroup, ["Normal","Subtype1","Subtype2","Subtype3"], 'Ordinal', true);


%% -------------------------
% 6.2 Clinical variables (numeric) + label dictionaries
% -------------------------
clin = struct();
clin.EpilepsyType    = double(T.EpilepsyType);      % 1=TLE, 2=EXE
clin.Lateralization  = double(T.Lateralization);    % 1=Left,2=Right,3=Unclear
clin.Pathology       = double(T.Pathology);         % 1=UHS,2=BHS,3=Lesion,4=Normal
clin.FTBTC           = double(T.FTBTC);             % 1=No,2=Yes
clin.AgeOnset        = double(T.AgeOnset);
clin.SeizureDuration = double(T.SeizureDuration);

LabelDict = struct();
LabelDict.EpilepsyType.codes   = [1 2];
LabelDict.EpilepsyType.text    = ["TLE","EXE"];

LabelDict.Lateralization.codes = [1 2 3];
LabelDict.Lateralization.text  = ["Left","Right","Unclear"];

LabelDict.Pathology.codes      = [1 2 3 4];
LabelDict.Pathology.text       = ["UHS","BHS","Lesion","Normal"];

LabelDict.FTBTC.codes          = [1 2];
LabelDict.FTBTC.text           = ["No","Yes"];

% =========================
% Color definitions (RGB)
% =========================

ColorDict = struct();

ColorDict.EpilepsyType = [
    0.15 0.35 0.70   % TLE  (deep blue)
    0.55 0.70 0.90   % EXE  (light blue)
    ];

ColorDict.Lateralization = [
    0.85 0.40 0.10   % Left   (deep orange)
    0.98 0.65 0.35   % Right  (light orange)
    0.70 0.70 0.70   % Unclear (gray)
    ];

ColorDict.Pathology = [
    0.45 0.20 0.60   % UHS    (deep purple)
    0.70 0.45 0.80   % BHS    (light purple)
    0.88 0.75 0.93   % Lesion (very light purple)
    0.75 0.75 0.75   % Normal (gray)
    ];
ColorDict.FTBTC = [
    0.70 0.15 0.20   % No   (deep red)
    0.95 0.55 0.55   % Yes  (light red)
    ];

%% -------------------------
% 6.3 Categorical variables: Chi-square + stacked bar
% -------------------------
catVars = ["EpilepsyType","Lateralization","Pathology","FTBTC"];

ResCorr = table();
ResGMV  = table();

for v = 1:numel(catVars)
    vn = catVars(v);
    x  = clin.(vn);

    % ----- Correspondence -----
    idx = ~isundefined(CorrGroup) & ~isnan(x);
    [tbl,~,~] = crosstab(CorrGroup(idx), x(idx));
    chi2 = sum((tbl - sum(tbl,2)*sum(tbl,1)/sum(tbl,'all')).^2 ./ ...
        (sum(tbl,2)*sum(tbl,1)/sum(tbl,'all')), 'all');
    dof  = (size(tbl,1)-1)*(size(tbl,2)-1);
    pval = 1 - chi2cdf(chi2,dof);

    ResCorr = [ResCorr; table("Corr_"+vn, pval, chi2, dof, ...
        join(string(LabelDict.(vn).codes)+"="+LabelDict.(vn).text,", "), ...
        'VariableNames',{'Test','p','chi2','dof','codebook'})];

    local_stacked_bar(CorrGroup, x, vn, ...
        fullfile(out_dir, "StackedBar_Corr_"+vn+".png"), ...
        LabelDict.(vn).codes, LabelDict.(vn).text,ColorDict.(vn));

    % ----- GMV -----
    idx = ~isundefined(GMVGroup) & ~isnan(x);
    [tbl,~,~] = crosstab(GMVGroup(idx), x(idx));
    chi2 = sum((tbl - sum(tbl,2)*sum(tbl,1)/sum(tbl,'all')).^2 ./ ...
        (sum(tbl,2)*sum(tbl,1)/sum(tbl,'all')), 'all');
    dof  = (size(tbl,1)-1)*(size(tbl,2)-1);
    pval = 1 - chi2cdf(chi2,dof);

    ResGMV = [ResGMV; table("GMV_"+vn, pval, chi2, dof, ...
        join(string(LabelDict.(vn).codes)+"="+LabelDict.(vn).text,", "), ...
        'VariableNames',{'Test','p','chi2','dof','codebook'})];

    local_stacked_bar(GMVGroup, x, vn, ...
        fullfile(out_dir, "StackedBar_GMV_"+vn+".png"), ...
        LabelDict.(vn).codes, LabelDict.(vn).text,ColorDict.(vn));
end


%% -------------------------
% 6.4 Continuous variables: Kruskal–Wallis + boxplot
% -------------------------
contVars = ["AgeOnset","SeizureDuration"];

for v = 1:numel(contVars)
    vn = contVars(v);
    y  = clin.(vn);

    % Correspondence
    idx = ~isundefined(CorrGroup) & ~isnan(y);
    [p, tbl] = kruskalwallis(y(idx), CorrGroup(idx), 'off');
    ResCorr = [ResCorr; table("Corr_"+vn+"_KW", p, tbl{2,5}, tbl{2,3}, ...
        "KruskalWallis", 'VariableNames',{'Test','p','chi2','dof','codebook'})];

    local_boxplot(CorrGroup, y, vn, fullfile(out_dir, "Box_Corr_"+vn+".png"));

    % GMV
    idx = ~isundefined(GMVGroup) & ~isnan(y);
    [p, tbl] = kruskalwallis(y(idx), GMVGroup(idx), 'off');
    ResGMV = [ResGMV; table("GMV_"+vn+"_KW", p, tbl{2,5}, tbl{2,3}, ...
        "KruskalWallis", 'VariableNames',{'Test','p','chi2','dof','codebook'})];

    local_boxplot(GMVGroup, y, vn, fullfile(out_dir, "Box_GMV_"+vn+".png"));
end


%% -------------------------
% 6.5 Save results
% -------------------------
writetable(ResCorr, fullfile(out_dir,'Clinical_Association.xlsx'), 'Sheet','Correspondence');
writetable(ResGMV,  fullfile(out_dir,'Clinical_Association.xlsx'), 'Sheet','GMV');


