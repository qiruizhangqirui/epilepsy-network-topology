% =========================================================================
% Script: S7_prepare_data_for_CCA_SusStain.m
% -------------------------------------------------------------------------
% Description:
%   Prepares data for CCA and SusStain analysis (S7).
%
%   1. Loads Subjects.xlsx and filters for TJU site and Patients (excludes HP).
%   2. Loads processed .mat files (W-scores, GM, Hubness).
%   3. Constructs a comprehensive table (T) with clinical info + features.
%   4. Generates dummy variables for categorical clinical features.
%   5. Saves result to 'Network_Wscore.xlsx' with ORIGINAL headers (spaces preserved).
%
% Inputs:
%   - data/Subjects.xlsx
%   - data/assignment_34.mat
%   - outputs/Hubness_Wscore.mat
%   - outputs/Correspondence_Wscore.mat
%   - outputs/TJU_GM_Wscore.mat
%
% Outputs:
%   - data/Network_Wscore.xlsx
%
% Author: Qirui Zhang,  Farber Institute for Neuroscience, Department of Neurology, Thomas Jefferson University
% Email: qirui.zhang@jefferson.edu;fmrizhangqr@126.com
% Date: 12/24/2025
% =========================================================================

clear;

%% =========================
% Part 1: Configuration & Inputs
% =========================
this_file = mfilename('fullpath');
if ~isempty(this_file)
    project_root = fileparts(fileparts(fileparts(this_file)));
else
    project_root = fileparts(fileparts(pwd));
end

data_dir     = fullfile(project_root, 'data');
result_dir   = fullfile(project_root, 'outputs');

if ~exist(result_dir, 'dir'); mkdir(result_dir); end
addpath(genpath(fullfile(project_root, 'src', 'utils')));

% Input Paths
pathSubjects       = fullfile(data_dir, 'Subjects.xlsx');
pathHubness        = fullfile(result_dir, "Hubness_Wscore.mat");
pathCorrespondence = fullfile(result_dir, "Correspondence_Wscore.mat");
pathGM             = fullfile(result_dir, "TJU_GM_Wscore.mat");
pathAssign34       = fullfile(data_dir, "assignment_34.mat");

% Output Path
pathOutTable       = fullfile(data_dir, 'Network_Wscore.xlsx');

%% =========================
% Part 2: Load & Filter Metadata
% =========================
fprintf('Loading metadata...\n');
subjectTable = readtable(pathSubjects, 'VariableNamingRule','preserve');

% 1. Filter for TJU Site
if ismember('site', subjectTable.Properties.VariableNames)
    idx_TJU = strcmpi(string(subjectTable.site), 'TJU');
else
    % Fallback if capitalization differs
    idx_TJU = strcmpi(string(subjectTable.Site), 'TJU');
end
subjectTable = subjectTable(idx_TJU, :);

% 2. Filter out Healthy Controls (HP)
group_raw_tju = string(subjectTable.Group);
isHP_tju      = strcmpi(strtrim(group_raw_tju), 'HP');
keep_mask_tju = ~isHP_tju;  % Mask for slicing external matrices

T = subjectTable(keep_mask_tju, :);
fprintf('Selected %d TJU Patients (excluded %d HPs).\n', height(T), sum(isHP_tju));

%% =========================
% Part 3: Process Clinical Variables
% =========================
% We will build OutT (table with valid names) AND OutputLabels (cell array of real names)
OutT = table();
OutputLabels = {};

% -- Base Variables --
OutT.SubID           = string(T.ID_new);          OutputLabels{end+1} = 'SubID';
OutT.Sex             = double(T.("sex_b(M1F0)")); OutputLabels{end+1} = 'Sex';
OutT.Age             = double(T.age);             OutputLabels{end+1} = 'Age';
OutT.HeadMotion      = double(T.mean_FD);         OutputLabels{end+1} = 'HeadMotion';
OutT.AgeOnset        = double(T.("age of sz onset")); OutputLabels{end+1} = 'AgeOnset';
OutT.SeizureDuration = double(T.("sz duration")); OutputLabels{end+1} = 'SeizureDuration';
OutT.T1notUse        = double(T.("T1w_failedQC")); OutputLabels{end+1} = 'T1notUse';

% -- Derived: Epilepsy Type --
group_clean = string(T.Group);
et = nan(height(T), 1);
et(contains(group_clean, 'TLE')) = 1;
et(contains(group_clean, 'EXE')) = 2;
OutT.EpilepsyType = et;
OutputLabels{end+1} = 'EpilepsyType';

% -- Derived: Lateralization --
% 1=Left, 2=Right, 3=Unclear
lat_L = double(T.lateralization_L);
lat_R = double(T.lateralization_R);
lat_UC= double(T.lateralization_UC);
lat = nan(height(T), 1);
lat(lat_L==1) = 1;
lat(lat_R==1) = 2;
lat(lat_UC==1) = 3;
OutT.Lateralization = lat;
OutputLabels{end+1} = 'Lateralization';

% -- Derived: Pathology --
% 1=UHS, 2=BHS, 3=Lesion, 4=Normal
path_UHS    = double(T.Pathology_UHS);
path_BHS    = double(T.Pathology_BHS);
path_Lesion = double(T.Pathology_Lesion);
path_Normal = double(T.Pathology_Normal);
patho = nan(height(T), 1);
patho(path_UHS==1)    = 1;
patho(path_BHS==1)    = 2;
patho(path_Lesion==1) = 3;
patho(path_Normal==1) = 4;
OutT.Pathology = patho;
OutputLabels{end+1} = 'Pathology';

% -- Derived: FTBTC --
% 0->1(No), 1->2(Yes)
ftbtc = double(T.("FTBTC(Y/N)")) + 1;
OutT.FTBTC = ftbtc;
OutputLabels{end+1} = 'FTBTC';

% -- Cognition (Preserve Spaces) --
cognition_labels = { ...
    'IQ','TMTA','TMTB','CVLT TL','CVLT LDFR','BNT', ...
    'Letter Fluency','Semantic fluency','Vocabulary','Similarities', ...
    'Matrix Reasoning','Digit Span','Block Design','Coding', ...
    'WCST PR','WCST CC','Logical Memory 1','Logical Memory 2', ...
    'ROCF Copy','Pegboard R','Pegboard L'};

for i = 1:numel(cognition_labels)
    varname = cognition_labels{i};
    validName = matlab.lang.makeValidName(varname);

    if ismember(varname, T.Properties.VariableNames)
        OutT.(validName) = double(T.(varname));
    else
        warning('Cognition variable "%s" not found. Filling with NaN.', varname);
        OutT.(validName) = nan(height(T), 1);
    end
    OutputLabels{end+1} = varname; % Original name (with spaces)
end

%% =========================
% Part 4: Load & Merge Features
% =========================
fprintf('Loading feature matrices...\n');
load(pathHubness, 'WscoreHubness');
load(pathCorrespondence, 'Correspondence_Wscore');
load(pathGM, 'WscoreGM');
load(pathAssign34, 'network_assignment');

extract_W = @(S) S.W(keep_mask_tju, :);

% --- GM Volumes ---
roi_labels34 = network_assignment.network_order_name;
roi_labels20 = network_assignment.network_order_name20;
labels_sub   = {'Accumbens','Amygdala','Caudate','Hippocampus','Pallidum','Putamen','Thalamus'};
labels_sub_LR= [strcat('L_', labels_sub), strcat('R_', labels_sub)];

GM_34 = extract_W(WscoreGM.TJU.AS200K34);
GM_20 = extract_W(WscoreGM.TJU.AS200K20);
GM_Sub= extract_W(WscoreGM.TJU.Subcortical);

for i = 1:numel(roi_labels34)
    nm = ['GMV_' roi_labels34{i}];
    vn = matlab.lang.makeValidName(nm);
    OutT.(vn) = GM_34(:, i);
    OutputLabels{end+1} = nm;
end

[isDup, ~] = ismember(roi_labels20, roi_labels34);
idx_add20  = find(~isDup);
for k = idx_add20(:)'
    nm = ['GMV_' roi_labels20{k}];
    vn = matlab.lang.makeValidName(nm);
    OutT.(vn) = GM_20(:, k);
    OutputLabels{end+1} = nm;
end

for i = 1:14
    nm = ['GMV_' labels_sub_LR{i}];
    vn = matlab.lang.makeValidName(nm);
    OutT.(vn) = GM_Sub(:, i);
    OutputLabels{end+1} = nm;
end

% --- Network Topology ---
Topology_labels = { ...
    'Normativity', ...
    'Non normativity', ...
    'TempPar','DefaultC','DefaultB','DefaultA','ControlC','ControlB','ControlA','LimbicA','LimbicB','Sal_VenAttnB','Sal_VenAttnA','DorsAttnB','DorsAttnA','SomatomotorB','SomatomotorA','VisualB','VisualA', ...
    'Accumbens','Amygdala','Caudate','Hippocampus','Pallidum','Putamen','Thalamus' ...
    };
dice_labels17  = Topology_labels(3:19);
dice_labelsSub = Topology_labels(20:end);

Normativity     = extract_W(Correspondence_Wscore.TJU.RawAtoms.Measure.Mean_Normativity);
Non_normativity = extract_W(Correspondence_Wscore.TJU.RawAtoms.Measure.Mean_Non_normativity);

Consensus_17    = extract_W(Correspondence_Wscore.TJU.RawAtoms.Consensus.AS200K17_17network);
CNR_Sub         = extract_W(Correspondence_Wscore.TJU.RawAtoms.Network_max_match.Subcortical);

% Normativity
OutT.Normativity = Normativity;
OutputLabels{end+1} = 'Normativity';

% Non normativity (Internal: Safe, Output: Space)
OutT.Non_normativity = Non_normativity;
OutputLabels{end+1} = 'Non normativity';

for i = 1:17
    nm = ['Correspondence_' dice_labels17{i}];
    vn = matlab.lang.makeValidName(nm);
    OutT.(vn) = Consensus_17(:, i);
    OutputLabels{end+1} = nm;
end
for i = 1:7
    nm = ['Correspondence_' dice_labelsSub{i}];
    vn = matlab.lang.makeValidName(nm);
    OutT.(vn) = CNR_Sub(:, i);
    OutputLabels{end+1} = nm;
end

% --- Hubness ---
Hub_34  = extract_W(WscoreHubness.TJU.RawAtoms.AS200K34);
Hub_Sub = extract_W(WscoreHubness.TJU.RawAtoms.Subcortical);

for i = 1:numel(roi_labels34)
    nm = ['Hubness_' roi_labels34{i}];
    vn = matlab.lang.makeValidName(nm);
    OutT.(vn) = Hub_34(:, i);
    OutputLabels{end+1} = nm;
end
for i = 1:14
    nm = ['Hubness_' labels_sub_LR{i}];
    vn = matlab.lang.makeValidName(nm);
    OutT.(vn) = Hub_Sub(:, i);
    OutputLabels{end+1} = nm;
end

%% =========================
% Part 5: Generate Dummy Variables
% =========================
N = height(OutT);
D1 = nan(N, 2);  % EpilepsyType
D2 = nan(N, 3);  % Lateralization
D3 = nan(N, 4);  % Pathology
D4 = nan(N, 2);  % FTBTC

% 1. EpilepsyType: 1=TLE, 2=EXE
idx = ~isnan(OutT.EpilepsyType);
if any(idx)
    D1(idx, :) = dummyvar(categorical(OutT.EpilepsyType(idx), [1 2]));
end
% 2. Lateralization: 1=Left, 2=Right, 3=Unclear
idx = ~isnan(OutT.Lateralization);
if any(idx)
    D2(idx, :) = dummyvar(categorical(OutT.Lateralization(idx), [1 2 3]));
end
% 3. Pathology: 1=UHS, 2=BHS, 3=Lesion, 4=Normal
idx = ~isnan(OutT.Pathology);
if any(idx)
    D3(idx, :) = dummyvar(categorical(OutT.Pathology(idx), [1 2 3 4]));
end
% 4. FTBTC: 1=No ("0"+1), 2=Yes ("1"+1)
idx = ~isnan(OutT.FTBTC);
if any(idx)
    D4(idx, :) = dummyvar(categorical(OutT.FTBTC(idx), [1 2]));
end

dummies = [D1 D2 D3 D4];
dummy_names = [ ...
    strcat('EpilepsyType_', {'TLE_EXE','EXE_TLE'}), ...
    strcat('Lateralization_', {'Left','Right','Unclear'}), ...
    strcat('Pathology_', {'UHS','BHS','Lesion','Normal'}), ...
    {'FTBTC_No', 'FTBTC_Yes_No'} ...
    ];

for i = 1:numel(dummy_names)
    vn = dummy_names{i};
    vns = matlab.lang.makeValidName(vn);
    OutT.(vns) = dummies(:, i);
    OutputLabels{end+1} = vn;
end

%% =========================
% Part 6: Save Output
% =========================
% Use writecell to preserve headers with spaces
fprintf('Saving to %s ...\n', pathOutTable);

% Convert table data to cell
DataCells = table2cell(OutT);

% Final Cell Matrix
FinalC = [OutputLabels; DataCells];

writecell(FinalC, pathOutTable);
fprintf('Done.\n');