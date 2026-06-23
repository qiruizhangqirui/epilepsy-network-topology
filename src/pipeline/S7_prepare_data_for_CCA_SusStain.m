% =========================================================================
% Script: S7_prepare_data_for_CCA_SusStain.m
% -------------------------------------------------------------------------
% Description:
%   Prepares data for CCA and SusStain analysis (S7).
%
%   1. Loads Subjects.xlsx and filters data for specified sites:
%      - 'TJU': All patients (exclude HP).
%      - 'JLH': Only 'TLE' and 'EXE' patients.
%   2. Loads processed .mat files (W-scores, GM, Hubness).
%   3. Constructs a comprehensive table (T) with clinical info + features.
%   4. Generates dummy variables for categorical clinical features.
%   5. Saves results to 'Network_Wscore.xlsx' into separate sheets ('TJU', 'JLH')
%      with IDENTICAL headers (JLH Missing GM features are filled with NaN).
%
% Inputs:
%   - data/Subjects.xlsx
%   - data/assignment_34.mat
%   - outputs/Hubness_Wscore.mat
%   - outputs/Correspondence_Wscore.mat
%   - outputs/TJU_GM_Wscore.mat
%
% Outputs:
%   - data/Network_Wscore.xlsx (Sheets: 'TJU', 'JLH')
%
% Author: Qirui Zhang,  Farber Institute for Neuroscience, Department of Neurology, Thomas Jefferson University
% Email: qirui.zhang@jefferson.edu;fmrizhangqr@126.com
% Date: 12/24/2025 (Updated: 02/07/2026)
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
% Part 2: Load Global Data
% =========================
fprintf('Loading metadata & feature matrices...\n');
subjectTable = readtable(pathSubjects, 'VariableNamingRule','preserve');

load(pathHubness, 'WscoreHubness');
load(pathCorrespondence, 'Correspondence_Wscore');
load(pathGM, 'WscoreGM');
load(pathAssign34, 'network_assignment');

% Helper to check site column
if ismember('site', subjectTable.Properties.VariableNames)
    SiteCol = string(subjectTable.site);
else
    SiteCol = string(subjectTable.Site);
end
GroupCol = string(subjectTable.Group);

% Define Sites to Process
SitesToProcess = {'TJU', 'JLH'};

for iSite = 1:numel(SitesToProcess)
    
    TargetSite = SitesToProcess{iSite};
    fprintf('\nProcessing Site: %s ...\n', TargetSite);
    
    %% =========================
    % Part 3: Filter Subjects per Site
    % =========================
    idx_site = strcmpi(SiteCol, TargetSite);
    
    if strcmp(TargetSite, 'TJU')
        % Rule: Exclude HP
        isHP = strcmpi(strtrim(GroupCol), 'HP');
        keep_mask = idx_site & ~isHP;
        fprintf('  [TJU] Keeping Non-HP subjects.\n');
        
    elseif strcmp(TargetSite, 'JLH')
        % Rule: Keep only TLE and EXE
        isTLE = contains(GroupCol, 'TLE', 'IgnoreCase', true);
        isEXE = contains(GroupCol, 'EXE', 'IgnoreCase', true);
        keep_mask = idx_site & (isTLE | isEXE);
        fprintf('  [JLH] Keeping only TLE and EXE subjects.\n');
    else
        warning('Unknown site rule for %s. Skipping.', TargetSite);
        continue;
    end
    
    T = subjectTable(keep_mask, :);
    N = height(T);
    
    if N == 0
        warning('  No subjects found for %s with current rules. Skipping.', TargetSite);
        continue;
    end
    fprintf('  Selected %d subjects.\n', N);

    %% =========================
    % Part 4: Process Clinical Variables
    % =========================
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
    OutT.ASM_count       = double(T.ASM_count);       OutputLabels{end+1} = 'ASM_count';

    % -- Derived: Epilepsy Type --
    group_clean = string(T.Group);
    et = nan(N, 1);
    et(contains(group_clean, 'TLE')) = 1;
    et(contains(group_clean, 'EXE')) = 2;
    OutT.EpilepsyType = et;
    OutputLabels{end+1} = 'EpilepsyType';

    % -- Derived: Lateralization --
    lat_L = double(T.lateralization_L);
    lat_R = double(T.lateralization_R);
    lat_UC= double(T.lateralization_UC);
    lat = nan(N, 1);
    lat(lat_L==1) = 1;
    lat(lat_R==1) = 2;
    lat(lat_UC==1)= 3;
    OutT.Lateralization = lat;
    OutputLabels{end+1} = 'Lateralization';

    % -- Derived: Pathology --
    path_UHS    = double(T.Pathology_UHS);
    path_BHS    = double(T.Pathology_BHS);
    path_Lesion = double(T.Pathology_Lesion);
    path_Normal = double(T.Pathology_Normal);
    patho = nan(N, 1);
    patho(path_UHS==1)    = 1;
    patho(path_BHS==1)    = 2;
    patho(path_Lesion==1) = 3;
    patho(path_Normal==1) = 4;
    OutT.Pathology = patho;
    OutputLabels{end+1} = 'Pathology';

    % -- Derived: FTBTC --
    ftbtc = double(T.("FTBTC(Y/N)")) + 1;
    OutT.FTBTC = ftbtc;
    OutputLabels{end+1} = 'FTBTC';

    % -- Cognition --
    cognition_labels = { ...
        'IQ','TMTA','TMTB','CVLT TL','CVLT LDFR','BNT', ...
        'Letter Fluency','Semantic fluency','Vocabulary','Similarities', ...
        'Matrix Reasoning','Digit Span','Block Design','Coding', ...
        'WCST PR','WCST CC','Logical Memory 1','Logical Memory 2', ...
        'ROCF Copy','Pegboard R','Pegboard L'};

    for k = 1:numel(cognition_labels)
        varname = cognition_labels{k};
        validName = matlab.lang.makeValidName(varname);

        if ismember(varname, T.Properties.VariableNames)
            OutT.(validName) = double(T.(varname));
        else
            % warning('  Cognition variable "%s" not found. Filling NaN.', varname);
            OutT.(validName) = nan(N, 1);
        end
        OutputLabels{end+1} = varname;
    end

    %% =========================
    % Part 5: Extract Features
    % =========================
    
    % --- Helper to extract w-scores ---
    % S: Struct containing W-scores (e.g. WscoreGM.TJU.AS200K34)
    % mask: Logical mask for current site (keep_mask)
    % Note: W-scores struct usually stores data for ALL subjects of that site in rows.
    % We need to match subjects again by index or trust the order if it matches 'keep_mask' logic subset.
    % BUT: W-scores are stored as 'W' (N_site x Feat). We filtered 'T' from 'subjectTable' (Global).
    % So we need to map Global Subject Index -> Site Subject Index.
    
    % Re-identify indices relative to the SITE specific data matrix
    % GLOBAL:
    %   idx_site: logical, all subjects in this site
    %   keep_mask: logical, subset of idx_site (e.g. TLE only)
    
    % LOCAL (inside .mat struct for that site):
    %   Data usually corresponds to `subjectTable(idx_site,:)`.
    %   So we need indices of `keep_mask` RELATIVE to `idx_site`.
    
    idx_global = find(keep_mask);
    idx_site_all = find(idx_site);
    
    [~, idx_local] = ismember(idx_global, idx_site_all);
    % idx_local is the row indices in the specific Site's W-score matrices
    
    extract_local = @(S) S.W(idx_local, :);
    
    % --- GM Volumes ---
    roi_labels34 = network_assignment.network_order_name;
    roi_labels20 = network_assignment.network_order_name20;
    labels_sub   = {'Accumbens','Amygdala','Caudate','Hippocampus','Pallidum','Putamen','Thalamus'};
    labels_sub_LR= [strcat('L_', labels_sub), strcat('R_', labels_sub)];
    
    % Check if GM data exists for this site
    process_GM = false;
    if isfield(WscoreGM, TargetSite)
        process_GM = true;
    else
       fprintf('  [Info] No GM W-scores found for site %s. Filling with NaN.\n', TargetSite); 
    end

    % 1. GM 34
    if process_GM && isfield(WscoreGM.(TargetSite), 'AS200K34')
        GM_34 = extract_local(WscoreGM.(TargetSite).AS200K34);
    else
        GM_34 = nan(N, numel(roi_labels34));
    end
    
    for k = 1:numel(roi_labels34)
        nm = ['GMV_' roi_labels34{k}];
        vn = matlab.lang.makeValidName(nm);
        OutT.(vn) = GM_34(:, k);
        OutputLabels{end+1} = nm;
    end

    % 2. GM 20 (Add missing from 34)
    [isDup, ~] = ismember(roi_labels20, roi_labels34);
    idx_add20  = find(~isDup);
    
    if process_GM && isfield(WscoreGM.(TargetSite), 'AS200K20')
        GM_20 = extract_local(WscoreGM.(TargetSite).AS200K20);
    else
        GM_20 = nan(N, numel(roi_labels20));
    end

    for k = idx_add20(:)'
        nm = ['GMV_' roi_labels20{k}];
        vn = matlab.lang.makeValidName(nm);
        OutT.(vn) = GM_20(:, k);
        OutputLabels{end+1} = nm;
    end

    % 3. GM Subcortical
    if process_GM && isfield(WscoreGM.(TargetSite), 'Subcortical')
        GM_Sub = extract_local(WscoreGM.(TargetSite).Subcortical);
    else
        GM_Sub = nan(N, 14);
    end

    for k = 1:14
        nm = ['GMV_' labels_sub_LR{k}];
        vn = matlab.lang.makeValidName(nm);
        OutT.(vn) = GM_Sub(:, k);
        OutputLabels{end+1} = nm;
    end

    % --- Network Topology (Correspondence) ---
    Topology_labels = { ...
        'Normativity', ...
        'Non normativity', ...
        'TempPar','DefaultC','DefaultB','DefaultA','ControlC','ControlB','ControlA','LimbicA','LimbicB','Sal_VenAttnB','Sal_VenAttnA','DorsAttnB','DorsAttnA','SomatomotorB','SomatomotorA','VisualB','VisualA', ...
        'Accumbens','Amygdala','Caudate','Hippocampus','Pallidum','Putamen','Thalamus' ...
        };
    dice_labels17  = Topology_labels(3:19);
    dice_labelsSub = Topology_labels(20:end);
    
    if isfield(Correspondence_Wscore, TargetSite)
        C_Site = Correspondence_Wscore.(TargetSite);
        Normativity     = extract_local(C_Site.RawAtoms.Measure.Mean_Normativity);
        Non_normativity = extract_local(C_Site.RawAtoms.Measure.Mean_Non_normativity);
        Consensus_17    = extract_local(C_Site.RawAtoms.Consensus.AS200K17_17network);
        
        if isfield(C_Site.RawAtoms.Network_max_match, 'Subcortical')
             CNR_Sub = extract_local(C_Site.RawAtoms.Network_max_match.Subcortical);
        else
             % Fallback/Warning if Subcortical field is missing
             CNR_Sub = nan(N, 7);
        end
    else
         fprintf('  [Warning] No Correspondence W-scores for site %s. Filling NaN.\n', TargetSite);
         Normativity     = nan(N, 1);
         Non_normativity = nan(N, 1);
         Consensus_17    = nan(N, 17);
         CNR_Sub         = nan(N, 7);
    end

    OutT.Normativity = Normativity;
    OutputLabels{end+1} = 'Normativity';

    OutT.Non_normativity = Non_normativity;
    OutputLabels{end+1} = 'Non normativity';

    for k = 1:17
        nm = ['Correspondence_' dice_labels17{k}];
        vn = matlab.lang.makeValidName(nm);
        OutT.(vn) = Consensus_17(:, k);
        OutputLabels{end+1} = nm;
    end
    for k = 1:7
        nm = ['Correspondence_' dice_labelsSub{k}];
        vn = matlab.lang.makeValidName(nm);
        OutT.(vn) = CNR_Sub(:, k);
        OutputLabels{end+1} = nm;
    end

    % --- Hubness ---
    roi_labels34 = network_assignment.network_order_name; 
    
    if isfield(WscoreHubness, TargetSite)
        H_Site = WscoreHubness.(TargetSite);
        Hub_34  = extract_local(H_Site.RawAtoms.AS200K34);
        Hub_Sub = extract_local(H_Site.RawAtoms.Subcortical);
    else
        fprintf('  [Warning] No Hubness W-scores for site %s. Filling NaN.\n', TargetSite);
        Hub_34  = nan(N, numel(roi_labels34));
        Hub_Sub = nan(N, 14);
    end

    for k = 1:numel(roi_labels34)
        nm = ['Hubness_' roi_labels34{k}];
        vn = matlab.lang.makeValidName(nm);
        OutT.(vn) = Hub_34(:, k);
        OutputLabels{end+1} = nm;
    end
    for k = 1:14
        nm = ['Hubness_' labels_sub_LR{k}];
        vn = matlab.lang.makeValidName(nm);
        OutT.(vn) = Hub_Sub(:, k);
        OutputLabels{end+1} = nm;
    end

    %% =========================
    % Part 6: Generate Dummy Variables
    % =========================
    D1 = nan(N, 2);  % EpilepsyType
    D2 = nan(N, 3);  % Lateralization
    D3 = nan(N, 4);  % Pathology
    D4 = nan(N, 2);  % FTBTC

    % 1. EpilepsyType: 1=TLE, 2=EXE
    idx = ~isnan(OutT.EpilepsyType);
    if any(idx); D1(idx, :) = dummyvar(categorical(OutT.EpilepsyType(idx), [1 2])); end

    % 2. Lateralization: 1=Left, 2=Right, 3=Unclear
    idx = ~isnan(OutT.Lateralization);
    if any(idx); D2(idx, :) = dummyvar(categorical(OutT.Lateralization(idx), [1 2 3])); end

    % 3. Pathology: 1=UHS, 2=BHS, 3=Lesion, 4=Normal
    idx = ~isnan(OutT.Pathology);
    if any(idx); D3(idx, :) = dummyvar(categorical(OutT.Pathology(idx), [1 2 3 4])); end

    % 4. FTBTC: 1=No ("0"+1), 2=Yes ("1"+1)
    idx = ~isnan(OutT.FTBTC);
    if any(idx); D4(idx, :) = dummyvar(categorical(OutT.FTBTC(idx), [1 2])); end

    dummies = [D1 D2 D3 D4];
    dummy_names = [ ...
        strcat('EpilepsyType_', {'TLE_EXE','EXE_TLE'}), ...
        strcat('Lateralization_', {'Left','Right','Unclear'}), ...
        strcat('Pathology_', {'UHS','BHS','Lesion','Normal'}), ...
        {'FTBTC_No', 'FTBTC_Yes_No'} ...
        ];

    for k = 1:numel(dummy_names)
        vn = dummy_names{k};
        vns = matlab.lang.makeValidName(vn);
        OutT.(vns) = dummies(:, k);
        OutputLabels{end+1} = vn;
    end

    %% =========================
    % Part 7: Save Output to Sheet
    % =========================
    fprintf('  Saving sheet "%s" to %s ...\n', TargetSite, pathOutTable);

    DataCells = table2cell(OutT);
    FinalC = [OutputLabels; DataCells];
    
    writecell(FinalC, pathOutTable, 'Sheet', TargetSite);

end

fprintf('All Done.\n');