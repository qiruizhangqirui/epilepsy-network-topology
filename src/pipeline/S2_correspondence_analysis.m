% =========================================================================
% Script: S2_correspondence_analysis.m
% -------------------------------------------------------------------------
% Description:
%   Calculates network correspondence measures between subject-specific
%   networks (derived from NCT/SPARK) and standard atlases.
%
%   - Loads aggregated derivatives (S1 output).
%   - Computes:
%       1. Normativity (Normv)
%       2. Non-normativity (NonNormv)
%       3. Network Max Match (NetMax)
%   - Generates consensus metrics across subjects.
%
% Inputs:
%   - outputs/All_derivatives_struct.mat
%   - data/Subjects.xlsx
%   - data/NCT_atlases/
%
% Outputs:
%   - outputs/Correspondence_measures.mat
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

% Load derivatives
% explicit assignment to avoid workspace pollution
derivativesData = load(fullfile(result_dir, 'All_derivatives_struct.mat'));
All_derivatives_struct = derivativesData.All_derivatives_struct;

%% =========================
% Part 2: Read Metadata
% =========================
subjectTable = readtable(fullfile(data_dir, 'Subjects.xlsx'), 'VariableNamingRule', 'preserve');

subjectIDs   = cellstr(string(subjectTable.("ID_new")));
Sex          = double(subjectTable.("sex_b(M1F0)")); %#ok<NASGU>
age          = double(subjectTable.("age"));         %#ok<NASGU>
HMotion      = double(subjectTable.("mean_FD"));     %#ok<NASGU>
Site         = cellstr(string(subjectTable.("site")));  %  JLH or TJU
group        = cellstr(string(subjectTable.("Group"))); %#ok<NASGU>

Nsub = numel(subjectIDs);

%% =========================
% Part 3: Atlas Configuration
% =========================
Atlas_type = {'AS200Y17','EG17','HCPICA','MG360J12','TY7','UKBICA','Subcortical'};
dataTypes  = {'RawAtoms', 'DenoisedAtoms'};

% Subcortical atlas: keep only first 7 columns (exclude Cerebellum)
include_cols = 1:7;

% Cortical atlases (exclude Subcortical) and their column indices
Atlas_type_cortical = {'AS200Y17','EG17','HCPICA','MG360J12','TY7','UKBICA'};
atlas_col_idx_cortical = 1:6;

% For across-atlas summary
is_cortical  = ~strcmp(Atlas_type, 'Subcortical');
cortical_idx = find(is_cortical);

%% =========================
% Part 4: Metrics Computation
% =========================
Measures = struct();
Natlas = numel(Atlas_type);

for iType = 1:numel(dataTypes)
    dataType = dataTypes{iType};

    Measures.(dataType).Network_max_match = cell(Nsub, Natlas);
    Measures.(dataType).Normativity       = nan(Nsub, Natlas);
    Measures.(dataType).Non_normativity   = nan(Nsub, Natlas);
    Measures.(dataType).Atom_number       = nan(Nsub, 1);
    Measures.(dataType).noise_number      = nan(Nsub, 1);
end

fprintf('\n[Part A] Computing normativity measures...\n');

for iAtlas = 1:Natlas
    atlasName = Atlas_type{iAtlas};
    is_subcortical = strcmp(atlasName, 'Subcortical');

    for iType = 1:numel(dataTypes)
        dataType = dataTypes{iType};

        for iSub = 1:Nsub
            subjID = subjectIDs{iSub};

            % -----------------------------
            % (1) Validate and extract per-subject / per-atlas data
            % -----------------------------
            if ~isfield(All_derivatives_struct, subjID) || ~isfield(All_derivatives_struct.(subjID), atlasName)
                Measures.(dataType).Network_max_match{iSub, iAtlas} = [];
                continue
            end

            Ssub = All_derivatives_struct.(subjID);

            Dice_coeff_raw = Ssub.(atlasName).Dice_coeff;  % [nAtoms x nTargets]

            % Atom count bookkeeping (raw atoms before cleaning)
            if isfield(Ssub.(atlasName), 'Atoms')
                atom_n = numel(Ssub.(atlasName).Atoms);
            else
                atom_n = size(Dice_coeff_raw, 1);
            end

            % Noise bookkeeping (subject-level)
            if isfield(Ssub, 'noise_idx')
                noise_idx = Ssub.noise_idx;
                noise_idx = noise_idx(~isnan(noise_idx));
                noise_n   = numel(noise_idx);
            else
                noise_idx = [];
                noise_n   = 0;
            end

            % -----------------------------
            % (2) Apply "DenoisedAtoms" (remove noisy atoms)
            % -----------------------------
            Dice_coeff = Dice_coeff_raw;
            if strcmpi(dataType, 'DenoisedAtoms') && ~isempty(noise_idx)
                keep_idx   = setdiff(1:size(Dice_coeff,1), noise_idx);
                Dice_coeff = Dice_coeff(keep_idx, :);
            end

            % -----------------------------
            % (3) Subcortical: exclude cerebellum columns
            % -----------------------------
            if is_subcortical
                Dice_coeff = Dice_coeff(:, include_cols);
            end

            % -----------------------------
            % (4) Compute metrics
            % -----------------------------
            [net_max, normv, nonnormv] = calc_normativity(Dice_coeff);

            % -----------------------------
            % (5) Store
            % -----------------------------
            Measures.(dataType).Network_max_match{iSub, iAtlas} = net_max;
            Measures.(dataType).Normativity(iSub, iAtlas)       = normv;
            Measures.(dataType).Non_normativity(iSub, iAtlas)   = nonnormv;

            Measures.(dataType).Atom_number(iSub,1)  = atom_n;
            Measures.(dataType).noise_number(iSub,1) = noise_n;
        end
    end
end

% Across-atlas summary (cortical only)
for iType = 1:numel(dataTypes)
    dataType = dataTypes{iType};
    Measures.(dataType).Mean_Normativity = mean( ...
        Measures.(dataType).Normativity(:, cortical_idx), 2, 'omitnan');
    Measures.(dataType).Mean_Non_normativity = mean( ...
        Measures.(dataType).Non_normativity(:, cortical_idx), 2, 'omitnan');
end

fprintf('[Part A] Done.\n');

%% =========================
% Part 5: Consensus Configuration
% =========================
fprintf('\n[Part B] Preparing consensus config...\n');

NCT_atlases_dir = fullfile(data_dir, 'NCT_atlases');
atlas_dir       = fullfile(NCT_atlases_dir, 'atlases', 'fs_LR_32k');
assign_dir      = fullfile(NCT_atlases_dir, 'network_assignment');

P = struct();
P.AS200Y17_mat        = fullfile(atlas_dir,  'YeoLab', 'AS200K17.mat');
P.AS200Y17_assign_mat = fullfile(assign_dir, 'AS200K17.mat');
P.TY7_mat             = fullfile(atlas_dir,  'YeoLab', 'TY7.mat');
P.EG17_mat            = fullfile(atlas_dir,  'WashU',  'EG17.mat');
P.MG360J12_mat        = fullfile(atlas_dir,  'Glasser','MG360J12.mat');

HCPICA_dir = fullfile(atlas_dir, 'HCPICA');
UKBICA_dir = fullfile(atlas_dir, 'UKBICA');

% ---------- Schaefer-200 labels in Conte69 (LH then RH) ----------
schaefer_csv = fullfile(data_dir, 'schaefer_200x17_conte69.csv');
labels_all   = csvread(schaefer_csv);
labels_all   = labels_all(:);

NPARC      = 200;
idx_nz     = labels_all > 0;
cnt_parcel = accumarray(labels_all(idx_nz), 1, [NPARC 1])';   % [1 x 200] vertices per parcel

% ---------- Load 200->17 mapping ----------
AS200Y17_map_file = fullfile(assign_dir, 'AS200Y17.mat');
S17 = load(AS200Y17_map_file);

map200to17   = S17.mapping(:)';              % [1 x 200]
net17_names  = S17.network_order_name(:)';   % [1 x 17]

% Pack into a single struct
cons = struct();
cons.P           = P;
cons.HCPICA_dir  = HCPICA_dir;
cons.UKBICA_dir  = UKBICA_dir;
cons.labels_all  = labels_all;
cons.NPARC       = NPARC;
cons.cnt_parcel  = cnt_parcel;
cons.map200to17  = map200to17;
cons.net17_names = net17_names;

fprintf('[Part B] Config ready.\n');

%% =========================
% Part 6: Consensus Computation
% =========================
fprintf('\n[Part B] Computing consensus match...\n');

for iType = 1:numel(dataTypes)
    dataType = dataTypes{iType};

    % Use only cortical atlases
    NetMaxCell = Measures.(dataType).Network_max_match(:, atlas_col_idx_cortical);

    [AS200K17_200ROI_Consensus, AS200K17_17network_Consensus, TemplateNames, nParcelPerNet17] = ...
        calc_consensus_normativity(NetMaxCell, Atlas_type_cortical, cons);

    % ---------- AS200K17 : 200 ROI level ----------
    Measures.(dataType).AS200K17_200ROI.Consensus = AS200K17_200ROI_Consensus;
    Measures.(dataType).AS200K17_200ROI.TemplateNames = TemplateNames;

    % ---------- AS200K17 : 17-network level ----------
    Measures.(dataType).AS200K17_17network.Consensus = AS200K17_17network_Consensus;
    Measures.(dataType).AS200K17_17network.NetID = (1:17)';
    Measures.(dataType).AS200K17_17network.NetAbbrev = cellstr(cons.net17_names(:));
    Measures.(dataType).AS200K17_17network.NParcelPerNet = nParcelPerNet17;

end

fprintf('[Part B] Done.\n');

%% =========================
% Part 7: Save
% =========================
% Rename back to Correspondence_measures for compatibility if needed,
% but using Measure logic internal to this script.
Correspondence_measures = Measures;

save(fullfile(result_dir, 'Correspondence_measures.mat'), 'Correspondence_measures', '-v7.3');
fprintf('\nSaved: %s\n', fullfile(result_dir, 'Correspondence_measures.mat'));
