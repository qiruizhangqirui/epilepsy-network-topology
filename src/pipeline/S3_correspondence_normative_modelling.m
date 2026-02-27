% =========================================================================
% Script: S3_correspondence_normative_modelling.m
% -------------------------------------------------------------------------
% Description:
%   Performs normative modeling on the network correspondence measures.
%   Uses PCNtoolkit (Python) to compute W-scores, using Healthy Controls
%   (HP) as the normative baseline.
%
% Inputs:
%   - outputs/Correspondence_measures.mat
%   - data/Subjects.xlsx
%
% Outputs:
%   - outputs/Correspondence_Wscore.mat
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

pathCorrespondence = fullfile(result_dir, 'Correspondence_measures.mat');
pathSubjects       = fullfile(data_dir, 'Subjects.xlsx');
pathOutput         = fullfile(result_dir, 'Correspondence_Wscore.mat');

if ~exist(result_dir,'dir'); mkdir(result_dir); end

%% =========================
% Part 2: Load Data
% =========================
S = load(pathCorrespondence);
assert(isfield(S,'Correspondence_measures'), ...
    'Missing Correspondence_measures in %s', pathCorrespondence);
Measures = S.Correspondence_measures;

dataTypes = fieldnames(Measures);   % {'RawAtoms','DenoisedAtoms'}

%% =========================
% Part 3: Read Metadata
% =========================
subjectTable = readtable(pathSubjects, 'VariableNamingRule','preserve');

subjectIDs = cellstr(string(subjectTable.("ID_new")));
site       = cellstr(string(subjectTable.("site")));
Sites      = {'JLH','TJU'};

sex = double(subjectTable.("sex_b(M1F0)"));
age = double(subjectTable.("age"));
hm  = double(subjectTable.("mean_FD"));

group_raw = cellstr(string(subjectTable.("Group")));

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

%% =========================
% Part 4: Normative Modeling (W-scores)
% =========================
% Computed by site (HP baseline within each site)
Correspondence_Wscore = struct();

for iSite = 1:numel(Sites)

    sitename = Sites{iSite};
    fprintf('\n=== Computing w-scores | Site: %s ===\n', sitename);

    idx_site = strcmp(site, sitename);
    idx_use  = idx_site;

    % HP baseline inside site
    idx_HP = idx_use & strcmp(group_ref,'HP');
    nHP = sum(idx_HP);
    fprintf('  HP baseline n = %d\n', nHP);

    for iType = 1:numel(dataTypes)

        dataType = dataTypes{iType};
        fprintf('  DataType: %s\n', dataType);

        % Atom number covariate
        if isfield(Measures.(dataType),'Atom_number')
            AtomN = double(Measures.(dataType).Atom_number);
        else
            AtomN = zeros(numel(subjectIDs),1);
        end

        X_cov  = [sex, age, hm, AtomN];
        X_site = X_cov(idx_use,:);
        idx_HP_site = idx_HP(idx_use);

        %% ---- 1) Normativity (per atlas) ----
        if isfield(Measures.(dataType),'Normativity')
            Y = Measures.(dataType).Normativity;
            Y_site = Y(idx_use,:);

            [W,Resid] = compute_wscore_pcn(X_site, Y_site, idx_HP_site);

            Correspondence_Wscore.(sitename).(dataType).Measure.Normativity.W     = W;
            Correspondence_Wscore.(sitename).(dataType).Measure.Normativity.Resid = Resid;
            Correspondence_Wscore.(sitename).(dataType).Measure.Normativity.HP_mean = mean(Y_site(idx_HP_site, :), 1, 'omitnan');
            Correspondence_Wscore.(sitename).(dataType).Measure.Normativity.HP_std  = std(Y_site(idx_HP_site, :), 0, 1, 'omitnan');
        end

        %% ---- 2) Non-normativity (per atlas) ----
        if isfield(Measures.(dataType),'Non_normativity')
            Y = Measures.(dataType).Non_normativity;
            Y_site = Y(idx_use,:);

            [W,Resid] = compute_wscore_pcn(X_site, Y_site, idx_HP_site);

            Correspondence_Wscore.(sitename).(dataType).Measure.Non_normativity.W     = W;
            Correspondence_Wscore.(sitename).(dataType).Measure.Non_normativity.Resid = Resid;
            Correspondence_Wscore.(sitename).(dataType).Measure.Non_normativity.HP_mean = mean(Y_site(idx_HP_site, :), 1, 'omitnan');
            Correspondence_Wscore.(sitename).(dataType).Measure.Non_normativity.HP_std  = std(Y_site(idx_HP_site, :), 0, 1, 'omitnan');
        end

        %% ---- 3) Mean_Normativity ----
        if isfield(Measures.(dataType),'Mean_Normativity')
            Y = Measures.(dataType).Mean_Normativity;
            Y_site = Y(idx_use,:);

            [W,Resid] = compute_wscore_pcn(X_site, Y_site, idx_HP_site);

            Correspondence_Wscore.(sitename).(dataType).Measure.Mean_Normativity.W     = W;
            Correspondence_Wscore.(sitename).(dataType).Measure.Mean_Normativity.Resid = Resid;
            Correspondence_Wscore.(sitename).(dataType).Measure.Mean_Normativity.HP_mean = mean(Y_site(idx_HP_site, :), 1, 'omitnan');
            Correspondence_Wscore.(sitename).(dataType).Measure.Mean_Normativity.HP_std  = std(Y_site(idx_HP_site, :), 0, 1, 'omitnan');
        end

        %% ---- 4) Mean_Non_normativity ----
        if isfield(Measures.(dataType),'Mean_Non_normativity')
            Y = Measures.(dataType).Mean_Non_normativity;
            Y_site = Y(idx_use,:);

            [W,Resid] = compute_wscore_pcn(X_site, Y_site, idx_HP_site);

            Correspondence_Wscore.(sitename).(dataType).Measure.Mean_Non_normativity.W     = W;
            Correspondence_Wscore.(sitename).(dataType).Measure.Mean_Non_normativity.Resid = Resid;
            Correspondence_Wscore.(sitename).(dataType).Measure.Mean_Non_normativity.HP_mean = mean(Y_site(idx_HP_site, :), 1, 'omitnan');
            Correspondence_Wscore.(sitename).(dataType).Measure.Mean_Non_normativity.HP_std  = std(Y_site(idx_HP_site, :), 0, 1, 'omitnan');
        end

        % ---- Network_max_match (per template) ----
        if isfield(Measures.(dataType),'Network_max_match') && ~isempty(Measures.(dataType).Network_max_match)

            NM = Measures.(dataType).Network_max_match;   % [Nsub × Ntpl] cell
            [Nsub, Ntpl] = size(NM);
            assert(Nsub == numel(subjectIDs), 'Network_max_match rows (%d) != #subjects (%d).', Nsub, numel(subjectIDs));

            % Template names
            Altas_type = {'AS200Y17','EG17','HCPICA','MG360J12','TY7','UKBICA','Subcortical'};
            tpl_names  = Altas_type;

            if isfield(Measures.(dataType),'TemplateNames')
                tn = Measures.(dataType).TemplateNames;
                if isstring(tn) || ischar(tn); tn = cellstr(string(tn)); end
                if iscell(tn) && numel(tn) == Ntpl
                    tpl_names = tn(:)';
                end
            end

            % Correspondence_Wscore.(sitename).(dataType).Network_max_match = struct();

            for iTpl = 1:Ntpl

                % --- infer network length for this template from first non-empty entry ---
                len = [];
                for ii0 = 1:Nsub
                    v0 = NM{ii0, iTpl};
                    if ~isempty(v0) && all(~isnan(v0(:)))
                        v0 = double(v0(:))';
                        len = numel(v0);
                        break;
                    end
                end
                if isempty(len)
                    warning('Network_max_match: Template %d has no valid entries. Skipped.', iTpl);
                    continue;
                end

                % --- build full matrix Y_full: [Nsub × len] ---
                Y_full = nan(Nsub, len);
                for ii = 1:Nsub
                    v = NM{ii, iTpl};
                    if isempty(v); continue; end
                    v = double(v(:))';   % force row
                    if numel(v) ~= len
                        % If mismatch occurs, pad or truncate
                        vv = nan(1, len);
                        vv(1:min(len,numel(v))) = v(1:min(len,numel(v)));
                        v = vv;
                    end
                    Y_full(ii,:) = v;
                end

                % --- site subset ---
                Y_site = Y_full(idx_use, :);

                % --- compute w-score using HP baseline within site ---
                [W, Resid] = compute_wscore_pcn(X_site, Y_site, idx_HP_site);

                % --- save ---
                tpl_field = matlab.lang.makeValidName(tpl_names{iTpl});
                Correspondence_Wscore.(sitename).(dataType).Network_max_match.(tpl_field).W     = W;
                Correspondence_Wscore.(sitename).(dataType).Network_max_match.(tpl_field).Resid = Resid;
                Correspondence_Wscore.(sitename).(dataType).Network_max_match.(tpl_field).HP_mean = mean(Y_site(idx_HP_site, :), 1, 'omitnan');
                Correspondence_Wscore.(sitename).(dataType).Network_max_match.(tpl_field).HP_std  = std(Y_site(idx_HP_site, :), 0, 1, 'omitnan');
            end
        end

        %% ---- 6) AS200K17 17-network consensus (N×17) ----
        if isfield(Measures.(dataType),'AS200K17_17network') && ...
                isfield(Measures.(dataType).AS200K17_17network,'Consensus')

            Y = Measures.(dataType).AS200K17_17network.Consensus;
            Y_site = Y(idx_use,:);

            [W,Resid] = compute_wscore_pcn(X_site, Y_site, idx_HP_site);

            Correspondence_Wscore.(sitename).(dataType).Consensus.AS200K17_17network.W     = W;
            Correspondence_Wscore.(sitename).(dataType).Consensus.AS200K17_17network.Resid = Resid;
            Correspondence_Wscore.(sitename).(dataType).Consensus.AS200K17_17network.HP_mean = mean(Y_site(idx_HP_site, :), 1, 'omitnan');
            Correspondence_Wscore.(sitename).(dataType).Consensus.AS200K17_17network.HP_std  = std(Y_site(idx_HP_site, :), 0, 1, 'omitnan');

            if isfield(Measures.(dataType).AS200K17_17network,'NetAbbrev')
                Correspondence_Wscore.(sitename).(dataType).Consensus.AS200K17_17network.FeatureNames = ...
                    cellstr(string(Measures.(dataType).AS200K17_17network.NetAbbrev));
            end
        end

        %% ---- 7) AS200K17 200-ROI consensus (N×200) ----
        if isfield(Measures.(dataType),'AS200K17_200ROI') && ...
                isfield(Measures.(dataType).AS200K17_200ROI,'Consensus')

            Y = Measures.(dataType).AS200K17_200ROI.Consensus;
            Y_site = Y(idx_use,:);

            [W,Resid] = compute_wscore_pcn(X_site, Y_site, idx_HP_site);

            Correspondence_Wscore.(sitename).(dataType).Consensus.AS200K17_200ROI.W     = W;
            Correspondence_Wscore.(sitename).(dataType).Consensus.AS200K17_200ROI.Resid = Resid;
            Correspondence_Wscore.(sitename).(dataType).Consensus.AS200K17_200ROI.HP_mean = mean(Y_site(idx_HP_site, :), 1, 'omitnan');
            Correspondence_Wscore.(sitename).(dataType).Consensus.AS200K17_200ROI.HP_std  = std(Y_site(idx_HP_site, :), 0, 1, 'omitnan');

            if isfield(Measures.(dataType).AS200K17_200ROI,'TemplateNames')
                Correspondence_Wscore.(sitename).(dataType).Consensus.AS200K17_200ROI.FeatureNames = ...
                    cellstr(string(Measures.(dataType).AS200K17_200ROI.TemplateNames));
            end
        end


        %% ---- Meta info ----
        Correspondence_Wscore.(sitename).(dataType).Meta.idx_use    = idx_use;
        Correspondence_Wscore.(sitename).(dataType).Meta.subjectIDs = subjectIDs(idx_use);
        Correspondence_Wscore.(sitename).(dataType).Meta.site       = site(idx_use);
        Correspondence_Wscore.(sitename).(dataType).Meta.group_ref  = group_ref(idx_use);
        Correspondence_Wscore.(sitename).(dataType).Meta.idx_HP     = idx_HP_site;
        Correspondence_Wscore.(sitename).(dataType).Meta.covariates = ...
            {'sex','age','mean_FD','Atom_number'};
    end
end

%% =========================
% Part 5: Save
% =========================
save(pathOutput, 'Correspondence_Wscore', '-v7.3');
fprintf('\nSaved w-score only: %s\n', pathOutput);
