% =========================================================================
% Script: S12_Validation.m
% -------------------------------------------------------------------------
% Description:
%   Validates the reproducibility of network findings between TJU and JLH
%   cohorts.
%
%   - Loads subject metadata (Subjects.xlsx) and feature W-scores.
%   - Splits data by site (TJU vs JLH).
%   - Groups subjects by epilepsy type (TLE, GGE, EXE, etc.).
%   - Computes group-level mean W-scores (Cortical + Subcortical).
%   - Performs cross-site correlation analysis (Scatter plots).
%   - Compares similarity and abnormality relative to Focal Epilepsy.
%   - Performs hierarchical clustering of epilepsy types.
%
% Note:
%   - Requires 'SpaitalCorrelation.m' in utils (Spatial permutation test)
%
% Inputs:
%   - data/Subjects.xlsx
%   - outputs/Correspondence_Wscore.mat
%   - outputs/Hubness_Wscore.mat
%
% Outputs:
%   - outputs/S12_Validation/*.png (Scatter plots, Dendrograms)
%
% Author: Qirui Zhang, Farber Institute for Neuroscience, Department of Neurology, Thomas Jefferson University
% Email: qirui.zhang@jefferson.edu; fmrizhangqr@126.com
% Date: 01/07/2026
% =========================================================================

clear; clc;

%% =========================
% Part 1: Configuration & Paths
% =========================
this_file = mfilename('fullpath');

if ~isempty(this_file)
    project_root = fileparts(fileparts(fileparts(this_file)));
else
    project_root = fileparts(fileparts(pwd));
end

data_dir    = fullfile(project_root, 'data');
result_dir  = fullfile(project_root, 'outputs');
out_dir     = fullfile(result_dir, 'Validation');

if ~exist(out_dir, 'dir'); mkdir(out_dir); end

addpath(genpath(fullfile(project_root, 'src', 'utils')));

% Input Paths
pathSubjects    = fullfile(data_dir, 'Subjects.xlsx');
pathCorr        = fullfile(result_dir, 'Correspondence_Wscore.mat');
pathHub         = fullfile(result_dir, 'Hubness_Wscore.mat');

%% =========================
% Part 2: Load Data
% =========================
fprintf('Loading data...\n');

% 1. Clinical Data
subjectTable = readtable(pathSubjects, 'VariableNamingRule', 'preserve');
subjectIDs   = cellstr(string(subjectTable.("ID_new")));
site         = cellstr(string(subjectTable.("site")));
Groups_raw   = cellstr(string(subjectTable.("Group")));

% 2. Feature Data (Correspondence & Hubness)
% Load full structures
tmpC = load(pathCorr);
Correspondence_Wscore = tmpC.Correspondence_Wscore;

tmpH = load(pathHub);
WscoreHubness = tmpH.WscoreHubness;

%% =========================
% Part 3: Define Groups (Logic from S4)
% =========================
% Refine TLE into Left/Right/Unclear/Unknown
Groups_ref = Groups_raw;
isTLE = strcmp(Groups_raw, 'TLE');

for i = 1:numel(Groups_raw)
    if isTLE(i)
        if subjectTable.("lateralization_L")(i) == 1
            Groups_ref{i} = 'TLE Left';
        elseif subjectTable.("lateralization_R")(i) == 1
            Groups_ref{i} = 'TLE Right';
        elseif subjectTable.("lateralization_UC")(i) == 1
            Groups_ref{i} = 'TLE Unclear';
        else
            Groups_ref{i} = 'TLE unknown';
        end
    end
end

% Define Site-specific Group Columns
siteNames = {'TJU', 'JLH'};
ColnamesBySite.TJU = {'Focal Epilepsy','TLE','EXE','TLE Left','TLE Right','HP'};
ColnamesBySite.JLH = {'Focal Epilepsy','TLE','EXE','GGE','SeLECTS','AE','TLE Left','TLE Right','HP'};

GroupTableBySite = struct();
SubjectIndicesBySite = struct();

for iSite = 1:numel(siteNames)
    sitename = siteNames{iSite};
    idx_site = strcmp(site, sitename);
    SubjectIndicesBySite.(sitename) = idx_site; % Store logical index for this site

    columns = ColnamesBySite.(sitename);
    Ns = sum(idx_site);

    M = zeros(Ns, numel(columns));

    % Get groups for this site's subjects
    g_raw = Groups_raw(idx_site);
    g_ref = Groups_ref(idx_site);

    for i = 1:Ns
        gr = g_raw{i};
        gf = g_ref{i};

        if strcmp(gr, 'HP')
            if any(strcmp(columns, 'HP')); M(i, strcmp(columns,'HP')) = 1; end
            continue
        end

        if strcmp(gr, 'TLE')
            if any(strcmp(columns,'TLE')); M(i, strcmp(columns,'TLE')) = 1; end
            if any(strcmp(columns,'Focal Epilepsy')); M(i, strcmp(columns,'Focal Epilepsy')) = 1; end

            if strcmp(gf,'TLE Left') && any(strcmp(columns,'TLE Left'))
                M(i, strcmp(columns,'TLE Left')) = 1;
            elseif strcmp(gf,'TLE Right') && any(strcmp(columns,'TLE Right'))
                M(i, strcmp(columns,'TLE Right')) = 1;
            end

        elseif strcmp(gr, 'EXE')
            if any(strcmp(columns,'EXE')); M(i, strcmp(columns,'EXE')) = 1; end
            if any(strcmp(columns,'Focal Epilepsy')); M(i, strcmp(columns,'Focal Epilepsy')) = 1; end

        elseif strcmp(gr, 'GGE') && any(strcmp(columns,'GGE'))
            M(i, strcmp(columns,'GGE')) = 1;

        elseif strcmp(gr, 'SeLECTS') && any(strcmp(columns,'SeLECTS'))
            M(i, strcmp(columns,'SeLECTS')) = 1;

        elseif strcmp(gr, 'AE') && any(strcmp(columns,'AE'))
            M(i, strcmp(columns,'AE')) = 1;
        end
    end

    GroupTableBySite.(sitename) = array2table(M, 'VariableNames', columns);
end

%% =========================
% Part 4: Calculate Group Means
% =========================
% Initialize storage for means
Means = struct();
% Structure level: Means.(Site).(FeatureType).(GroupName) -> [1 x (200+8)] vector
% FeatureType: 'Correspondence', 'Hubness'

DataTypeSets = {'Correspondence', 'Hubness'};

for iSite = 1:numel(siteNames)
    sitename = siteNames{iSite};
    GT = GroupTableBySite.(sitename);
    GroupNames = GT.Properties.VariableNames;

    % Extract raw data for this site (assumes Wscores are aligned with Subjects.xlsx order or we need to index)
    % NOTE: The Wscore structures usually contain data for specific atlases.
    % We need to access the correct field for the current site.
    % Based on S4, wscores has fields like 'TJU' and 'JLH'.

    for iType = 1:numel(DataTypeSets)
        dtype = DataTypeSets{iType};

        if strcmp(dtype, 'Correspondence')
            SourceStruct = Correspondence_Wscore.(sitename);
            % Path: .RawAtoms.Consensus.AS200K17_200ROI.W (Cortical)
            % Path: .RawAtoms.Network_max_match.Subcortical.W (Subcortical)

            if isfield(SourceStruct, 'RawAtoms') && isfield(SourceStruct.RawAtoms, 'Consensus')
                DataCort = SourceStruct.RawAtoms.Consensus.AS200K17_200ROI.W;
            else
                % Fallback or error if structure differs (check S11/S4)
                % Assuming standard structure based on prompt
                warning('Check structure for %s %s Cortical', sitename, dtype);
                DataCort = [];
            end

            if isfield(SourceStruct, 'RawAtoms') && isfield(SourceStruct.RawAtoms, 'Network_max_match')
                DataSub = SourceStruct.RawAtoms.Network_max_match.Subcortical.W;
            else
                DataSub = [];
            end

        elseif strcmp(dtype, 'Hubness')
            SourceStruct = WscoreHubness.(sitename);
            % Path: .RawAtoms.AS200K17_200ROI.W (Cortical)
            % Path: .RawAtoms.Subcortical.W (Subcortical)

            if isfield(SourceStruct, 'RawAtoms') && isfield(SourceStruct.RawAtoms, 'AS200K17_200ROI')
                DataCort = SourceStruct.RawAtoms.AS200K17_200ROI.W;
            else
                DataCort = [];
            end

            if isfield(SourceStruct, 'RawAtoms') && isfield(SourceStruct.RawAtoms, 'Subcortical')
                DataSub = SourceStruct.RawAtoms.Subcortical.W;
            else
                DataSub = [];
            end
        end

        % Combine (Concatenate along features)
        % DataCort: Ns x 200, DataSub: Ns x 8
        if isempty(DataCort) || isempty(DataSub)
            warning('Data missing for %s %s. Skipping.', sitename, dtype);
            continue;
        end

        DataAll = [DataCort, DataSub]; % Ns x 208

        % Calculate Means per Group
        for g = 1:numel(GroupNames)
            gname = GroupNames{g};
            gname_clean = strrep(gname, ' ', '_'); % standardizing field names

            idx_g = logical(GT.(gname));

            if sum(idx_g) > 0
                Means.(sitename).(dtype).(gname_clean) = mean(DataAll(idx_g, :), 1, 'omitnan');
            else
                Means.(sitename).(dtype).(gname_clean) = nan(1, size(DataAll,2));
            end
        end
    end
end

%% =========================
% Part 5: Validation (TJU vs JLH Scatter)
% =========================
% Overlapping groups to compare
CompareGroups = {'Focal_Epilepsy','TLE','EXE','TLE_Left','TLE_Right'};

% Plot Settings
metrics = {'Correspondence', 'Hubness'};
colors  = {[0.4, 0.6, 1.0], [ 1.0, 0.4, 0.4]}; % blue for corr, red for hub
StatsCell = {}; % To collect results

for iM = 1:numel(metrics)
    metric = metrics{iM};
    col    = colors{iM};

    for ig = 1:numel(CompareGroups)
        gname = CompareGroups{ig};

        % Get Means
        if isfield(Means.TJU.(metric), gname) && isfield(Means.JLH.(metric), gname)
            vec_TJU = Means.TJU.(metric).(gname);
            vec_JLH = Means.JLH.(metric).(gname);

            % --- 1. Standard Correlation (Full 208 regions) ---
            valid = ~isnan(vec_TJU) & ~isnan(vec_JLH);
            x = vec_TJU(valid);
            y = vec_JLH(valid);

            if isempty(x); continue; end

            [r_full, p_full] = corr(x', y', 'type', 'Pearson');

            % --- 2. Spin Test (using SpaitalCorrelation) ---
            % Used for both Cortical + Subcortical (Full brain)
            try
                % Ensure column vectors
                x_in = x(:);
                y_in = y(:);
                mask = ones(length(x_in), 1);
                nPer = 10000;

                [~, ~, ~, p_spin, ~, ~] = SpatialCorrelation(x_in, y_in, mask, nPer);
            catch ME
                warning('SpatialCorrelation failed for %s %s: %s', metric, gname, ME.message);
                p_spin = NaN;
            end

            % --- 3. Record Stats ---
            StatsCell{end+1, 1} = metric;
            StatsCell{end, 2}   = gname;
            StatsCell{end, 3}   = r_full;
            StatsCell{end, 4}   = p_full;
            StatsCell{end, 5}   = p_spin;

            % --- 4. Plot ---
            f = figure('Color','w','Position',[100 100 560 480]);
            scatter(x, y, 50, col, 'filled', 'MarkerFaceAlpha', 0.8); hold on;
            h = lsline; set(h, 'Color', 'k', 'LineWidth', 2.5);

            grid on;
            xlabel(sprintf('TJU %s (Mean W)', metric), 'FontSize', 24, 'FontWeight', 'bold');
            ylabel(sprintf('JLH %s (Mean W)', metric), 'FontSize', 24, 'FontWeight', 'bold');

            set(gca, 'FontSize', 18, 'LineWidth', 1.5, 'Box', 'on');

            % Title with Spin p-value
            if isnan(p_spin)
                title_str = sprintf('%s: %s (r=%.2f, p=%.3e)', strrep(gname,'_',' '), metric, r_full, p_full);
            else
                title_str = sprintf('%s: %s (r=%.2f, p_{spin}=%.3f)', strrep(gname,'_',' '), metric, r_full, p_spin);
            end
            title(title_str, 'FontSize', 14, 'FontWeight', 'bold');

            disk_name = sprintf('Scatter_TJUvsJLH_%s_%s.png', metric, gname);
            saveas(f, fullfile(out_dir, disk_name));
            close(f);

            % Save scatter plot raw data to CSV (TJU vs JLH)
            T_scatter = table(x(:), y(:), 'VariableNames', {['TJU_' metric], ['JLH_' metric]});
            csv_scatter = sprintf('Scatter_TJUvsJLH_%s_%s.csv', metric, gname);
            writetable(T_scatter, fullfile(out_dir, csv_scatter));
            fprintf('  Saved scatter plot raw data to: %s\n', csv_scatter);
        else
            fprintf('Skipping %s %s (Group missing)\n', metric, gname);
        end
    end
end

% Save Statistics Table
if ~isempty(StatsCell)
    T_Stats = cell2table(StatsCell, 'VariableNames', {'Metric', 'Group', 'R_Pearson', 'P_Parametric', 'P_Spin'});
    writetable(T_Stats, fullfile(out_dir, 'Validation_Statistics.xlsx'));
    disp('Statistics saved to Validation_Statistics.xlsx');
end

%% =========================
% Part 6: Comparison to Focal Epilepsy (JLH Only Example)
% =========================
% Reference: JLH Focal Epilepsy
% Compare other JLH groups to this reference
RefSite = 'JLH';
RefGroup = 'Focal_Epilepsy';
TargetGroups = {'TLE','EXE','GGE','SeLECTS','AE','TLE_Left','TLE_Right'};

if isfield(Means.(RefSite).Correspondence, RefGroup)
    ref_vec_corr = Means.(RefSite).Correspondence.(RefGroup);
    ref_vec_hub  = Means.(RefSite).Hubness.(RefGroup);

    r_corr = []; abn_corr = [];
    r_hub  = []; abn_hub  = [];
    labels = {};

    for ig = 1:numel(TargetGroups)
        gname = TargetGroups{ig};

        if isfield(Means.(RefSite).Correspondence, gname)
            % Correspondence
            vec = Means.(RefSite).Correspondence.(gname);
            valid = ~isnan(vec) & ~isnan(ref_vec_corr);
            r_corr(end+1)   = corr(vec(valid)', ref_vec_corr(valid)', 'type','Pearson');
            abn_corr(end+1) = mean(abs(vec(valid) - ref_vec_corr(valid)));

            % Hubness
            vec = Means.(RefSite).Hubness.(gname);
            valid = ~isnan(vec) & ~isnan(ref_vec_hub);
            r_hub(end+1)   = corr(vec(valid)', ref_vec_hub(valid)', 'type','Pearson');
            abn_hub(end+1) = mean(abs(vec(valid) - ref_vec_hub(valid)));

            labels{end+1} = strrep(gname, '_', ' ');
        end
    end

    % Plot Combined
    f = figure('Color','w', 'Position', [100 100 700 500]); hold on; box on; grid on;

    % Correspondence
    scatter(r_corr, abn_corr, 100, 'o', 'filled', 'MarkerFaceColor',[0 0.45 0.74], 'MarkerEdgeColor','k');
    % Hubness
    scatter(r_hub, abn_hub, 100, 'o', 'filled', 'MarkerFaceColor',[0.85 0.33 0.1], 'MarkerEdgeColor','k');

    % Labels
    for i = 1:numel(labels)
        text(r_corr(i)+0.02, abn_corr(i), labels{i}, 'Color',[0 0.45 0.74], 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
        text(r_hub(i), abn_hub(i)-0.012, labels{i}, 'Color',[0.85 0.33 0.1], 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    end

    xlabel('Correlation with Focal Epilepsy (Spatial Pattern)');
    ylabel('Mean Absolute Difference from Focal Epilepsy (Magnitude)');
    title('Similarity & Abnormality relative to Focal Epilepsy (JLH)');
    legend({'Correspondence','k-hubness'}, 'Location','best');
    set(gca,'FontSize',11);

    saveas(f, fullfile(out_dir, 'Similarity_Abnormality_JLH.png'));
    close(f);

    % Save Similarity & Abnormality data to CSV
    T_sa = table(labels', r_corr', abn_corr', r_hub', abn_hub', ...
        'VariableNames', {'Group', 'Correspondence_Similarity', 'Correspondence_Abnormality', 'k_hubness_Similarity', 'k_hubness_Abnormality'});
    writetable(T_sa, fullfile(out_dir, 'Similarity_Abnormality_JLH.csv'));
    fprintf('  Saved Similarity & Abnormality data to: Similarity_Abnormality_JLH.csv\n');
else
    warning('Reference group Focal_Epilepsy not found in Means.JLH');
end


% Groups to correlate
ClusterGroups = {'TLE_Left','TLE_Right','EXE','GGE','SeLECTS','AE'};
FeatureMat = [];
Labels = {};

for ig = 1:numel(ClusterGroups)
    gname = ClusterGroups{ig};
    if isfield(Means.JLH.Correspondence, gname)
        % Concatenate Correspondence + Hubness
        vec_c = Means.JLH.Correspondence.(gname);
        vec_h = Means.JLH.Hubness.(gname);

        if ~any(isnan(vec_c)) && ~any(isnan(vec_h))
            FeatureMat(end+1, :) = [vec_c, vec_h];
            Labels{end+1} = strrep(gname, '_', ' ');
        end
    end
end

if ~isempty(FeatureMat)
    %% =========================
    % Part 7: Feature Correlation Matrix
    % =========================
    % Calculate Correlation (Groups x Groups) with Spin Test
    nGroups = length(Labels);
    R_mat = zeros(nGroups);
    P_mat = zeros(nGroups); % This will store P_spin

    nPer = 10000;

    fprintf('Computing pairwise spin tests for %d groups...\n', nGroups);

    for i = 1:nGroups
        for j = 1:nGroups
            if i == j
                R_mat(i,j) = 1;
                P_mat(i,j) = 0;
            else
                % Extract vectors (Features x 1)
                vec_i = FeatureMat(i, :)';
                vec_j = FeatureMat(j, :)';
                mask = ones(length(vec_i), 1);

                % Call SpaitalCorrelation
                try
                    [r_val, ~, ~, p_spin, ~, ~] = SpatialCorrelation(vec_i, vec_j, mask, nPer);
                    R_mat(i,j) = r_val;
                    P_mat(i,j) = p_spin;
                catch
                    warning('SpatialCorrelation failed for %s vs %s', Labels{i}, Labels{j});
                    R_mat(i,j) = corr(vec_i, vec_j);
                    P_mat(i,j) = NaN;
                end
            end
        end
    end

    % Plot
    f = figure('Color','w','Position',[100 100 600 500]);

    % Set diagonal to 0 for visualization
    R_plot = R_mat;
    R_plot(logical(eye(size(R_mat)))) = 0;

    imagesc(R_plot);
    col = Reds;
    colormap(col); colorbar;
    caxis([0 1]);

    set(gca, 'XTick', 1:nGroups, 'XTickLabel', Labels, 'XTickLabelRotation', 45);
    set(gca, 'YTick', 1:nGroups, 'YTickLabel', Labels);
    set(gca, 'FontSize', 12, 'FontWeight', 'bold');
    title('Correlation Matrix of Epilepsy Types (JLH)');

    % Add Text (R lower, P upper)
    for i = 1:nGroups
        for j = 1:nGroups
            % Determine text color based on R value
            % r > 0.7 -> White, else Black
            val_for_color = R_mat(i,j);
            if val_for_color > 0.7
                txtCol = 'w';
            else
                txtCol = 'k';
            end

            if i > j % Lower Triangle: R
                text(j, i, sprintf('%.2f', R_mat(i,j)), 'HorizontalAlignment','center', 'Color',txtCol, 'FontWeight','bold');
            elseif i < j % Upper Triangle: P
                p_val = P_mat(i,j);
                if p_val < 0.001
                    txt = '<.001';
                else
                    txt = sprintf('%.3f', p_val);
                end
                text(j, i, txt, 'HorizontalAlignment','center', 'Color',txtCol, 'FontWeight','bold');
            else % Diagonal
                % Diagonal kept empty
            end
        end
    end

    saveas(f, fullfile(out_dir, 'Cluster_Correlation_Matrix_JLH.png'));
    close(f);

    % Save Result Table
    % Create a table with R and P values
    % Columns: Group1, Group2, R, P
    CorrResults = {};
    for i = 1:nGroups
        for j = i+1:nGroups
            CorrResults{end+1, 1} = Labels{i};
            CorrResults{end, 2}   = Labels{j};
            CorrResults{end, 3}   = R_mat(i,j);
            CorrResults{end, 4}   = P_mat(i,j);
        end
    end

    if ~isempty(CorrResults)
        T_Corr = cell2table(CorrResults, 'VariableNames', {'Group1', 'Group2', 'R', 'P'});
        writetable(T_Corr, fullfile(out_dir, 'Group_Correlation_JLH.xlsx'));
    end



    fprintf('Done. Results saved to %s\n', out_dir);
end
