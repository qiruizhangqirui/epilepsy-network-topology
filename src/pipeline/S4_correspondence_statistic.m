% =========================================================================
% Script: S4_correspondence_statistic.m
% -------------------------------------------------------------------------
% Description:
%   Performs group-level statistical analysis on the W-scores computed in S3.
%
%   - Loads W-scores and subject metadata.
%   - Defines analysis groups (TLE, Focal, etc.).
%   - Generates:
%       1. Normativity/Non-normativity boxplots per atlas.
%       2. W-score Matrix plots throughout cohorts and atlases.
%       3. 17-network radar plots.
%       4. Subcortical surface maps.
%       5. 200-ROI surface maps (Permutation tests).
%
% Note:
%   - Requires 'permuztest' (permutation testing): https://github.com/mickcrosse/PERMUTOOLS
%   - Requires 'spider_plot' (radar plots): https://github.com/NewGuy012/spider_plot
%   - Requires 'ENIGMA' toolbox (plot_cortical, parcel_to_surface): https://github.com/MICA-MNI/ENIGMA
%   - Requires 'mafdr' (Bioinformatics Toolbox).
%
% Inputs:
%   - outputs/Correspondence_Wscore.mat
%   - outputs/All_derivatives_struct.mat
%   - data/Subjects.xlsx
%
% Outputs:
%   - outputs/correspondence_statistic/*.tiff (Figures)
%   - outputs/correspondence_statistic/AS200K17_200ROI_stats.mat
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

pathWscore      = fullfile(result_dir, 'Correspondence_Wscore.mat');
pathDerivatives = fullfile(result_dir, 'All_derivatives_struct.mat');
pathSubjects    = fullfile(data_dir, 'Subjects.xlsx');

out_dir = fullfile(result_dir, 'correspondence_statistic');
if ~exist(out_dir,'dir'); mkdir(out_dir); end

%% =========================
% Part 2: Load Data
% =========================
wscoreData = load(pathWscore);
wscores    = wscoreData.Correspondence_Wscore;

D = load(pathDerivatives);
Atlas_network_name = D.All_derivatives_struct.JLH0001; % Use first subject as template

subjectTable = readtable(pathSubjects,'VariableNamingRule','preserve');

subjectIDs = cellstr(string(subjectTable.("ID_new")));
site       = cellstr(string(subjectTable.("site")));
group_raw  = cellstr(string(subjectTable.("Group")));

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
% Part 4: Define Site Groups & Tables
% =========================
siteNames = {'JLH','TJU'};
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

        if strcmp(gr, 'HP')
            M(i, strcmp(columns,'HP')) = 1;
            continue
        end

        if strcmp(gr, 'TLE')
            if any(strcmp(columns,'TLE'))
                M(i, strcmp(columns,'TLE')) = 1;
            end
            if any(strcmp(columns,'Focal Epilepsy'))
                M(i, strcmp(columns,'Focal Epilepsy')) = 1;
            end
            if strcmp(gf,'TLE Left') && any(strcmp(columns,'TLE Left'))
                M(i, strcmp(columns,'TLE Left')) = 1;
            elseif strcmp(gf,'TLE Right') && any(strcmp(columns,'TLE Right'))
                M(i, strcmp(columns,'TLE Right')) = 1;
            end

        elseif strcmp(gr, 'EXE')
            if any(strcmp(columns,'EXE'))
                M(i, strcmp(columns,'EXE')) = 1;
            end
            if any(strcmp(columns,'Focal Epilepsy'))
                M(i, strcmp(columns,'Focal Epilepsy')) = 1;
            end

        elseif strcmp(gr, 'GGE') && any(strcmp(columns,'GGE'))
            M(i, strcmp(columns,'GGE')) = 1;

        elseif strcmp(gr, 'SeLECTS') && any(strcmp(columns,'SeLECTS'))
            M(i, strcmp(columns,'SeLECTS')) = 1;

        elseif strcmp(gr, 'AE') && any(strcmp(columns,'AE'))
            M(i, strcmp(columns,'AE')) = 1;
        end
    end

    groupTableBySite.(sitename) = array2table(M, 'VariableNames', columns);
end

% Subset group lists for plotting
PlotGroupsBySite = struct();
PlotGroupsBySite.TJU = {'Focal Epilepsy','TLE','EXE','TLE Left','TLE Right'};
PlotGroupsBySite.JLH = {'Focal Epilepsy','TLE','EXE','GGE','SeLECTS','AE'};

% Helper: take a site-specific GT and return a reduced GT containing only selected columns
get_GT_subset = @(GT, cols) GT(:, intersect(cols, GT.Properties.VariableNames, 'stable'));

%% =========================
% Part 5: Parameters & Types
% =========================
dataTypes  = fieldnames(wscores.TJU);
Atlas_type = {'AS200Y17','EG17','HCPICA','MG360J12','TY7','UKBICA','Subcortical'};
Natlas = numel(Atlas_type);

alpha_fdr   = 0.05;
alpha_maxT  = 0.05;
nperm       = 10000;
color_range = [-0.6 0.6];
color_range_hp =  [0.3 0.6];
cmap_name_hp =  'Reds';
cmap_name   = 'RdBu_r';

% Initialize Stats Storage
Correspondence_stats = struct();

%% =========================================================
% PART 1: Normativity / Non-normativity
% =========================================================
fprintf('\n=== PART 1: Normativity / Non-normativity (per atlas + mean) ===\n');

for iSite = 1:numel(siteNames)

    sitename = siteNames{iSite};
    idx_site = strcmp(site, sitename);

    columns = ColnamesBySite.(sitename);
    Ns = sum(idx_site);
    Ng = numel(columns);

    groupTable = groupTableBySite.(sitename);
    Group_names = groupTable.Properties.VariableNames;

    for iType = 1:numel(dataTypes)
        dataType = dataTypes{iType};
        fprintf('  Site=%s | DataType=%s\n', sitename, dataType);

        % -------------------------------------------------
        % 1A) Per-atlas Normativity boxplots
        % -------------------------------------------------
        if isfield(wscores.(sitename).(dataType), 'Measure') && ...
                isfield(wscores.(sitename).(dataType).Measure, 'Normativity') && ...
                isfield(wscores.(sitename).(dataType).Measure.Normativity, 'W')

            W = wscores.(sitename).(dataType).Measure.Normativity.W;  % Ns x Natlas
            for iAtlas = 1:Natlas
                w_vec = W(:,iAtlas);
                out = fullfile(out_dir, sprintf('%s_%s_Normativity_%s_box.tiff', sitename, dataType, Atlas_type{iAtlas}));
                label_str = sprintf('%s | %s | Normativity (%s)', sitename, dataType, Atlas_type{iAtlas});
                [w_med, p_val, labs, w_std] = plot_group_wscore(w_vec, groupTable, label_str, out);
                Correspondence_stats.(sitename).(dataType).Normativity.(Atlas_type{iAtlas}).w_median = w_med;
                Correspondence_stats.(sitename).(dataType).Normativity.(Atlas_type{iAtlas}).p_val = p_val;
                Correspondence_stats.(sitename).(dataType).Normativity.(Atlas_type{iAtlas}).w_std = w_std;
                Correspondence_stats.(sitename).(dataType).Normativity.(Atlas_type{iAtlas}).labels = labs;
                close all;
            end
        else
            warning('Missing %s.%s.Measure.Normativity.W', sitename, dataType);
        end

        % Mean Normativity boxplot
        if isfield(wscores.(sitename).(dataType), 'Measure') && ...
                isfield(wscores.(sitename).(dataType).Measure, 'Mean_Normativity') && ...
                isfield(wscores.(sitename).(dataType).Measure.Mean_Normativity, 'W')

            w_vec = wscores.(sitename).(dataType).Measure.Mean_Normativity.W; % Ns x 1
            out = fullfile(out_dir, sprintf('%s_%s_MeanNorm_box.tiff', sitename, dataType));
            label_str = sprintf('%s | %s | Mean Normativity', sitename, dataType);
            [w_med, p_val, labs, w_std] = plot_group_wscore(w_vec, groupTable, label_str, out);
            Correspondence_stats.(sitename).(dataType).Mean_Normativity.w_median = w_med;
            Correspondence_stats.(sitename).(dataType).Mean_Normativity.p_val = p_val;
            Correspondence_stats.(sitename).(dataType).Mean_Normativity.w_std = w_std;
            Correspondence_stats.(sitename).(dataType).Mean_Normativity.labels = labs;
            close all;
        else
            warning('Missing %s.%s.Measure.Mean_Normativity.W', sitename, dataType);
        end

        % -------------------------------------------------
        % 1B) Per-atlas Non-normativity boxplots
        % -------------------------------------------------
        if isfield(wscores.(sitename).(dataType), 'Measure') && ...
                isfield(wscores.(sitename).(dataType).Measure, 'Non_normativity') && ...
                isfield(wscores.(sitename).(dataType).Measure.Non_normativity, 'W')

            W = wscores.(sitename).(dataType).Measure.Non_normativity.W;  % Ns x Natlas
            for iAtlas = 1:Natlas
                w_vec = W(:,iAtlas);
                out = fullfile(out_dir, sprintf('%s_%s_NonNorm_%s_box.tiff', sitename, dataType, Atlas_type{iAtlas}));
                label_str = sprintf('%s | %s | Non-normativity (%s)', sitename, dataType, Atlas_type{iAtlas});
                [w_med, p_val, labs, w_std] = plot_group_wscore(w_vec, groupTable, label_str, out);
                Correspondence_stats.(sitename).(dataType).Non_normativity.(Atlas_type{iAtlas}).w_median = w_med;
                Correspondence_stats.(sitename).(dataType).Non_normativity.(Atlas_type{iAtlas}).p_val = p_val;
                Correspondence_stats.(sitename).(dataType).Non_normativity.(Atlas_type{iAtlas}).w_std = w_std;
                Correspondence_stats.(sitename).(dataType).Non_normativity.(Atlas_type{iAtlas}).labels = labs;
                close all;
            end
        else
            warning('Missing %s.%s.Measure.Non_normativity.W', sitename, dataType);
        end

        % Mean Non-normativity boxplot
        if isfield(wscores.(sitename).(dataType), 'Measure') && ...
                isfield(wscores.(sitename).(dataType).Measure, 'Mean_Non_normativity') && ...
                isfield(wscores.(sitename).(dataType).Measure.Mean_Non_normativity, 'W')

            w_vec = wscores.(sitename).(dataType).Measure.Mean_Non_normativity.W; % Ns x 1
            out = fullfile(out_dir, sprintf('%s_%s_MeanNonNorm_box.tiff', sitename, dataType));
            label_str = sprintf('%s | %s | Mean Non-normativity', sitename, dataType);
            [w_med, p_val, labs, w_std] = plot_group_wscore(w_vec, groupTable, label_str, out);
            Correspondence_stats.(sitename).(dataType).Mean_Non_normativity.w_median = w_med;
            Correspondence_stats.(sitename).(dataType).Mean_Non_normativity.p_val = p_val;
            Correspondence_stats.(sitename).(dataType).Mean_Non_normativity.w_std = w_std;
            Correspondence_stats.(sitename).(dataType).Mean_Non_normativity.labels = labs;
            close all;
        else
            warning('Missing %s.%s.Measure.Mean_Non_normativity.W', sitename, dataType);
        end

        % -------------------------------------------------
        % 1C) Wscore MATRIX PLOT (Normativity + Mean)
        % -------------------------------------------------
        fprintf('    -> Matrix plot: Normativity (+Mean)\n');

        w_mat = nan(Ng, Natlas+1);
        p_mat = nan(Ng, Natlas+1);

        if isfield(wscores.(sitename).(dataType), 'Measure') && ...
                isfield(wscores.(sitename).(dataType).Measure, 'Normativity') && ...
                isfield(wscores.(sitename).(dataType).Measure.Normativity, 'W')

            W = wscores.(sitename).(dataType).Measure.Normativity.W; % Ns x Natlas
            for iAtlas = 1:Natlas
                [w_median, p_val] = plot_group_wscore(W(:,iAtlas), groupTable, '', '');
                w_mat(:,iAtlas) = w_median;
                p_mat(:,iAtlas) = p_val;
                close all;
            end
        end

        if isfield(wscores.(sitename).(dataType), 'Measure') && ...
                isfield(wscores.(sitename).(dataType).Measure, 'Mean_Normativity') && ...
                isfield(wscores.(sitename).(dataType).Measure.Mean_Normativity, 'W')

            [w_median, p_val] = plot_group_wscore(wscores.(sitename).(dataType).Measure.Mean_Normativity.W, groupTable, '', '');
            w_mat(:,Natlas+1) = w_median;
            p_mat(:,Natlas+1) = p_val;
            close all;
        end

        % BH-FDR across all groups x (atlases+Mean)
        p_all = p_mat(:);
        q_all = nan(size(p_all));
        valid = ~isnan(p_all);
        q_all(valid) = mafdr(max(p_all(valid), realmin('double')), 'BHFDR', true);
        q_mat = reshape(q_all, size(p_mat));

        fig = figure('Color','w','Position',[100 100 900 650]);
        imagesc(w_mat, [-0.8 0.8]);

        if exist('RdBu_r','file') == 2
            colormap(RdBu_r);
        else
            colormap(parula);
            warning('RdBu_r not found. Using parula.');
        end
        colorbar;

        xlabels = [Atlas_type, {'Cortical Average'}];
        set(gca,'XTick',1:numel(xlabels),'XTickLabel',xlabels,'XTickLabelRotation',45,'FontSize',11, 'FontWeight','bold');
        set(gca,'YTick',1:Ng,'YTickLabel',Group_names,'FontSize',11, 'FontWeight','bold');
        title(sprintf('%s | %s | Normativity (median W-score)', sitename, dataType), 'Interpreter','none');

        for r = 1:Ng
            for c = 1:(Natlas+1)
                if isnan(w_mat(r,c)); continue; end
                stars = '';
                if ~isnan(q_mat(r,c)) && q_mat(r,c) < 0.01
                    stars = '**';
                elseif ~isnan(q_mat(r,c)) && q_mat(r,c) < 0.05
                    stars = '*';
                end

                text_color = 'w';
                if abs(w_mat(r,c)) < 0.4; text_color = [0.2, 0.2, 0.2]; end

                text(c, r, sprintf('%.2f%s', w_mat(r,c), stars), ...
                    'HorizontalAlignment','center','FontSize',11, 'FontWeight','bold','Color', text_color);
            end
        end

        out = fullfile(out_dir, sprintf('%s_%s_Normativity_WscoreMatrix.tiff', sitename, dataType));
        print(fig, '-dtiff', '-r300', out);
        close(fig);

        % -------------------------------------------------
        % 1D) Wscore MATRIX PLOT (Non-normativity + Mean)
        % -------------------------------------------------
        fprintf('    -> Matrix plot: Non-normativity (+Mean)\n');

        w_mat = nan(Ng, Natlas+1);
        p_mat = nan(Ng, Natlas+1);

        if isfield(wscores.(sitename).(dataType), 'Measure') && ...
                isfield(wscores.(sitename).(dataType).Measure, 'Non_normativity') && ...
                isfield(wscores.(sitename).(dataType).Measure.Non_normativity, 'W')

            W = wscores.(sitename).(dataType).Measure.Non_normativity.W;
            for iAtlas = 1:Natlas
                [w_median, p_val] = plot_group_wscore(W(:,iAtlas), groupTable, '', '');
                w_mat(:,iAtlas) = w_median;
                p_mat(:,iAtlas) = p_val;
                close all;
            end
        end

        if isfield(wscores.(sitename).(dataType), 'Measure') && ...
                isfield(wscores.(sitename).(dataType).Measure, 'Mean_Non_normativity') && ...
                isfield(wscores.(sitename).(dataType).Measure.Mean_Non_normativity, 'W')

            [w_median, p_val] = plot_group_wscore(wscores.(sitename).(dataType).Measure.Mean_Non_normativity.W, groupTable, '', '');
            w_mat(:,Natlas+1) = w_median;
            p_mat(:,Natlas+1) = p_val;
            close all;
        end

        p_all = p_mat(:);
        q_all = nan(size(p_all));
        valid = ~isnan(p_all);
        q_all(valid) = mafdr(max(p_all(valid), realmin('double')), 'BHFDR', true);
        q_mat = reshape(q_all, size(p_mat));

        fig = figure('Color','w','Position',[100 100 900 650]);
        imagesc(w_mat, [-0.8 0.8]);
        if exist('RdBu_r','file') == 2; colormap(RdBu_r); else; colormap(parula); end
        colorbar;

        xlabels = [Atlas_type, {'Cortical Average'}];
        set(gca,'XTick',1:numel(xlabels),'XTickLabel',xlabels,'XTickLabelRotation',45,'FontSize',11, 'FontWeight','bold');
        set(gca,'YTick',1:Ng,'YTickLabel',Group_names,'FontSize',11, 'FontWeight','bold');
        title(sprintf('%s | %s | Non_normativity (median W-score)', sitename, dataType), 'Interpreter','none');

        for r = 1:Ng
            for c = 1:(Natlas+1)
                if isnan(w_mat(r,c)); continue; end
                stars = '';
                if ~isnan(q_mat(r,c)) && q_mat(r,c) < 0.01; stars = '**';
                elseif ~isnan(q_mat(r,c)) && q_mat(r,c) < 0.05; stars = '*'; end

                text_color = 'w';
                if abs(w_mat(r,c)) < 0.4; text_color = [0.2, 0.2, 0.2]; end

                text(c, r, sprintf('%.2f%s', w_mat(r,c), stars), ...
                    'HorizontalAlignment','center','FontSize',11, 'FontWeight','bold','Color', text_color);
            end
        end

        out = fullfile(out_dir, sprintf('%s_%s_NonNorm_WscoreMatrix.tiff', sitename, dataType));
        print(fig, '-dtiff', '-r300', out);
        close(fig);
    end
end

%% =========================================================
% PART 2: 17-networks Radar & Boxplots
% =========================================================
fprintf('\n=== PART 2: 17-networks (box + radar) ===\n');

alpha_maxT_radar = 0.05;
nperm_radar      = 10000;
r_lim            = [-1.2 0.5];

for iSite = 1:numel(siteNames)

    sitename = siteNames{iSite};
    groupTable = groupTableBySite.(sitename);
    Group_names = groupTable.Properties.VariableNames;

    for iType = 1:numel(dataTypes)
        dataType = dataTypes{iType};

        if ~isfield(wscores.(sitename).(dataType),'Consensus') || ...
                ~isfield(wscores.(sitename).(dataType).Consensus,'AS200K17_17network') || ...
                ~isfield(wscores.(sitename).(dataType).Consensus.AS200K17_17network,'W')
            continue
        end

        Data17 = wscores.(sitename).(dataType).Consensus.AS200K17_17network.W;   % Ns x 17
        NetLabels = Atlas_network_name.AS200Y17.Networks;

        % A) 17-network boxplot
        GT_plot = get_GT_subset(groupTable, PlotGroupsBySite.(sitename));
        out = fullfile(out_dir, sprintf('%s_%s_Consensus17net_box.tiff', sitename, dataType));

        [p_val, p_maxT, w_med, w_std] = plot_group_wscore_by_network(Data17, GT_plot, NetLabels, out);
        Correspondence_stats.(sitename).(dataType).Consensus.AS200K17_17network.p_val = p_val;
        Correspondence_stats.(sitename).(dataType).Consensus.AS200K17_17network.p_maxT = p_maxT;
        Correspondence_stats.(sitename).(dataType).Consensus.AS200K17_17network.w_median = w_med;
        Correspondence_stats.(sitename).(dataType).Consensus.AS200K17_17network.w_std = w_std;
        Correspondence_stats.(sitename).(dataType).Consensus.AS200K17_17network.NetLabels = NetLabels;
        close all;

        % B) Radar plot per group
        for kk = 1:numel(Group_names)
            gname = Group_names{kk};
            idx_g = logical(groupTable.(gname));

            if sum(idx_g) < 3; continue; end

            G  = Data17(idx_g,:);
            mu = mean(G, 1, 'omitnan');

            sigma = std(G, 0, 1);
            sigma(sigma == 0) = eps;

            [~, p_perm] = permuztest(G, 0, sigma, ...
                'nperm', nperm_radar, 'alpha', alpha_maxT_radar, ...
                'tail', 'both', 'correct', true);

            NetStar = NetLabels;
            for nn = 1:numel(NetStar)
                if p_perm(nn) < 0.01; NetStar{nn} = [NetStar{nn}, '**'];
                elseif p_perm(nn) < 0.05; NetStar{nn} = [NetStar{nn}, '*']; end
            end
            NetStar = cellfun(@(x) strrep(x, '_', ''), NetStar, 'UniformOutput', false);

            mu_plot = [zeros(1, numel(mu)); mu];
            f = figure('Color','w','Position',[100 100 900 640]);
            spider_plot(mu_plot, ...
                'AxesLabels', NetStar, ...
                'AxesLimits', repmat(r_lim, numel(mu), 1)', ...
                'AxesInterval', 5, 'AxesPrecision', 2, 'AxesDisplay', 'one', ...
                'AxesFontSize', 10, 'LabelFontSize', 18, 'LineWidth', [1.5 2]);

            title(sprintf('%s | %s | %s | 17-net radar (max-T)', sitename, dataType, gname), 'Interpreter','none');
            out = fullfile(out_dir, sprintf('%s_%s_%s_Radar17.tiff', sitename, dataType, strrep(gname,' ','_')));
            print(f, '-dtiff', '-r300', out);
            close(f);
        end
    end
end

%% =========================================================
% PART 3: Canonical Network Representation (CNR)
% =========================================================
fprintf('\n=== PART 3: Canonical Network Representation ===\n');

out_dir_cnr = fullfile(out_dir, 'Canonical_Network_Representation');
if ~exist(out_dir_cnr, 'dir'); mkdir(out_dir_cnr); end

for iSite = 1:numel(siteNames)
    sitename = siteNames{iSite};
    groupTable = groupTableBySite.(sitename);

    for iType = 1:numel(dataTypes)
        dataType = dataTypes{iType};

        CNR = wscores.(sitename).(dataType).Network_max_match;
        tpl_fields = fieldnames(CNR);

        for tt = 1:numel(tpl_fields)
            tpl_field = tpl_fields{tt};
            DataNet = CNR.(tpl_field).W;

            NetLabels = Atlas_network_name.(tpl_fields{tt}).Networks;
            if strcmp(tpl_field, 'Subcortical'); NetLabels = NetLabels(1:7); end

            out = fullfile(out_dir_cnr, sprintf('%s_%s_CanonicalNetworkRepresentation_%s_box.tiff', sitename, dataType, tpl_field));
            GT_plot = get_GT_subset(groupTable, PlotGroupsBySite.(sitename));
            NetLabels = cellfun(@(x) strrep(x, '_', ''), NetLabels, 'UniformOutput', false);

            [p_val, pmaxT_matrix, w_med, w_std] = plot_group_wscore_by_network(DataNet, GT_plot, NetLabels, out);
            Correspondence_stats.(sitename).(dataType).CNR.(tpl_field).p_val = p_val;
            Correspondence_stats.(sitename).(dataType).CNR.(tpl_field).p_maxT = pmaxT_matrix;
            Correspondence_stats.(sitename).(dataType).CNR.(tpl_field).w_median = w_med;
            Correspondence_stats.(sitename).(dataType).CNR.(tpl_field).w_std = w_std;
            close all;

           % Subcortical surface plot
            if strcmpi(tpl_field, 'Subcortical')
                out_dir_sub = fullfile(out_dir_cnr, 'Subcortical_Surface');
                if ~exist(out_dir_sub, 'dir'); mkdir(out_dir_sub); end

                group_list = GT_plot.Properties.VariableNames;
                if size(DataNet,2) > numel(NetLabels); DataNet_use = DataNet(:, 1:numel(NetLabels));
                else; DataNet_use = DataNet; end

                for gg = 1:numel(group_list)
                    gname = group_list{gg};
                    idx_g = GT_plot.(gname) > 0;
                    if sum(idx_g) < 3; continue; end

                    mu = mean(DataNet_use(idx_g, :), 1, 'omitnan');

                    out_mu = fullfile(out_dir_sub, sprintf('%s_%s_CNR_%s_%s_meanW.tiff', sitename, dataType, tpl_field, strrep(gname,' ','_')));
                    parcel_to_Subsurface([mu,mu], out_mu, color_range, cmap_name);

                    pmaxT = pmaxT_matrix(gg, 1:numel(NetLabels));
                    h_maxT = pmaxT < alpha_maxT;
                    mu_maxT = mu; mu_maxT(~h_maxT) = 0;

                    out_mu_maxT = fullfile(out_dir_sub, sprintf('%s_%s_CNR_%s_%s_meanW_Tmax.tiff', sitename, dataType, tpl_field, strrep(gname,' ','_')));
                    parcel_to_Subsurface([mu_maxT,mu_maxT], out_mu_maxT, color_range, cmap_name);
                end

                % --- HP raw mean plot ---
                if isfield(CNR.(tpl_field), 'HP_mean')
                    mu_hp = CNR.(tpl_field).HP_mean;
                    if size(mu_hp,2) > numel(NetLabels); mu_hp = mu_hp(1:numel(NetLabels)); end
                    out_hp = fullfile(out_dir_sub, sprintf('%s_%s_CNR_%s_HP_meanW.tiff', sitename, dataType, tpl_field));
                    parcel_to_Subsurface([mu_hp, mu_hp], out_hp, color_range_hp, cmap_name_hp);
                end
            end
            close all;
        end
    end
end

%% =========================================================
% PART 4: 200 ROI consensus (AS200K17_200ROI)
% =========================================================
fprintf('\n=== PART 4: AS200K17_200ROI (permuztest + cortical plots) ===\n');

for iSite = 1:numel(siteNames)
    sitename = siteNames{iSite};
    groupTable = groupTableBySite.(sitename);
    Group_names = groupTable.Properties.VariableNames;

    for iType = 1:numel(dataTypes)
        dataType = dataTypes{iType};

        if ~isfield(wscores.(sitename).(dataType),'Consensus') || ...
                ~isfield(wscores.(sitename).(dataType).Consensus,'AS200K17_200ROI') || ...
                ~isfield(wscores.(sitename).(dataType).Consensus.AS200K17_200ROI,'W')
            continue
        end

        W200 = wscores.(sitename).(dataType).Consensus.AS200K17_200ROI.W;

        % --- HP raw mean plot ---
        if isfield(wscores.(sitename).(dataType).Consensus.AS200K17_200ROI, 'HP_mean')
            mu_hp = wscores.(sitename).(dataType).Consensus.AS200K17_200ROI.HP_mean;
            pv_hp = parcel_to_surface(mu_hp, 'schaefer_200x17_conte69');
            f_hp = figure('Color','w','Position',[100 100 960 720]);
            plot_cortical(pv_hp, 'surface_name','conte69', 'color_range', color_range_hp, 'cmap', cmap_name_hp);
            out_hp = fullfile(out_dir, sprintf('%s_%s_S200_HP_mean.tiff', sitename, dataType));
            print(f_hp, '-dtiff', '-r300', out_hp);
            close(f_hp);
        end

        for kk = 1:numel(Group_names)
            gname = Group_names{kk};
            idx_g = logical(groupTable.(gname));
            n = sum(idx_g);

            if n < 3; continue; end

            G  = W200(idx_g,:);
            mu = mean(G, 1, 'omitnan');

            % t-test + FDR
            [~, pval, ~, stats] = ttest(G, 0, 'Alpha', 0.05);
            pval = max(pval, realmin('double'));
            qval = mafdr(pval, 'BHFDR', true);

            % Permutation
            sigma = std(G, 0, 1);
            sigma(sigma == 0) = eps;

            [Z_perm, p_perm] = permuztest(G, 0, sigma, ...
                'nperm', nperm, 'alpha', alpha_maxT, 'tail', 'both', 'correct', true);

            h_perm = p_perm < alpha_maxT;

            % Plots
            pv = parcel_to_surface(mu, 'schaefer_200x17_conte69');
            f = figure('Color','w','Position',[100 100 960 720]);
            plot_cortical(pv, 'surface_name','conte69', 'color_range', color_range, 'cmap', cmap_name);
            out1 = fullfile(out_dir, sprintf('%s_%s_%s_S200_mean.tiff', sitename, dataType, strrep(gname,' ','_')));
            print(f, '-dtiff', '-r300', out1);
            close(f);

            mu_tmax = mu; mu_tmax(~h_perm) = 0;
            pv = parcel_to_surface(mu_tmax, 'schaefer_200x17_conte69');
            f = figure('Color','w','Position',[100 100 960 720]);
            plot_cortical(pv, 'surface_name','conte69', 'color_range', color_range, 'cmap', cmap_name);
            out3 = fullfile(out_dir, sprintf('%s_%s_%s_S200_Tmax.tiff', sitename, dataType, strrep(gname,' ','_')));
            print(f, '-dtiff', '-r300', out3);
            close(f);

            % Stats preservation
            gfield = matlab.lang.makeValidName(strrep(gname,' ','_'));
            StatOut.(sitename).(dataType).(gfield).n = n;
            StatOut.(sitename).(dataType).(gfield).mu = mu;
            StatOut.(sitename).(dataType).(gfield).t  = stats.tstat;
            StatOut.(sitename).(dataType).(gfield).p  = pval;
            StatOut.(sitename).(dataType).(gfield).q  = qval;
            StatOut.(sitename).(dataType).(gfield).perm.Z = Z_perm;
            StatOut.(sitename).(dataType).(gfield).perm.p = p_perm;
            StatOut.(sitename).(dataType).(gfield).perm.h = h_perm;
        end
    end
end
% 
% save(fullfile(out_dir, 'AS200K17_200ROI_stats.mat'), 'StatOut', '-v7.3');
Correspondence_stats.AS200K17_200ROI = StatOut;
save(fullfile(result_dir, 'Correspondence_stats.mat'), 'Correspondence_stats', '-v7.3');
fprintf('\nDONE. Outputs saved to:\n  %s\n', out_dir);

%% =========================================================
% PART 5: Standalone Binary Threshold Maps (TJU RawAtoms 200ROI)
% =========================================================
fprintf('\n=== PART 5: Binary Threshold Maps (TJU RawAtoms) ===\n');

st = 'TJU';
dt = 'RawAtoms';

if isfield(wscores.(st).(dt).Consensus, 'AS200K17_200ROI')
    W200 = wscores.(st).(dt).Consensus.AS200K17_200ROI.W;
    
    % 1) Focal Epilepsy group W < -0.7 (Blue)
    groupTable = groupTableBySite.(st);
    if any(strcmp(groupTable.Properties.VariableNames, 'Focal Epilepsy'))
        idx_focal = logical(groupTable.("Focal Epilepsy"));
        if sum(idx_focal) >= 3
            mu_focal = median(W200(idx_focal, :), 1, 'omitnan');
            mask_w = double(mu_focal < -0.5);
            
            pv = parcel_to_surface(mask_w, 'schaefer_200x17_conte69');
            f = figure('Color','w','Position',[100 100 960 720]);
            plot_cortical(pv, 'surface_name','conte69', 'color_range', [0 1], 'cmap', 'Blues');
            out_file = fullfile(out_dir, sprintf('%s_%s_FocalEpi_Wlt-05_mask.tiff', st, dt));
            print(f, '-dtiff', '-r300', out_file);
%             close(f);
        end
    end

    % 2) HP group mean > 0.53 (Red)
    if isfield(wscores.(st).(dt).Consensus.AS200K17_200ROI, 'HP_mean')
        mu_hp = wscores.(st).(dt).Consensus.AS200K17_200ROI.HP_mean;
        mask_hp = double(mu_hp > 0.50);
        
        pv = parcel_to_surface(mask_hp, 'schaefer_200x17_conte69');
        f = figure('Color','w','Position',[100 100 960 720]);
        plot_cortical(pv, 'surface_name','conte69', 'color_range', [0 1], 'cmap', 'Reds');
        out_file = fullfile(out_dir, sprintf('%s_%s_HPmean_gt05_mask.tiff', st, dt));
        print(f, '-dtiff', '-r300', out_file);
%         close(f);
    end
end
