function [pval_matrix, pmaxT_matrix, median_matrix, std_matrix] = plot_group_wscore_by_network( ...
    Data, GroupTable, Atlas_network_name, output_name)
% =========================================================================
% Function: plot_group_wscore_by_network
% -------------------------------------------------------------------------
% Description:
%   Draws W-score boxplots for each network across multiple groups, performing:
%     1) One-sample t-test vs 0 -> raw p-values
%     2) Max-T permutation correction (Family-wise correction) -> p_perm
%
% Inputs:
%   Data               - Nsub x Nnet W-score matrix
%   GroupTable         - Nsub x Ng Group Table (Columns = group names; >0 = member)
%   Atlas_network_name - 1 x Nnet cell array of network names
%   output_name        - Output file path (TIFF export recommended)
%
% Outputs:
%   pval_matrix        - Ng_valid x Nnet matrix of raw p-values
%   pmaxT_matrix       - Ng_valid x Nnet matrix of Tmax-corrected p-values
%   median_matrix      - Ng_valid x Nnet matrix of medians
%   std_matrix         - Ng_valid x Nnet matrix of standard deviations
%
% Notes:
%   - Significance annotation: p<0.05(*); Tmax<0.05(**); Tmax<0.01(***)
%   - Tmax settings: alpha_maxT=0.05, nperm=10000, two-tailed, family-wise correction
%   - Dependency: permuztest
%
% Author: Qirui Zhang, Farber Institute for Neuroscience, Thomas Jefferson University
% Date: 12/30/2025
% =========================================================================

% ---------------- Configuration ----------------
remove_outliers = true;     % Robust clipping of Y-axis (does not affect stats)
outlier_pct     = [0.1, 99.9];

alpha_maxT = 0.05;
nperm      = 10000;

% ---------------- Check Inputs ----------------
[Nsub, Nnet] = size(Data);
if height(GroupTable) ~= Nsub
    error('Data rows (%d) must match GroupTable rows (%d).', Nsub, height(GroupTable));
end
if numel(Atlas_network_name) ~= Nnet
    error('Atlas_network_name length (%d) must match Data columns (%d).', numel(Atlas_network_name), Nnet);
end

% ---------------- Select Valid Groups (Non-empty) ----------------
all_names   = GroupTable.Properties.VariableNames;
Ng_all      = numel(all_names);
valid_mask  = false(1, Ng_all);
for j = 1:Ng_all
    col = GroupTable{:, j};
    if ~isnumeric(col) && ~islogical(col)
        error('GroupTable column "%s" is not numeric or logical.', all_names{j});
    end
    m = (col > 0) & isfinite(col);
    valid_mask(j) = any(m);
end

if ~any(valid_mask)
    error('No valid groups found in GroupTable (all columns are 0/NaN).');
end

Group_labels = all_names(valid_mask);
GroupMat     = GroupTable{:, valid_mask};
Ngroup       = numel(Group_labels);

% ---------------- Initialization ----------------
colors         = lines(max(Ngroup, 5));
pval_matrix    = nan(Ngroup, Nnet);   % Raw p vs 0
pmaxT_matrix   = nan(Ngroup, Nnet);   % Tmax corrected p
group_medians  = nan(Ngroup, Nnet);   % For plotting
group_stds     = nan(Ngroup, Nnet);   % For stats
group_data     = cell(Ngroup,1);      % Store group data for Tmax

% X-axis setup
x_spacing   = 1.2;
x_positions = (0:(Nnet-1)) * x_spacing;
x_shift     = linspace(-0.4, 0.4, Ngroup);  % Group offset

% ---------------- Figure Canvas ----------------
figure('Units','normalized','Position',[0.05 0.1 0.9 0.5]); hold on;
box_width   = 0.15;
star_offset = 0.2;

% Grey reference line y = 0
plot([min(x_positions)-1, max(x_positions)+1], [0 0], ...
    'Color', [0.7 0.7 0.7], 'LineStyle', '-', 'LineWidth', 1.2);

% Initialize handles for legend
box_handles = gobjects(Ngroup,1);

% ---------------- Plotting & T-Test ----------------
for i = 1:Nnet
    for g = 1:Ngroup
        idx = GroupMat(:, g) > 0;           % >0 implies membership
        this_data = Data(idx, i);
        if isempty(this_data), continue; end

        % Store full network data for Tmax
        if i == 1
            group_data{g} = Data(idx, :);   % n X Nnet
        end

        % Plot
        x_pos = x_positions(i) + x_shift(g);
        h = boxchart(x_pos * ones(sum(idx),1), this_data, ...
            'BoxFaceColor', colors(g,:), ...
            'LineWidth', 1.2, ...
            'BoxWidth', box_width, ...
            'MarkerStyle', 'none');

        if i == 1
            box_handles(g) = h;
        end

        % One-sample t-test vs 0 (Raw p)
        [~, p] = ttest(this_data, 0);
        pval_matrix(g, i) = p;

        % Median for positioning stars
        group_medians(g, i) = median(this_data);
        group_stds(g, i)    = std(this_data);
    end
end

% ---------------- Tmax (max-T) Permutation (Group-wise) ----------------
for g = 1:Ngroup
    Xg = group_data{g};          % n x Nnet
    if isempty(Xg), continue; end

    % Standard deviation for sigma (avoid 0/inf)
    sigma = std(Xg, 0, 1);
    sigma(sigma == 0 | ~isfinite(sigma)) = eps;

    try
        % Two-tailed, Family-wise correction
        % Returns p_perm (1 x Nnet)
        [~, p_perm, ~, ~] = permuztest(Xg, 0, sigma, ...
            'nperm', nperm, 'alpha', alpha_maxT, ...
            'tail', 'both', 'correct', true);
    catch ME
        warning('Tmax permutation failed for group %s: %s. Check permuztest.', Group_labels{g}, ME.message);
        p_perm = nan(1, size(Xg,2));
    end

    pmaxT_matrix(g, :) = p_perm;
end

% ---------------- Significance Annotation (*, **, ***) ----------------
for i = 1:Nnet
    for g = 1:Ngroup
        p_raw  = pval_matrix(g, i);
        p_tmax = pmaxT_matrix(g, i);
        if isnan(p_raw) && isnan(p_tmax), continue; end

        if ~isnan(p_tmax) && p_tmax < 0.01
            stars = '**';                % Tmax < 0.01
        elseif ~isnan(p_tmax) && p_tmax < 0.05
            stars = '*';                  % 仅原始 p < 0.05
        else
            continue;
        end

        x_pos  = x_positions(i) + x_shift(g);
        y_star = group_medians(g, i) + star_offset;
        text(x_pos, y_star, stars, ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 11, ...
            'FontWeight', 'bold');
    end
end

% ---------------- Y-axis Clipping (Visualization Only) ----------------
all_vals = Data(:);
all_vals = all_vals(~isnan(all_vals));
if ~isempty(all_vals)
    if remove_outliers
        q_low  = prctile(all_vals, outlier_pct(1));
        q_high = prctile(all_vals, outlier_pct(2));
    else
        q_low  = min(all_vals);
        q_high = max(all_vals);
    end
    y_margin_ratio = 0.20;
    y_range = max(eps, q_high - q_low);
    ylim([q_low - y_range*y_margin_ratio, q_high + y_range*y_margin_ratio]);
end

% ---------------- Axes & Labels ----------------
xticks(x_positions);
xticklabels(Atlas_network_name);
xtickangle(45);
ylabel('W-score');
legend(box_handles, Group_labels, 'Location', 'best');
title('Group W-score per Network (t-test & max-T)');
xlim([min(x_positions)-0.8, max(x_positions)+0.8]);
box on;

% ---------------- Save Figure ----------------
if nargin >= 4 && ~isempty(output_name)
    [~,~,ext] = fileparts(output_name);
    if isempty(ext)
        print(gcf, output_name, '-dtiff', '-r600');
    else
        print(gcf, output_name, '-r600', ['-d' strrep(ext,'.','')]);
    end
end
% close(gcf);

% Outputs
median_matrix = group_medians;
std_matrix    = group_stds;

end
