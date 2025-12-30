function [w_vec, p_vec, labels, std_vec] = plot_group_wscore(w_scores, GroupTable, label_str, output_name)
% =========================================================================
% Function: plot_group_wscore
% -------------------------------------------------------------------------
% Description:
%   Generates a boxplot of a single W-score metric across multiple groups.
%   Performs one-sample t-tests against zero and BH-FDR correction.
%   Annotates significance with stars (** q<0.01, * q<0.05).
%
% Inputs:
%   w_scores    - N x 1 column vector of W-scores
%   GroupTable  - Table where each column is a group indicator
%   label_str   - (Optional) String for plot title annotation
%   output_name - (Optional) Output filename for saving the figure
%
%   w_vec       - Per-group median (vector)
%   p_vec       - Per-group raw p-value (vector)
%   labels      - Group names (cell array)
%   std_vec     - Per-group standard deviation (vector)
%
% Author: Qirui Zhang, Farber Institute for Neuroscience, Thomas Jefferson University
% Date: 12/30/2025
% =========================================================================

% ===== Configuration =====
box_width       = 0.4;
star_offset     = 0.15;
y_margin_ratio  = 0.20;
outlier_pct_lim = [0.1, 99.9];

% ===== Input checks =====
if size(w_scores, 2) ~= 1
    error('w_scores must be an N×1 column vector.');
end
if height(GroupTable) ~= size(w_scores, 1)
    error('Row count mismatch: w_scores and GroupTable must have the same number of rows.');
end

% ===== Select valid columns (not all 0/NaN) =====
varNames = GroupTable.Properties.VariableNames;
is_valid_col = false(1, numel(varNames));

for j = 1:numel(varNames)
    col = GroupTable{:, j};

    if ~isnumeric(col) && ~islogical(col)
        error('Column "%s" in GroupTable is not numeric or logical.', varNames{j});
    end

    mask = (col > 0) & isfinite(col);  % >0 indicates membership; require finite
    is_valid_col(j) = any(mask);
end

if ~any(is_valid_col)
    error('No valid groups found in GroupTable (all columns are all-0 / NaN / empty).');
end

varNames = varNames(is_valid_col);
GroupMat = GroupTable{:, is_valid_col};
Ngroup   = numel(varNames);

colors = lines(max(Ngroup, 5));

% ===== Stats containers =====
w_med       = nan(1, Ngroup);
w_std       = nan(1, Ngroup);
p_vals      = nan(1, Ngroup);
box_handles = gobjects(Ngroup, 1);

% ===== Figure canvas =====
figure('Units', 'normalized', 'Position', [0.15 0.2 0.38 0.56]); hold on;
plot([0.5, Ngroup + 0.5], [0 0], 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);

% ===== Plot each group =====
for g = 1:Ngroup
    idx = GroupMat(:, g) > 0;   % allow non-binary membership; >0 treated as in-group
    xg  = w_scores(idx);
    if isempty(xg), continue; end

    h = boxchart(g * ones(sum(idx), 1), xg, ...
        'BoxFaceColor', colors(g, :), ...
        'LineWidth', 1.2, ...
        'BoxWidth', box_width, ...
        'MarkerStyle', 'none');
    box_handles(g) = h;

    % One-sample t-test vs 0
    [~, p] = ttest(xg, 0);
    p_vals(g) = p;
    p_vals(g) = p;
    w_med(g)  = median(xg);
    w_std(g)  = std(xg);

    % Outliers (hollow circles): Tukey 1.5×IQR
    q1  = prctile(xg, 25);
    q3  = prctile(xg, 75);
    iqr = q3 - q1;
    lo  = q1 - 1.5 * iqr;
    hi  = q3 + 1.5 * iqr;

    is_out = (xg < lo) | (xg > hi);  % element-wise logical
    if any(is_out)
        yo = xg(is_out);
        xo = g + (rand(size(yo)) - 0.5) * 0.08;
        scatter(xo, yo, 24, 'o', ...
            'MarkerEdgeColor', colors(g, :), ...
            'MarkerFaceColor', 'none', ...
            'LineWidth', 1.1);
    end
end

% ===== Y-axis limits (robust to extremes) =====
all_vals = w_scores(~isnan(w_scores));
if ~isempty(all_vals)
    q_low   = prctile(all_vals, outlier_pct_lim(1));
    q_high  = prctile(all_vals, outlier_pct_lim(2));
    y_range = max(eps, q_high - q_low);
    ylim([q_low - y_range * y_margin_ratio, q_high + y_range * y_margin_ratio]);
end

% ===== BH-FDR (requires Bioinformatics Toolbox: mafdr) =====
q_vals = nan(size(p_vals));
q_vals = mafdr(p_vals, 'BHFDR', true);

% ===== Significance stars (FDR-first) =====
for g = 1:Ngroup
    p = p_vals(g);
    q = q_vals(g);
    if isnan(p), continue; end

    if ~isnan(q) && q < 0.01
        stars = '**';
    elseif ~isnan(q) && q < 0.05
        stars = '*';
    else
        stars = '';
    end

    if isempty(stars), continue; end
    y_star = w_med(g) + star_offset;
    text(g, y_star, stars, ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 11, 'FontWeight', 'bold');
end

% ===== Axes / title =====
xlim([0.5, Ngroup + 0.5]);
xticks(1:Ngroup);
xticklabels(varNames);
ylabel('W-score', 'FontSize', 18); % Increased font size for Y-label
set(gca, 'FontSize', 18);          % Increased font size for axis ticks (X and Y)

if exist('label_str', 'var') && ~isempty(label_str)
    title(sprintf('W-score by group: %s', label_str), 'Interpreter', 'none', 'FontSize', 14);
else
    title('W-score by group', 'FontSize', 14);
end
box on;

% ===== Export =====
if nargin >= 4 && ~isempty(output_name)
    print(gcf, output_name, '-dtiff', '-r600');
end

% ===== Outputs =====
w_vec   = w_med(:);
p_vec   = p_vals(:);
labels  = varNames(:);
std_vec = w_std(:);
end
