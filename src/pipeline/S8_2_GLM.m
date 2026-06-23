% =========================================================================
% Script: S8_2_GLM_noASM.m
% -------------------------------------------------------------------------
% Description:
%   GLM-based analysis pipeline between Network Features (X) and Clinical (Y).
%   Performs:
%     1. GLM (linear regression) for each pair:
%        X_Feature ~ 1 + Y_Variable + Age + Sex
%     2. Benjamini-Hochberg FDR correction within each clinical variable (Y).
%     3. Stat summary table export to CSV.
%     4. Generates high-quality heatmaps of T-values (T-statistic)
%        separately for FDR-corrected significance.
%
% Inputs:
%   - data/Network_Wscore.xlsx (Sheet 1)
%
% Outputs:
%   - outputs/GLM_noASM/GLM_stats_summary.csv
%   - outputs/GLM_noASM/Heatmap_<Subset>_FDR.png (showing FDR sig)
%
% Author: Antigravity
% Date: 06/18/2026
% =========================================================================

clear;
clc;
set(0, 'DefaultFigureVisible', 'off'); % Run in headless plotting mode

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
out_dir     = fullfile(project_root, 'outputs', 'GLM');
xlsx_file   = fullfile(data_dir, 'Network_Wscore.xlsx');

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

fprintf('Project Root: %s\n', project_root);
fprintf('Loading Excel file: %s (Sheet 1)\n', xlsx_file);

if ~exist(xlsx_file, 'file')
    error('Input file not found: %s. Please check paths or run S7 first.', xlsx_file);
end

% Read Sheet 1 (TJU data)
T = readtable(xlsx_file, 'Sheet', 1, 'VariableNamingRule', 'preserve');

%% =========================
% Part 2: Variable Definitions
% =========================
all_varnames = T.Properties.VariableNames;

% X variables (Imaging)
correspondence_mask = strncmp('Correspondence_', all_varnames, 15);
x1_cols = [{'Normativity', 'Non normativity'}, all_varnames(correspondence_mask)];  % Network Correspondence
hubness_mask = strncmp('Hubness_', all_varnames, 8);
x2_cols = all_varnames(hubness_mask);                                               % Hubness

% Y variables (Clinical)
y1_cols = {'IQ','TMTA','TMTB','CVLT TL','CVLT LDFR','BNT','Letter Fluency',...
           'Semantic fluency','Vocabulary','Similarities','Matrix Reasoning',...
           'Digit Span','Block Design','Coding','WCST PR','WCST CC',...
           'Logical Memory 1','Logical Memory 2','ROCF Copy','Pegboard R','Pegboard L'}; % Cognition

y2_cols = {'EpilepsyType_TLE_EXE','FTBTC_Yes_No','Lateralization_Left','Lateralization_Right','Lateralization_Unclear',...
           'Pathology_UHS','Pathology_BHS','Pathology_Lesion','Pathology_Normal', 'AgeOnset', 'SeizureDuration'}; % Epilepsy Clinicals (ASM_count removed)

% Validate columns
x_all = unique([x1_cols, x2_cols], 'stable');
y_all = unique([y1_cols, y2_cols], 'stable');
missing_x = setdiff(x_all, T.Properties.VariableNames);
missing_y = setdiff(y_all, T.Properties.VariableNames);

if ~isempty(missing_x) || ~isempty(missing_y)
    if ~isempty(missing_x)
        fprintf('\nError: The following X features are missing from the excel file:\n');
        disp(missing_x');
    end
    if ~isempty(missing_y)
        fprintf('\nError: The following Y clinical variables are missing from the excel file:\n');
        disp(missing_y');
    end
    error('Excel columns validation failed. Please make sure the Excel file has been updated by running S7.');
end

% Define 4 subsets
subsets = { ...
    struct('name', 'Correspondence_vs_Cognition', 'x', {x1_cols}, 'y', {y1_cols}), ...
    struct('name', 'Correspondence_vs_Clinicals', 'x', {x1_cols}, 'y', {y2_cols}), ...
    struct('name', 'Hubness_vs_Cognition',        'x', {x2_cols}, 'y', {y1_cols}), ...
    struct('name', 'Hubness_vs_Clinicals',        'x', {x2_cols}, 'y', {y2_cols})  ...
};

all_stats_cells = cell(length(subsets), 1);

%% =========================
% Part 3: Subset Analyses Loop
% =========================
for s = 1:length(subsets)
    sub = subsets{s};
    fprintf('\n------------------------------------------------------------\n');
    fprintf('Processing Subset: %s (%d X, %d Y)...\n', sub.name, length(sub.x), length(sub.y));
    fprintf('------------------------------------------------------------\n');
    
    num_x_sub = length(sub.x);
    num_y_sub = length(sub.y);
    sub_results = cell(num_x_sub * num_y_sub, 10);
    sub_row = 1;
    
    for i = 1:num_x_sub
        x_col = sub.x{i};
        x_val = double(T.(x_col));
        for j = 1:num_y_sub
            y_col = sub.y{j};
            y_val = double(T.(y_col));
            
            % Covariates
            age_val = double(T.Age);
            sex_val = double(T.Sex);
            
            % Check complete cases
            cc_mask = ~isnan(x_val) & ~isnan(y_val) & ~isnan(age_val) & ~isnan(sex_val);
            n_used = sum(cc_mask);
            
            if n_used < 10
                sub_results{sub_row, 1} = x_col;
                sub_results{sub_row, 2} = y_col;
                sub_results{sub_row, 3} = 'GLM';
                sub_results{sub_row, 4} = NaN;
                sub_results{sub_row, 5} = NaN;
                sub_results{sub_row, 6} = NaN;
                sub_results{sub_row, 7} = n_used;
                sub_results{sub_row, 8} = NaN;
                sub_results{sub_row, 9} = NaN;
                sub_results{sub_row, 10} = NaN;
                sub_row = sub_row + 1;
                continue;
            end
            
            tbl = table(x_val(cc_mask), y_val(cc_mask), age_val(cc_mask), sex_val(cc_mask), ...
                'VariableNames', {'X', 'Y', 'Age', 'Sex'});
            formula = 'X ~ 1 + Y + Age + Sex';
            
            % Fit GLM
            try
                mdl = fitlm(tbl, formula);
                coeff_names = mdl.CoefficientNames;
                y_idx = find(strcmp(coeff_names, 'Y'));
                
                if ~isempty(y_idx)
                    beta_val = mdl.Coefficients.Estimate(y_idx);
                    t_val = mdl.Coefficients.tStat(y_idx);
                    p_val = mdl.Coefficients.pValue(y_idx);
                    df_val = mdl.DFE;
                    r_effect = t_val / sqrt(t_val^2 + df_val);
                else
                    beta_val = NaN;
                    t_val = NaN;
                    p_val = NaN;
                    r_effect = NaN;
                end
            catch
                beta_val = NaN;
                t_val = NaN;
                p_val = NaN;
                r_effect = NaN;
            end
            
            sub_results{sub_row, 1} = x_col;
            sub_results{sub_row, 2} = y_col;
            sub_results{sub_row, 3} = 'GLM';
            sub_results{sub_row, 4} = t_val;      % Statistic_Value
            sub_results{sub_row, 5} = beta_val;   % Raw_Effect
            sub_results{sub_row, 6} = r_effect;   % Effect_Size_r (still saved in table for reference)
            sub_results{sub_row, 7} = n_used;
            sub_results{sub_row, 8} = p_val;
            sub_results{sub_row, 9} = NaN;        % Placeholder for FDR_P_Value
            sub_results{sub_row, 10} = NaN;       % Placeholder for MaxT_P_Value (not run, filled with NaN for compatibility)
            sub_row = sub_row + 1;
        end
    end
    
    sub_table = cell2table(sub_results, 'VariableNames', ...
        {'X_Feature', 'Y_Variable', 'Test_Type', 'Statistic_Value', 'Raw_Effect', 'Effect_Size_r', 'N', 'P_Value', 'FDR_P_Value', 'MaxT_P_Value'});
    
    % 2. Perform FDR Correction separately for each Y variable (cognitive/clinical parameter) within this subset
    p_vals = sub_table.P_Value;
    fdr_p_all = nan(size(p_vals));
    
    for j = 1:num_y_sub
        y_col = sub.y{j};
        y_mask = strcmp(sub_table.Y_Variable, y_col);
        
        p_to_correct = sub_table.P_Value(y_mask);
        non_nan_mask = ~isnan(p_to_correct);
        p_to_correct_valid = p_to_correct(non_nan_mask);
        
        if ~isempty(p_to_correct_valid)
            % Benjamini-Hochberg FDR
            [sorted_p, sort_idx] = sort(p_to_correct_valid);
            m_correct = length(p_to_correct_valid);
            adj_p_sorted = sorted_p .* (m_correct ./ (1:m_correct)');
            for k = m_correct-1:-1:1
                adj_p_sorted(k) = min(adj_p_sorted(k), adj_p_sorted(k+1));
            end
            adj_p_sorted = min(adj_p_sorted, 1);
            adj_p = zeros(size(p_to_correct_valid));
            adj_p(sort_idx) = adj_p_sorted;
            
            fdr_p_y = nan(size(p_to_correct));
            fdr_p_y(non_nan_mask) = adj_p;
            fdr_p_all(y_mask) = fdr_p_y;
        end
    end
    sub_table.FDR_P_Value = fdr_p_all;
    
    % Add Subset Identifier
    sub_table.Subset_Name = repmat({sub.name}, height(sub_table), 1);
    all_stats_cells{s} = sub_table;
    
    % 3. Extract matrices for heatmap plotting (using Statistic_Value/T-value)
    [t_mat, raw_p_mat, fdr_p_mat] = extract_submatrix_t(sub_table, sub.x, sub.y, 'FDR_P_Value');
    
    % Define custom figure size based on dimensions
    fig_w = max(600, min(1400, 100 + 40 * length(sub.y)));
    fig_h = max(500, min(1400, 150 + 20 * length(sub.x)));
    fig_size = [fig_w, fig_h];
    
    % 4. Plot FDR Heatmap (showing T-values)
    fdr_title = sprintf('%s (GLM T-values, FDR corrected within variable)', strrep(sub.name, '_', ' '));
    fdr_save_path = fullfile(out_dir, sprintf('Heatmap_%s_FDR.png', sub.name));
    plot_heatmap_t(sub.x, sub.y, t_mat, raw_p_mat, fdr_p_mat, fdr_title, fdr_save_path, fig_size);
end

%% =========================
% Part 4: Export Combined Table
% =========================
res_table_all = vertcat(all_stats_cells{:});
summary_file = fullfile(out_dir, 'GLM_stats_summary.csv');
writetable(res_table_all, summary_file);
fprintf('\n============================================================\n');
fprintf('Exported combined GLM statistics summary to:\n  %s\n', summary_file);
fprintf('All Done! Heatmaps and statistics saved to: %s\n', out_dir);


%% =========================================================================
% Helper Functions
% =========================================================================

function [t_mat, p_mat, fdr_mat] = extract_submatrix_t(res_table, rows_subset, cols_subset, p_col_name)
    num_r = length(rows_subset);
    num_c = length(cols_subset);
    t_mat = nan(num_r, num_c);
    p_mat = nan(num_r, num_c);
    fdr_mat = nan(num_r, num_c);
    
    for r = 1:num_r
        r_name = rows_subset{r};
        for c = 1:num_c
            c_name = cols_subset{c};
            idx = strcmp(res_table.X_Feature, r_name) & strcmp(res_table.Y_Variable, c_name);
            if any(idx)
                t_mat(r, c) = res_table.Statistic_Value(idx); % T-value
                p_mat(r, c) = res_table.P_Value(idx);
                fdr_mat(r, c) = res_table.(p_col_name)(idx);
            end
        end
    end
end

function plot_heatmap_t(X_names, Y_names, data_matrix, p_matrix, corrected_p_matrix, title_str, save_path, fig_size)
    fig = figure('Position', [100, 100, fig_size(1), fig_size(2)], 'Visible', 'off');
    
    % Diverging colormap: Blue (negative) -> White (zero) -> Red (positive)
    c_map = [
        linspace(0, 1, 64)', linspace(0.4, 1, 64)', linspace(1, 1, 64)';   % Blue to white
        linspace(1, 1, 64)', linspace(1, 0.4, 64)', linspace(1, 0, 64)'    % White to red
    ];
    c_map(64, :) = []; % Remove duplicate white row in the center
    
    % Plot
    imagesc(data_matrix);
    colormap(c_map);
    
    % Determine color axis limits symmetrically
    valid_data = data_matrix(~isnan(data_matrix));
    if isempty(valid_data)
        max_val = 4;
    else
        max_val = max(abs(valid_data(:)));
    end
    max_val = max(max_val, 2); % Prevent division by zero or too small limits
    caxis([-max_val, max_val]);
    
    % Visual styling
    cb = colorbar;
    set(cb, 'FontSize', 12);
    ylabel(cb, 't-statistic (T-value)', 'FontSize', 13);
    
    % Axis labels
    set(gca, 'XTick', 1:length(Y_names), 'XTickLabel', Y_names);
    set(gca, 'YTick', 1:length(X_names), 'YTickLabel', X_names);
    xtickangle(45);
    set(gca, 'TickLabelInterpreter', 'none');
    set(gca, 'FontSize', 12);
    
    % Add text annotations (* for p < 0.05, ** for corrected p < 0.05)
    hold on;
    [num_rows, num_cols] = size(data_matrix);
    for r = 1:num_rows
        for c = 1:num_cols
            p_val = p_matrix(r, c);
            corr_val = corrected_p_matrix(r, c);
            
            if ~isnan(p_val)
                if corr_val < 0.05
                    text(c, r, '**', 'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'middle', 'Color', 'black', 'FontSize', 16, 'FontWeight', 'bold');
                elseif p_val < 0.05
                    text(c, r, '*', 'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'middle', 'Color', 'black', 'FontSize', 16, 'FontWeight', 'bold');
                end
            end
        end
    end
    hold off;
    
    % Print high-res image
    print(fig, save_path, '-dpng', '-r150');
    close(fig);
end
