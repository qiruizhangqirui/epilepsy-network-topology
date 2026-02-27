function local_plot_bar(vals_L, pvals_L, vals_R, pvals_R, suffix_name, out_dir)
    % Define Ordered Labels (S11 Order)
    % 17 Networks
    % Ord mapping from Provided Input Order to S11 preference
    % Input Order: 
    % 1.TempPar 2.DefC 3.DefB 4.DefA 5.ContC 6.ContB 7.ContA 8.LimA 9.LimB 
    % 10.SalVenB 11.SalVenA 12.DorsAttnB 13.DorsAttnA 14.SomMotB 15.SomMotA 16.VisPeri 17.VisCent
    
    ord = [17, 16, 15, 14, 12, 11, 13, 1, 10, 7, 8, 5, 9, 2, 6, 3, 4];
    
    % Full Labels (S11 Order)
    labels_s11 = {
        'VisualA', 'VisualB', 'SomatomotorA', 'SomatomotorB', ...
        'DorsAttnB', 'SalVenAttnA', 'DorsAttnA', 'TempPar', ...
        'SalVenAttnB', 'ControlA', 'LimbicA', 'ControlC', ...
        'LimbicB', 'DefaultC', 'ControlB', 'DefaultB', 'DefaultA'
    };
    
    % Function to get colors
    get_colors = @(v, p) get_bar_colors(v, p);
    
    % Reorder Data
    if numel(vals_L) ~= 17 || numel(vals_R) ~= 17
        warning('local_plot_bar_S11: Input must be 17 elements per hemisphere.');
        return;
    end
    
    vals_L_ord = vals_L(ord); pvals_L_ord = pvals_L(ord);
    vals_R_ord = vals_R(ord); pvals_R_ord = pvals_R(ord);
    
    colors_L = get_colors(vals_L_ord, pvals_L_ord);
    colors_R = get_colors(vals_R_ord, pvals_R_ord);

    % Determine Y Label
    if contains(suffix_name, 'Corr', 'IgnoreCase', true)
        ylab = 'r value';
    elseif contains(suffix_name, 'Ttest', 'IgnoreCase', true)
        ylab = 't value';
    else
        ylab = 'Statistic';
    end

    % Create Compound Figure
    f = figure('Color','w','Position',[100, 100, 1000, 630]); 
    
    % --- Top: Left Hemisphere ---
    subplot(2,1,1);
    set(gca, 'Position', [0.10, 0.60, 0.85, 0.35]); % Manual Position (Top)
    
    b1 = bar(1:17, vals_L_ord, 'FaceColor','flat');
    b1.CData = colors_L;
    b1.EdgeColor = 'none';
    
    apply_y_settings(gca, ylab);
    
    xticks(1:17);
    xticklabels({}); % No labels for top plot
    ylabel(ylab, 'FontSize', 15);
    % title('Left Hemisphere'); % Removed
    box off; 
    set(gca, 'TickDir','out', 'LineWidth',1.2, 'FontSize',13, 'Color','none');
    
    % --- Bottom: Right Hemisphere ---
    subplot(2,1,2);
    set(gca, 'Position', [0.10, 0.22, 0.85, 0.35]); % Manual Position (Bottom)
    
    b2 = bar(1:17, vals_R_ord, 'FaceColor','flat');
    b2.CData = colors_R;
    b2.EdgeColor = 'none';

    apply_y_settings(gca, ylab);
    
    xticks(1:17);
    xticklabels(labels_s11);
    xtickangle(45);
    ylabel(ylab, 'FontSize', 15);
    % title('Right Hemisphere'); % Removed
    box off; 
    set(gca, 'TickDir','out', 'LineWidth',1.2, 'FontSize',13, 'Color','none');

    % Save as Transparent PNG
    scan_name = fullfile(out_dir, [suffix_name '.png']);
    exportgraphics(f, scan_name, 'Resolution', 300, 'BackgroundColor', 'none');
    close(f);
end

function colors = get_bar_colors(vals, pvals)
    colors = repmat([0.8, 0.8, 0.8], 17, 1); % Default Gray
    
    % Blue: Significance Negative
    idx_neg = pvals < 0.05 & vals < 0;
    if any(idx_neg)
        colors(idx_neg, :) = repmat([0.2, 0.4, 0.8], sum(idx_neg), 1);
    end

    % Red: Significance Positive
    idx_pos = pvals < 0.05 & vals > 0;
    if any(idx_pos)
        colors(idx_pos, :) = repmat([0.8, 0.2, 0.2], sum(idx_pos), 1);
    end
end

function apply_y_settings(ax, ylab)
    % Fixed Y-Axis Limits and Ticks
    if strcmp(ylab, 'r value')
        ylim(ax, [-0.5, 0]); 
        yticks(ax, -0.5:0.1:0);
    elseif strcmp(ylab, 't value')
        ylim(ax, [-5, 2]);
        yticks(ax, -5:1:2);
    end
end
