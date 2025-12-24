function [pval_matrix, pmaxT_matrix] = plot_group_wscore_by_network( ...
    Data, GroupTable, Atlas_network_name, output_name)
% 绘制每组在每个网络的 W-score 箱线图，并进行：
%   1) 单样本 t 检验 vs 0 → 原始 p
%   2) 族内 max-T 置换校正（Tmax）→ p_perm
%
% 输入:
%   - Data: Nsub × Nnet 的 W-score 矩阵
%   - GroupTable: Nsub × Ng 的分组表（列为任意组名；>0 视作属于该组）
%   - Atlas_network_name: 1 × Nnet 的 cell，网络名称
%   - output_name: 输出图像文件路径（建议含扩展名；若不含，本函数以TIFF导出）
%
% 输出:
%   - pval_matrix : Ng_valid × Nnet 的原始 t 检验 p 值
%   - pmaxT_matrix: Ng_valid × Nnet 的 Tmax 校正后 p 值
%
% 备注：
%   - 显著性标注采用：p<0.05(*)；Tmax<0.05(**)；Tmax<0.01(***)
%   - Tmax 置换参数：alpha_maxT=0.05, nperm=10000, 双尾，族内校正
%
% 依赖：
%   - permuztest(G, 0, sigma, 'nperm', 10000, 'alpha', 0.05, 'tail','both','correct',true)

% ---------------- 配置 ----------------
remove_outliers = true;     % y 轴范围是否做稳健裁剪（不影响统计）
outlier_pct     = [0.1, 99.9];

alpha_maxT = 0.05;
nperm      = 10000;

% ---------------- 检查输入 ----------------
[Nsub, Nnet] = size(Data);
if height(GroupTable) ~= Nsub
    error('Data 的行数 (%d) 与 GroupTable 的行数 (%d) 不一致。', Nsub, height(GroupTable));
end
if numel(Atlas_network_name) ~= Nnet
    error('Atlas_network_name 的长度 (%d) 必须等于 Data 的列数 (%d)。', numel(Atlas_network_name), Nnet);
end

% ---------------- 选择有效组列（非全0/NaN）----------------
all_names   = GroupTable.Properties.VariableNames;
Ng_all      = numel(all_names);
valid_mask  = false(1, Ng_all);
for j = 1:Ng_all
    col = GroupTable{:, j};
    if ~isnumeric(col) && ~islogical(col)
        error('GroupTable 的列 "%s" 不是数值或逻辑类型。', all_names{j});
    end
    m = (col > 0) & isfinite(col);
    valid_mask(j) = any(m);
end

if ~any(valid_mask)
    error('GroupTable 中没有含有成员的组（所有列均为全0/NaN）。');
end

Group_labels = all_names(valid_mask);
GroupMat     = GroupTable{:, valid_mask};
Ngroup       = numel(Group_labels);

% ---------------- 初始化 ----------------
colors         = lines(max(Ngroup, 5));
pval_matrix    = nan(Ngroup, Nnet);   % 每组 vs 0 的原始 p
pmaxT_matrix   = nan(Ngroup, Nnet);   % 每组 × 网络的 Tmax 校正 p
group_medians  = nan(Ngroup, Nnet);   % 中位数记录
group_data     = cell(Ngroup,1);      % 存每组的 n×Nnet 数据，供 Tmax 使用

% x 位置设置
x_spacing   = 1.2;
x_positions = (0:(Nnet-1)) * x_spacing;
x_shift     = linspace(-0.4, 0.4, Ngroup);  % group 偏移

% ---------------- 画布 ----------------
figure('Units','normalized','Position',[0.05 0.1 0.9 0.7]); hold on;
box_width   = 0.15;
star_offset = 0.2;

% 灰色参考线 y = 0
plot([min(x_positions)-1, max(x_positions)+1], [0 0], ...
    'Color', [0.7 0.7 0.7], 'LineStyle', '-', 'LineWidth', 1.2);

% 初始化 boxchart 句柄（用于 legend）
box_handles = gobjects(Ngroup,1);

% ---------------- 绘制（并做 t 检验）----------------
for i = 1:Nnet
    for g = 1:Ngroup
        idx = GroupMat(:, g) > 0;           % >0 视为属于该组
        this_data = Data(idx, i);
        if isempty(this_data), continue; end

        % 记录该组的全网络数据（用于 Tmax）
        if i == 1
            group_data{g} = Data(idx, :);   % n × Nnet
        end

        % 绘图
        x_pos = x_positions(i) + x_shift(g);
        h = boxchart(x_pos * ones(sum(idx),1), this_data, ...
            'BoxFaceColor', colors(g,:), ...
            'LineWidth', 1.2, ...
            'BoxWidth', box_width, ...
            'MarkerStyle', 'none');

        if i == 1
            box_handles(g) = h;  % 首个网络记录句柄用于 legend
        end

        % 单样本 t 检验 vs 0（原始 p）
        [~, p] = ttest(this_data, 0);
        pval_matrix(g, i) = p;

        % 中位数（用于星标放置）
        group_medians(g, i) = median(this_data);
    end
end

% ---------------- Tmax（max-T）置换校正（逐组，一次性覆盖所有网络）----------------
for g = 1:Ngroup
    Xg = group_data{g};          % n × Nnet
    if isempty(Xg), continue; end

    % 列标准差作为 sigma，防止 0 与非有限
    sigma = std(Xg, 0, 1);
    sigma(sigma == 0 | ~isfinite(sigma)) = eps;

    try
        % 双尾，族内校正
        % 返回 p_perm（1×Nnet）：在家族内经 max-T 校正后的 p 值
        [~, p_perm, ~, ~] = permuztest(Xg, 0, sigma, ...
            'nperm', nperm, 'alpha', alpha_maxT, ...
            'tail', 'both', 'correct', true);
    catch ME
        warning('组 %s 的 Tmax 置换失败：%s。请确认 permuztest 在路径上。', Group_labels{g}, ME.message);
        p_perm = nan(1, size(Xg,2));
    end

    pmaxT_matrix(g, :) = p_perm;
end

% ---------------- 显著性标注（*, **, ***）----------------
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

% ---------------- Y 轴范围设置（仅影响显示）----------------
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

% ---------------- 坐标轴与标签 ----------------
xticks(x_positions);
xticklabels(Atlas_network_name);
xtickangle(45);
ylabel('W-score');
legend(box_handles, Group_labels, 'Location', 'best');
title('Group W-score per Network (t-test & max-T)');
xlim([min(x_positions)-0.8, max(x_positions)+0.8]);
box on;

% ---------------- 保存图像 ----------------
if nargin >= 4 && ~isempty(output_name)
    [~,~,ext] = fileparts(output_name);
    if isempty(ext)
        print(gcf, output_name, '-dtiff', '-r600');
    else
        print(gcf, output_name, '-r600', ['-d' strrep(ext,'.','')]);
    end
end
% close(gcf);

end
