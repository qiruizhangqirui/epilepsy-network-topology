function local_stacked_bar(groupCat, x, varName, savePath, codes, labels, colors, figWidth, yLimMax)
% local_stacked_bar
% colors: [nCategory x 3] RGB, rows correspond to codes/labels
% figWidth: (optional) width of the figure in pixels. Default 520.
% yLimMax: (optional) maximum y-axis value in percentage. Default 100.

if nargin < 9 || isempty(yLimMax)
    yLimMax = 100;
end
if nargin < 8 || isempty(figWidth)
    figWidth = 520;
end


% -------------------------
% Prepare data
% -------------------------
idx = ~isundefined(groupCat) & ~isnan(x);
g = groupCat(idx);
x = x(idx);

G = categories(g);
M = zeros(numel(G), numel(codes));

for i = 1:numel(G)
    xi = x(g == G{i});
    for j = 1:numel(codes)
        M(i,j) = sum(xi == codes(j));
    end
end

% Convert to percentage
P = M ./ max(sum(M,2),1) * 100;

% -------------------------
% Plot
% -------------------------
fig = figure('Color','w','Units','pixels','Position',[100 100 figWidth 420]);

bh = bar(P, 'stacked');
set(bh, 'EdgeColor', 'none');

% Apply colors
for j = 1:numel(bh)
    bh(j).FaceColor = 'flat';
    bh(j).CData     = repmat(colors(j,:), size(P,1), 1);
end

% Axes settings
ax = gca;
ax.FontSize = 14;
ax.LineWidth = 1.2;
ax.TickDir = 'out';

xticks(1:numel(G));
xticklabels(G);
xlim([0.5, numel(G)+0.5]);

ylim([0 yLimMax]);

% Smart yticks
if yLimMax <= 25
    yticks(0:5:yLimMax);
    yticklabels(string(0:5:yLimMax) + "%");
elseif yLimMax <= 50
    yticks(0:10:yLimMax);
    yticklabels(string(0:10:yLimMax) + "%");
else
    yticks(0:25:yLimMax);
    yticklabels(string(0:25:yLimMax) + "%");
end

ylabel('');  % remove "Proportion"
title(varName, 'Interpreter','none', 'FontSize', 16, 'FontWeight','bold');

% Legend
lgd = legend(labels, 'Location','eastoutside');
lgd.FontSize = 13;
lgd.Box = 'off';

% Tight layout feel
box off
set(ax, 'LooseInset', max(get(ax,'TightInset'), 0.02))

% Save
exportgraphics(fig, savePath, 'Resolution', 300);


end
