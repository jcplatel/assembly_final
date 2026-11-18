function fig = figure_coeurs(dossier)
%FIGURE_COEURS  Étape 3d : cœur clonal de chaque assemblée, sans rien
%   recalculer : lit coeurs_agregation.mat dans <PathSave>\soeurs_agregation\.
%     A. score de chaque assemblée testable (z contre des fausses
%        assemblées de voisines), par animal ; rouge = cœur clonal (FDR),
%        orange = contient des sœurs seulement, gris = pas de preuve ;
%     B. les deux p de chaque assemblée (-log10) ; pointillés = 0,05 ;
%     C. par animal : part des assemblées testables avec un cœur clonal.
%   Écrit coeurs_assemblees.png.

dossier = regexprep(char(dossier), '[\\/]+$', '');
S = load(fullfile(dossier, 'coeurs_agregation.mat'));
C = S.C(S.C.testable, :);
P = S.P;
coul = struct('gris', [0.65 0.65 0.65], 'orange', [0.93 0.55 0.15], 'rouge', [0.80 0.10 0.10]);
v = C.verdict;
cc = repmat(coul.gris, height(C), 1);
cc(strcmp(v, 'contient des sœurs'), :) = repmat(coul.orange, sum(strcmp(v, 'contient des sœurs')), 1);
cc(strcmp(v, 'cœur clonal'), :) = repmat(coul.rouge, sum(strcmp(v, 'cœur clonal')), 1);
animaux = P.animal;
nA = numel(animaux);
[~, ia] = ismember(C.animal, animaux);
rng(1);
fig = figure('Color', 'w', 'Position', [40 60 1600 520]);
tl = tiledlayout(fig, 1, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Clonal core of each assembly (%d testable assemblies, %d animals)', height(C), nA), ...
    'FontWeight', 'bold', 'FontSize', 12);

% A. score par animal
ax = nexttile(tl, [1 2]); hold(ax, 'on');
x = ia + 0.3 * (rand(height(C), 1) - 0.5);
scatter(ax, x, C.z_voisines, 36, cc, 'filled', 'MarkerEdgeColor', [0.3 0.3 0.3]);
yline(ax, 0, '--', 'Color', [0.6 0.6 0.6]);
set(ax, 'XTick', 1:nA, 'XTickLabel', animaux, 'XLim', [0.4 nA + 0.6], 'TickLabelInterpreter', 'none', ...
    'FontSize', 9, 'TickDir', 'out');
xtickangle(ax, 30);
ylabel(ax, 'clonal core score (z vs neighbour cells)');
title(ax, 'A. Each assembly', 'FontSize', 10);
box(ax, 'off');

% B. les deux p
ax = nexttile(tl); hold(ax, 'on');
scatter(ax, -log10(C.p_champ), -log10(C.p_voisines), 30, cc, 'filled', 'MarkerEdgeColor', [0.3 0.3 0.3]);
xline(ax, -log10(0.05), ':', 'Color', [0.5 0.5 0.5]);
yline(ax, -log10(0.05), ':', 'Color', [0.5 0.5 0.5]);
xlabel(ax, '-log_{10} p (shuffled colours)');
ylabel(ax, '-log_{10} p (neighbour cells)');
axis(ax, 'square');
set(ax, 'FontSize', 9, 'TickDir', 'out');
title(ax, 'B. Two tests', 'FontSize', 10);
box(ax, 'off');

% C. part par animal
ax = nexttile(tl, [1 2]); hold(ax, 'on');
tous = 100 * P.n_coeur ./ max(1, P.n_testables);
b1 = bar(ax, 1:nA, tous, 0.6, 'FaceColor', coul.rouge, 'EdgeColor', 'none');
for i = 1:nA
    text(ax, i, tous(i), sprintf('%d/%d', P.n_coeur(i), P.n_testables(i)), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 8);
end
legend(ax, b1, {'clonal core (more sisters than neighbour cells)'}, 'Location', 'northeast', ...
    'Box', 'off', 'FontSize', 8);
set(ax, 'XTick', 1:nA, 'XTickLabel', animaux, 'XLim', [0.4 nA + 0.6], 'TickLabelInterpreter', 'none', ...
    'FontSize', 9, 'TickDir', 'out');
xtickangle(ax, 30);
ylim(ax, [0, max([tous; 10]) * 1.3]);
ylabel(ax, '% of testable assemblies');
title(ax, sprintf('C. Per animal (FDR %g%%)', 100 * S.opts.fdr_q), 'FontSize', 10);
box(ax, 'off');
exportgraphics(fig, fullfile(dossier, 'coeurs_assemblees.png'), 'Resolution', 200);
end
