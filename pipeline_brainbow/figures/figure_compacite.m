function fig = figure_compacite(dossier)
%FIGURE_COMPACITE  Compacité spatiale des assemblées (étape 4, sans couleur),
%   sans rien recalculer : lit compacite_assemblees.mat.
%
%   figure_compacite('E:\...\agregation\controles_champ')
%
%   A. chaque assemblée : étalement observé (distance moyenne entre ses
%      cellules, µm) contre l'étalement moyen de groupes de même taille tirés
%      au hasard parmi les cellules actives ; sous la diagonale = plus
%      compacte. Rond plein : plus compacte que 95 % des tirages.
%   B. rapport observé / hasard, un point par assemblée rangé par animal
%      (plein : p <= 0,05), trait noir = l'animal (moyenne).
%   C. répartition des p des assemblées : plates si les assemblées ne sont
%      pas plus compactes que le hasard ; pic à gauche si elles le sont.
%   Écrit assembly_compactness.png.

dossier = regexprep(char(dossier), '[\\/]+$', '');
C = load(fullfile(dossier, 'compacite_assemblees.mat'));
T = C.T;
animaux = C.animaux;
nA = numel(animaux);
[~, ia] = ismember(T.animal, animaux);
couleurs_animaux = lines(max(nA, 1));
sig = T.p <= 0.05;

fig = figure('Color', 'w', 'Position', [60 60 1600 560]);
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['Are assemblies spatially compact? (no colour used; %d assemblies, %d animals) | ' ...
    'ratio %.2f [%.2f-%.2f], mixed model p = %s | Wilcoxon on animals p = %s'], height(T), nA, ...
    C.rapport_lmm, C.rapport_ic(1), C.rapport_ic(2), texte_p(C.p_lmm), texte_p(C.p_wilcoxon)), ...
    'FontWeight', 'bold', 'FontSize', 11);

% A. observé contre hasard
ax = nexttile(tl); hold(ax, 'on');
h = gobjects(nA, 1);
for a = 1:nA
    m = ia == a;
    h(a) = scatter(ax, T.etalement_hasard_um(m & ~sig), T.etalement_um(m & ~sig), 30, couleurs_animaux(a, :), ...
        'LineWidth', 1);
    scatter(ax, T.etalement_hasard_um(m & sig), T.etalement_um(m & sig), 30, couleurs_animaux(a, :), 'filled');
end
lim = [0, max([T.etalement_um; T.etalement_hasard_um]) * 1.05];
plot(ax, lim, lim, 'k--');
axis(ax, 'square'); xlim(ax, lim); ylim(ax, lim);
xlabel(ax, 'random groups of the same size: mean distance (\mum)');
ylabel(ax, 'assembly: mean distance between its cells (\mum)');
legend(ax, h, animaux, 'Location', 'southeast', 'Box', 'off', 'FontSize', 7, 'Interpreter', 'none');
title(ax, {'A. Each assembly (below the line = more compact)', ...
    sprintf('filled: more compact than 95%% of random groups (%d/%d)', sum(sig), height(T))}, 'FontSize', 10);
box(ax, 'off');

% B. rapport par animal
ax = nexttile(tl); hold(ax, 'on');
rng(1);
for a = 1:nA
    m = ia == a;
    x = a + 0.16 * (rand(sum(m), 1) - 0.5) * 2;
    r = T.rapport(m);
    s = sig(m);
    scatter(ax, x(~s), r(~s), 30, [0.45 0.45 0.45], 'LineWidth', 1);
    scatter(ax, x(s), r(s), 30, [0.45 0.45 0.45], 'filled');
    plot(ax, a + [-0.25 0.25], (1 + C.ecart_a(a)) * [1 1], 'k-', 'LineWidth', 2);
end
yline(ax, 1, '--', 'Color', [0.6 0.6 0.6]);
set(ax, 'XTick', 1:nA, 'XTickLabel', animaux, 'XLim', [0.4 nA + 0.6], 'TickLabelInterpreter', 'none', ...
    'FontSize', 9, 'TickDir', 'out');
xtickangle(ax, 30);
ylabel(ax, 'spatial extent: assembly / random groups');
title(ax, {'B. Ratio < 1: more compact than chance', ...
    sprintf('dot = assembly (filled: p <= 0.05); line = animal | Wilcoxon p = %s', texte_p(C.p_wilcoxon))}, ...
    'FontSize', 10);
box(ax, 'off');

% C. répartition des p
ax = nexttile(tl); hold(ax, 'on');
bords = 0:0.05:1;
n = histcounts(min(T.p, 1 - 1e-9), bords);
bar(ax, bords(1:end - 1) + 0.025, n, 1, 'FaceColor', [0.55 0.55 0.55], 'EdgeColor', 'w');
yline(ax, height(T) / (numel(bords) - 1), 'k--', 'expected if not compact', 'LabelHorizontalAlignment', 'right', ...
    'FontSize', 8);
xlabel(ax, 'p of each assembly (fraction of random groups as compact)');
ylabel(ax, 'assemblies');
xlim(ax, [0 1]);
title(ax, {'C. Distribution of p', 'flat = not compact; peak on the left = compact'}, 'FontSize', 10);
box(ax, 'off');

exportgraphics(fig, fullfile(dossier, 'assembly_compactness.png'), 'Resolution', 200);
end


function t = texte_p(p)
if isnan(p)
    t = 'NA';
elseif p < 0.001
    t = '< 0.001';
else
    t = sprintf('%.3f', p);
end
end
