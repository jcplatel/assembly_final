function fig = figure_sequences(dossier)
%FIGURE_SEQUENCES  Étape 5, toutes les sessions, sans rien recalculer : lit
%   <PathSave>\sequences_agregation\sequences_agregation.mat.
%     A. chaque session, rangée par animal : rho de Spearman entre le
%        décalage des cellules de la séquence sur une moitié des courses et
%        sur l'autre (> 0 : même ordre) ; rond plein : p < 0,05 ; trait
%        noir : moyenne de l'animal.
%     B. tous les animaux : moyenne +/- sem ; p = Wilcoxon unilatéral sur
%        les animaux.
%     C. part des cellules dans la séquence ; D. nombre de courses ;
%     E. cellules de la séquence retrouvées sur l'autre moitié, contre le
%        nombre attendu si les deux choix étaient indépendants.
%   Écrit sequences.png.

dossier = regexprep(char(dossier), '[\\/]+$', '');
S = load(fullfile(dossier, 'sequences_agregation.mat'));
T = S.T(isfinite(S.T.rho), :);
b = S.bilan;
GRIS = [0.45 0.45 0.45];
animaux = unique(T.animal);
nA = numel(animaux);
[~, ia] = ismember(T.animal, animaux);
rng(1);
fig = figure('Color', 'w', 'Position', [40 40 1500 900]);
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Sequences during running (%d animals, %d sessions)', nA, height(T)), ...
    'FontWeight', 'bold', 'FontSize', 12);

% A. chaque session
ax = nexttile(tl, [1 2]); hold(ax, 'on');
for a = 1:nA
    m = find(ia == a);
    x = a + 0.16 * (rand(numel(m), 1) - 0.5) * 2;
    sig = T.p_rho(m) < 0.05;
    scatter(ax, x(~sig), T.rho(m(~sig)), 40, GRIS, 'LineWidth', 1);
    scatter(ax, x(sig), T.rho(m(sig)), 40, GRIS, 'filled');
    plot(ax, a + [-0.25 0.25], mean(T.rho(m)) * [1 1], 'k-', 'LineWidth', 2);
end
yline(ax, 0, '--', 'Color', [0.6 0.6 0.6]);
set(ax, 'XTick', 1:nA, 'XTickLabel', animaux, 'XLim', [0.4 nA + 0.6], 'TickLabelInterpreter', 'none', ...
    'FontSize', 9, 'TickDir', 'out');
xtickangle(ax, 30);
ylabel(ax, 'rho (order, one half vs other half of the runs)');
title(ax, {'A. Each session, by animal', sprintf('%d/%d sessions p < 0.05', b.n_sig, height(T))}, 'FontSize', 10);
box(ax, 'off');
yl = ylim(ax);

% B. tous les animaux
ax = nexttile(tl); hold(ax, 'on');
x = 1 + 0.10 * (rand(height(S.A), 1) - 0.5) * 2;
scatter(ax, x, S.A.rho, 40, GRIS, 'filled');
errorbar(ax, 1, b.rho_moyen, b.rho_sem, 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 7, ...
    'LineWidth', 2, 'CapSize', 8);
yline(ax, 0, '--', 'Color', [0.6 0.6 0.6]);
set(ax, 'XTick', 1, 'XTickLabel', {'all animals'}, 'XLim', [0.5 1.5], 'YLim', yl, 'FontSize', 9, 'TickDir', 'out');
title(ax, {'B. All animals', texte_p(b.p_wilcoxon)}, 'FontSize', 10);
box(ax, 'off');

% C, D : part des cellules dans la séquence, nombre de courses
valeurs = {T.pct_sequence, 'cells in the sequence (%)', 'C. Sequence cells'; ...
    T.n_courses, 'runs in the session', 'D. Runs'};
for k = 1:2
    ax = nexttile(tl); hold(ax, 'on');
    v = valeurs{k, 1};
    for a = 1:nA
        m = find(ia == a);
        x = a + 0.16 * (rand(numel(m), 1) - 0.5) * 2;
        scatter(ax, x, v(m), 36, GRIS, 'filled');
        plot(ax, a + [-0.25 0.25], mean(v(m)) * [1 1], 'k-', 'LineWidth', 2);
    end
    set(ax, 'XTick', 1:nA, 'XTickLabel', animaux, 'XLim', [0.4 nA + 0.6], 'TickLabelInterpreter', 'none', ...
        'FontSize', 9, 'TickDir', 'out');
    xtickangle(ax, 30);
    ylim(ax, [0, max(v) * 1.1 + eps]);
    ylabel(ax, valeurs{k, 2});
    title(ax, valeurs{k, 3}, 'FontSize', 10);
    box(ax, 'off');
end

% E. retrouvées sur l'autre moitié contre hasard
ax = nexttile(tl); hold(ax, 'on');
scatter(ax, T.n_retrouvees_hasard, T.n_retrouvees, 36, GRIS, 'filled');
lim = [0, max([T.n_retrouvees; T.n_retrouvees_hasard; 1]) * 1.1];
plot(ax, lim, lim, '--', 'Color', [0.6 0.6 0.6]);
set(ax, 'XLim', lim, 'YLim', lim, 'FontSize', 9, 'TickDir', 'out');
axis(ax, 'square');
xlabel(ax, 'expected if independent');
ylabel(ax, 'found again on the other half');
title(ax, 'E. Sequence cells found again', 'FontSize', 10);
box(ax, 'off');
exportgraphics(fig, fullfile(dossier, 'sequences.png'), 'Resolution', 200);
end


function t = texte_p(p)
if isnan(p)
    t = 'p = NA';
elseif p < 0.001
    t = 'p < 0.001';
else
    t = sprintf('p = %.3f', p);
end
end
