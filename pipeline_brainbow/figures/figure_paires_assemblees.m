function fig = figure_paires_assemblees(dossier)
%FIGURE_PAIRES_ASSEMBLEES  Figure de main_paires_assemblees, sans rien
%   recalculer : lit paires_assemblees.mat dans dossier (10 nov.).
%     A. test principal (assemblées valides, mélange global) : excès de
%        paires de même couleur dans les assemblées, pondérées par la rareté
%        de la couleur ; un rond par session (plein : p < 0,05), trait =
%        moyenne de l'animal ;
%     B. contrôle spatial (mélange par blocs), même présentation ;
%     C. chaque session : observé contre attendu par hasard ;
%     D. les quatre analyses : moyenne des animaux +- sem et p (Wilcoxon).
%   Écrit paires_assemblees.png dans dossier.

dossier = regexprep(char(dossier), '[\\/]+$', '');
S = load(fullfile(dossier, 'paires_assemblees.mat'));
T = S.T; A = S.A; B = S.B;
fig = figure('Color', 'w', 'Position', [40 60 1500 820]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Same-colour pairs inside assemblies, all assemblies of a session pooled (rarity-weighted)', ...
    'FontWeight', 'bold', 'FontSize', 12);
animaux = unique(T.animal, 'stable');
titres = {'A. Main test: valid assemblies, colours shuffled in the field', ...
    'B. Spatial control: colours shuffled within blocks'};
modes = {'globale', 'blocs'};
for g = 1:2
    ax = nexttile(tl); hold(ax, 'on');
    t = T(strcmp(T.assemblees, 'valides') & strcmp(T.melange, modes{g}), :);
    for k = 1:numel(animaux)
        m = strcmp(t.animal, animaux{k});
        if ~any(m), continue; end
        y = 100 * t.exces(m);
        sig = t.p_session(m) < 0.05;
        x = k + 0.12 * randn(sum(m), 1);
        scatter(ax, x(~sig), y(~sig), 36, [0.4 0.4 0.4], 'LineWidth', 1);
        scatter(ax, x(sig), y(sig), 36, [0.4 0.4 0.4], 'filled');
        plot(ax, k + [-0.3 0.3], mean(y) * [1 1], 'k-', 'LineWidth', 2);
    end
    yline(ax, 0, '--', 'Color', [0.6 0.6 0.6]);
    set(ax, 'XTick', 1:numel(animaux), 'XTickLabel', animaux, 'XTickLabelRotation', 30, 'FontSize', 9, ...
        'TickDir', 'out', 'TickLabelInterpreter', 'none');
    xlim(ax, [0.4, numel(animaux) + 0.6]);
    ylabel(ax, 'excess over shuffled colours (%)');
    b = B(strcmp(B.assemblees, 'valides') & strcmp(B.melange, modes{g}), :);
    if isempty(b)
        title(ax, titres{g}, 'FontSize', 10);
    else
        title(ax, {titres{g}, sprintf('%+.1f%%, %d/%d animals > 0, p = %s', 100 * b.exces_moyen, ...
            b.n_animaux_positifs, b.n_animaux, texte_p(b.p_wilcoxon))}, 'FontSize', 10);
    end
    box(ax, 'off');
end

% C. observé contre attendu
ax = nexttile(tl); hold(ax, 'on');
t = T(strcmp(T.assemblees, 'valides') & strcmp(T.melange, 'globale'), :);
sig = t.p_session < 0.05;
scatter(ax, t.Q_attendu(~sig), t.Q_obs(~sig), 30, [0.4 0.4 0.4], 'LineWidth', 1);
scatter(ax, t.Q_attendu(sig), t.Q_obs(sig), 30, [0.8 0.1 0.1], 'filled');
lim = [0, max([t.Q_attendu; t.Q_obs; 1]) * 1.05];
plot(ax, lim, lim, '--', 'Color', [0.6 0.6 0.6]);
xlim(ax, lim); ylim(ax, lim);
xlabel(ax, 'weighted same-colour pairs expected by chance');
ylabel(ax, 'observed');
title(ax, sprintf('C. Each session (red: p < 0.05, %d/%d)', sum(sig), height(t)), 'FontSize', 10);
box(ax, 'off'); set(ax, 'FontSize', 9, 'TickDir', 'out');

% D. les quatre analyses
ax = nexttile(tl); hold(ax, 'on');
etiq = cell(1, height(B));
for r = 1:height(B)
    a = A(strcmp(A.assemblees, B.assemblees{r}) & strcmp(A.melange, B.melange{r}), :);
    mo = 100 * mean(a.exces);
    se = 100 * std(a.exces) / sqrt(max(height(a), 1));
    coul = [0.5 0.5 0.5];
    if r == 1, coul = [0 0 0]; end
    plot(ax, mo + [-se se], [r r], '-', 'Color', coul, 'LineWidth', 2);
    plot(ax, mo, r, 'o', 'Color', coul, 'MarkerFaceColor', coul, 'MarkerSize', 8);
    text(ax, mo + se, r, sprintf('   p = %s', texte_p(B.p_wilcoxon(r))), 'FontSize', 9);
    if strcmp(B.assemblees{r}, 'valides'), n1 = 'valid'; else, n1 = 'all'; end
    if strcmp(B.melange{r}, 'globale'), n2 = 'field shuffle'; else, n2 = 'block shuffle'; end
    etiq{r} = sprintf('%s assemblies, %s', n1, n2);
end
xline(ax, 0, ':', 'Color', [0.5 0.5 0.5]);
set(ax, 'YTick', 1:height(B), 'YTickLabel', etiq, 'YDir', 'reverse', 'FontSize', 9, 'TickDir', 'out');
ylim(ax, [0.5, height(B) + 0.5]);
xlabel(ax, 'excess (%, mean of animals \pm sem)');
title(ax, 'D. Sensitivity', 'FontSize', 10);
box(ax, 'off');
exportgraphics(fig, fullfile(dossier, 'paires_assemblees.png'), 'Resolution', 200);
end


function t = texte_p(p)
if isnan(p), t = 'NA'; elseif p < 0.001, t = '< 0.001'; else, t = sprintf('%.3f', p); end
end
