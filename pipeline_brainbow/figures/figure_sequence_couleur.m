function fig = figure_sequence_couleur(dossier)
%FIGURE_SEQUENCE_COULEUR  Étape 5 : les cellules de même couleur sont-elles
%   plus proches dans la séquence ? Sans rien recalculer : lit
%   <PathSave>\sequences_agregation\sequences_agregation.mat.
%     A. distance dans la séquence (fraction de sa longueur), paires de même
%        couleur contre couleurs différentes : un trait par animal (moyenne
%        de ses sessions), noir = moyenne des animaux +/- sem ;
%        p = Wilcoxon unilatéral sur les animaux.
%     B. chaque session, rangée par animal : écart en % (même / différentes
%        - 1) ; rond plein : p < 0,05 (couleurs mélangées) ; trait noir :
%        moyenne de l'animal.
%     C. contrôle : rho de Spearman entre distance dans la séquence et
%        distance dans l'espace (paires de chaque session).
%   Écrit sequence_colour.png.

dossier = regexprep(char(dossier), '[\\/]+$', '');
S = load(fullfile(dossier, 'sequences_agregation.mat'));
fig = figure('Color', 'w', 'Position', [60 60 1500 480]);
if ~isfield(S, 'Ac') || isempty(S.Ac)
    annotation(fig, 'textbox', [0.1 0.4 0.8 0.2], 'String', ...
        'No session with enough same- and different-colour pairs in the sequence', ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);
    exportgraphics(fig, fullfile(dossier, 'sequence_colour.png'), 'Resolution', 200);
    return
end
Ac = S.Ac; b = S.bc;
T = S.T(isfinite(S.T.delta_couleur), :);
GRIS = [0.7 0.7 0.7];
GRIS_F = [0.45 0.45 0.45];
tl = tiledlayout(fig, 1, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Same colour closer in the running sequence? (%d animals, %d sessions)', b.n_animaux, ...
    b.n_sessions), 'FontWeight', 'bold', 'FontSize', 12);

% A. même couleur contre couleurs différentes
ax = nexttile(tl); hold(ax, 'on');
y = [Ac.dist_meme, Ac.dist_diff];
for a = 1:size(y, 1)
    plot(ax, [1 2], y(a, :), '-', 'Color', GRIS);
end
errorbar(ax, [1 2], mean(y, 1), std(y, 0, 1) / sqrt(size(y, 1)), '-o', 'Color', 'k', 'LineWidth', 2, ...
    'MarkerFaceColor', 'k', 'MarkerSize', 6, 'CapSize', 6);
set(ax, 'XTick', [1 2], 'XTickLabel', {'same colour', 'different colour'}, 'XLim', [0.6 2.4], ...
    'FontSize', 9, 'TickDir', 'out');
ylabel(ax, 'distance in the sequence (fraction)');
title(ax, {'A. Distance in the sequence', sprintf('%+.1f%%, %s', b.pct_moyen, texte_p(b.p_wilcoxon))}, ...
    'FontSize', 10);
box(ax, 'off');

% B. chaque session
animaux = Ac.animal;
nA = numel(animaux);
[~, ia] = ismember(T.animal, animaux);
rng(1);
ax = nexttile(tl, [1 2]); hold(ax, 'on');
for a = 1:nA
    m = find(ia == a);
    x = a + 0.16 * (rand(numel(m), 1) - 0.5) * 2;
    sig = T.p_couleur(m) < 0.05;
    scatter(ax, x(~sig), T.pct_couleur(m(~sig)), 40, GRIS_F, 'LineWidth', 1);
    scatter(ax, x(sig), T.pct_couleur(m(sig)), 40, GRIS_F, 'filled');
    plot(ax, a + [-0.25 0.25], Ac.pct(a) * [1 1], 'k-', 'LineWidth', 2);
end
yline(ax, 0, '--', 'Color', [0.6 0.6 0.6]);
set(ax, 'XTick', 1:nA, 'XTickLabel', animaux, 'XLim', [0.4 nA + 0.6], 'TickLabelInterpreter', 'none', ...
    'FontSize', 9, 'TickDir', 'out');
xtickangle(ax, 30);
ylabel(ax, 'same vs different colour (%, < 0: closer)');
title(ax, {'B. Each session, by animal', sprintf('%d/%d sessions p < 0.05', b.n_sig, b.n_sessions)}, ...
    'FontSize', 10);
box(ax, 'off');

% C. contrôle : la séquence suit-elle l'espace ?
ax = nexttile(tl, [1 2]); hold(ax, 'on');
for a = 1:nA
    m = find(ia == a);
    x = a + 0.16 * (rand(numel(m), 1) - 0.5) * 2;
    scatter(ax, x, T.rho_espace(m), 36, GRIS_F, 'filled');
    plot(ax, a + [-0.25 0.25], Ac.rho_espace(a) * [1 1], 'k-', 'LineWidth', 2);
end
yline(ax, 0, '--', 'Color', [0.6 0.6 0.6]);
set(ax, 'XTick', 1:nA, 'XTickLabel', animaux, 'XLim', [0.4 nA + 0.6], 'TickLabelInterpreter', 'none', ...
    'FontSize', 9, 'TickDir', 'out');
xtickangle(ax, 30);
ylabel(ax, 'rho (sequence distance vs space)');
title(ax, {'C. Does the sequence follow space?', sprintf('mean rho = %.3f, %s', b.rho_espace_moyen, ...
    texte_p(b.p_espace))}, 'FontSize', 10);
box(ax, 'off');
exportgraphics(fig, fullfile(dossier, 'sequence_colour.png'), 'Resolution', 200);
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
