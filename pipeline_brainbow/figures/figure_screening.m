function figure_screening(S, dossier)
%FIGURE_SCREENING  Figure du screening des réglages (13 nov.).
%   figure_screening(S, dossier)   S : sortie de screening_agregation
%   figure_screening(dossier)      relit <dossier>\screening.mat
%   A. excès de toutes les combinaisons éligibles, et meilleur excès de
%      chaque mélange des couleurs (ce que donne le tri seul)
%   B à F. excès moyen de chaque niveau de chaque réglage (point noir), ce
%      que donne le hasard (barre grise, 95 %) et la meilleure combinaison
%      qui a ce niveau (triangle rouge)
%   Écrit <dossier>\figure_screening.png

if nargin < 2 && (ischar(S) || isstring(S))
    dossier = char(S);
    Sf = load(fullfile(dossier, 'screening.mat'), 'S');
    S = Sf.S;
end
C = S.C; F = S.F; G = S.G;
el = C.eligible > 0;
if ~any(el), return; end
B = numel(S.max_nul);

fig = figure('Color', 'w', 'Position', [40 60 1600 900]);
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Settings screening: %s (%d sessions, %d eligible settings)', strjoin(S.animaux, ', '), ...
    numel(S.sessions), sum(el)), 'FontSize', 11, 'Interpreter', 'none');

% A. toutes les combinaisons, et le meilleur de chaque mélange
ax = nexttile(tl); hold(ax, 'on');
x = 100 * C.exces(el);
xn = 100 * S.max_nul(:);
xn = xn(isfinite(xn));
bords = linspace(min([x; xn]) - 1, max([x; xn]) + 1, 45);
histogram(ax, x, bords, 'Normalization', 'probability', 'FaceColor', [0.25 0.45 0.75], 'EdgeColor', 'none');
histogram(ax, xn, bords, 'Normalization', 'probability', 'FaceColor', [0.55 0.55 0.55], 'EdgeColor', 'none', ...
    'FaceAlpha', 0.6);
xline(ax, max(x), '-', 'Color', [0.8 0.1 0.1], 'LineWidth', 1.5);
xline(ax, 0, '--', 'Color', [0.5 0.5 0.5]);
legend(ax, {'all tested settings (real colours)', 'best setting of each colour shuffle', ...
    'best setting (real colours)'}, 'Location', 'northoutside', 'Box', 'off', 'FontSize', 8);
xlabel(ax, 'excess of biased assemblies (points above chance)');
ylabel(ax, 'fraction');
p_sel = (1 + sum(S.max_nul >= max(C.exces(el)) - 1e-12)) / (B + 1);
title(ax, sprintf('A. Best setting vs selection alone: p = %.3g', p_sel), 'FontSize', 10);
set(ax, 'FontSize', 9, 'TickDir', 'out'); box(ax, 'off');

% B. prétraitement
f1 = F(1);
noms = cell(1, numel(G.pre));
for p = 1:numel(G.pre)
    if strcmp(G.pre(p).facteur, 'base')
        noms{p} = 'base';
    else
        noms{p} = sprintf('%s: %s', G.pre(p).facteur, G.pre(p).niveau);
    end
end
f1.niveaux = noms;
panneau(nexttile(tl), f1, 'B. Preprocessing (one change at a time)', true);
% C à F. réglages des assemblées
panneau(nexttile(tl), F(2), 'C. SCE clusters kept (rule, percentile)', false);
panneau(nexttile(tl), F([3 4]), 'D. Members rule and assembly version', false);
panneau(nexttile(tl), F(5), 'E. Choice of K', false);
panneau(nexttile(tl), F([6 7]), 'F. Step-4 criteria (min cells, min SCE)', false);

exportgraphics(fig, fullfile(dossier, 'figure_screening.png'), 'Resolution', 150);
end


function panneau(ax, Fs, titre, avec_legende)
% une ligne par niveau : hasard (barre grise), moyenne (point), meilleure combinaison (triangle)
hold(ax, 'on');
y = 0; yt = []; yl = {};
h = gobjects(1, 3);
for f = 1:numel(Fs)
    for q = 1:numel(Fs(f).niveaux)
        if Fs(f).n(q) == 0, continue; end
        y = y + 1;
        h(1) = plot(ax, 100 * [Fs(f).nul_bas(q), Fs(f).nul_haut(q)], [y y], '-', 'Color', [0.78 0.78 0.78], ...
            'LineWidth', 6);
        h(2) = plot(ax, 100 * Fs(f).moyenne(q), y, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', ...
            'MarkerSize', 5);
        h(3) = plot(ax, 100 * Fs(f).meilleur(q), y, '^', 'MarkerFaceColor', [0.8 0.1 0.1], ...
            'MarkerEdgeColor', [0.8 0.1 0.1], 'MarkerSize', 4);
        yt(end+1) = y; %#ok<AGROW>
        if numel(Fs) > 1
            yl{end+1} = sprintf('%s: %s', Fs(f).nom, Fs(f).niveaux{q}); %#ok<AGROW>
        else
            yl{end+1} = Fs(f).niveaux{q}; %#ok<AGROW>
        end
    end
    y = y + 0.6;
end
xline(ax, 0, '--', 'Color', [0.5 0.5 0.5]);
if isempty(yt)
    axis(ax, 'off');
    return
end
set(ax, 'YTick', yt, 'YTickLabel', yl, 'YDir', 'reverse', 'TickLabelInterpreter', 'none', 'FontSize', 8, ...
    'TickDir', 'out', 'YLim', [0.3, y]);
xlabel(ax, 'excess (points above chance)');
title(ax, titre, 'FontSize', 10);
if avec_legende && all(isgraphics(h))
    legend(ax, h, {'chance (95%)', 'mean of all settings with this level', 'best setting with this level'}, ...
        'Location', 'southoutside', 'Box', 'off', 'FontSize', 7, 'NumColumns', 3);
end
box(ax, 'off');
end
