function fig = figure_couleurs(colorcell, Donnees_Clustering, titre)
%% =========================================================================
%  FIGURE_COULEURS (17 nov., demandée par l'utilisateur)
%  Une seule figure, deux graphes :
%    à gauche : triangle de Maxwell, chaque cellule colorée par sa classe
%               (survol : proportions R V B, numéro, classe) ;
%    à droite : proportion de chaque couleur (cellules des classes 1 à 6),
%               à côté des proportions de référence : médiane des champs de
%               référence, barre = du plus petit au plus grand champ
%               (relevés à l'œil sur le violin plot « Distribution des
%               pourcentages par couleur », ~17 champs).
%               * = en dehors de l'étendue des champs de référence.
%
%  Lancer :  figure_couleurs(colorcell)                       (barres seules)
%            figure_couleurs(colorcell, Donnees_Clustering)   (+ triangle)
%            figure_couleurs(colorcell, Donnees_Clustering, 'titre')
%  Classes : 1 rouge, 2 vert, 3 bleu, 4 jaune, 5 magenta, 6 cyan, 7 blanc, 8 aucune.
%
%  Les valeurs sont lues sur l'image (précision ~0.5 %) : les remplacer par
%  les vraies si on a les données du violin plot.
%  (props_theo_raw de check_colorcell.m, jaune 29.4 / magenta 20.8, a le
%  jaune et le magenta inversés par rapport au violin plot.)
%  =========================================================================

if nargin < 2, Donnees_Clustering = []; end
if nargin < 3, titre = ''; end

%              Red   Green Blue  Yellow Magenta Cyan     (en %)
REF_MED   = [34.0,  5.5,  5.0,  19.5,  27.5,  10.0];
REF_MIN   = [29.3,  3.2,  1.1,  12.7,  19.6,   5.3];
REF_MAX   = [38.0, 13.4,  8.3,  24.2,  37.2,  16.2];

NOMS = {'Red', 'Green', 'Blue', 'Yellow', 'Magenta', 'Cyan', 'White', 'None'};
PALETTE = [
    1.0 0.0 0.0;    % 1 red
    0.0 0.8 0.0;    % 2 green
    0.0 0.0 1.0;    % 3 blue
    0.9 0.9 0.0;    % 4 yellow
    1.0 0.0 1.0;    % 5 magenta
    0.0 0.8 0.8;    % 6 cyan
    0.92 0.92 0.92; % 7 white
    0.1 0.1 0.1];   % 8 none

c = double(colorcell(:));
c(~ismember(c, 1:8)) = 8;           % 0 ou NaN (anciens fichiers) -> aucune

fig = figure('Name', ['Colours ' titre], 'Color', 'w', 'Position', [80 120 1250 540]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
if ~isempty(titre), title(tl, titre, 'Interpreter', 'none', 'FontWeight', 'bold'); end

%% ---------- triangle de Maxwell ----------
ax1 = nexttile(tl);
hold(ax1, 'on');
if isempty(Donnees_Clustering) || size(Donnees_Clustering, 1) ~= numel(c)
    axis(ax1, 'off');
    if isempty(Donnees_Clustering)
        msg = 'No triangle data (run step 7, or no couleurs\_details file)';
    else
        msg = sprintf('Triangle data: %d cells, colorcell: %d', size(Donnees_Clustering, 1), numel(c));
    end
    text(ax1, 0.5, 0.5, msg, 'HorizontalAlignment', 'center', 'FontSize', 11);
else
    X = Donnees_Clustering(:, 1);
    Y = Donnees_Clustering(:, 2);
    h = sqrt(3) / 2;
    plot(ax1, [0 1 0.5 0], [0 0 h 0], 'k-', 'LineWidth', 1.5);
    % lignes de grille tous les 0.25
    for f = 0.25:0.25:0.75
        plot(ax1, [f, 0.5 + f / 2], [0, (1 - f) * h], ':', 'Color', [0.8 0.8 0.8]);    % r constant
        plot(ax1, [f / 2, 1 - f / 2], [f * h, f * h], ':', 'Color', [0.8 0.8 0.8]);   % b constant
        plot(ax1, [f, f / 2], [0, f * h], ':', 'Color', [0.8 0.8 0.8]);               % g constant
    end
    text(ax1, 1.02, 0, 'Red', 'FontWeight', 'bold', 'FontSize', 12, 'Color', [0.8 0 0]);
    text(ax1, -0.02, 0, 'Green', 'FontWeight', 'bold', 'FontSize', 12, ...
        'HorizontalAlignment', 'right', 'Color', [0 0.6 0]);
    text(ax1, 0.5, h + 0.04, 'Blue', 'FontWeight', 'bold', 'FontSize', 12, ...
        'HorizontalAlignment', 'center', 'Color', [0 0 0.8]);
    hs = scatter(ax1, X, Y, 30, PALETTE(c, :), 'filled', 'MarkerEdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.5);
    % infobulles : proportions, numéro, classe
    B = Y / h; R = X - 0.5 * B; V = 1 - R - B;
    hs.DataTipTemplate.DataTipRows = [ ...
        dataTipTextRow('Cell', (1:numel(c))', '%d'), ...
        dataTipTextRow('Class', NOMS(c)'), ...
        dataTipTextRow('R', R, '%.3f'), ...
        dataTipTextRow('G', V, '%.3f'), ...
        dataTipTextRow('B', B, '%.3f')];
    axis(ax1, 'equal');
    axis(ax1, 'off');
    xlim(ax1, [-0.15 1.15]); ylim(ax1, [-0.05 h + 0.1]);
    title(ax1, sprintf('Maxwell triangle (%d cells)', numel(c)));
end
hold(ax1, 'off');

%% ---------- proportions vs référence ----------
ax2 = nexttile(tl);
n6 = arrayfun(@(k) sum(c == k), 1:6);
N = sum(n6);
if N == 0
    axis(ax2, 'off');
    text(ax2, 0.5, 0.5, 'No cell in classes 1 to 6', 'HorizontalAlignment', 'center');
    return
end
pct = 100 * n6 / N;
hors = pct < REF_MIN | pct > REF_MAX;      % en dehors de l'étendue des champs de référence

hold(ax2, 'on');
b = bar(ax2, 1:6, [pct(:), REF_MED(:)], 'grouped');
b(1).FaceColor = 'flat'; b(1).CData = PALETTE(1:6, :); b(1).EdgeColor = 'k';
b(2).FaceColor = 'flat'; b(2).CData = PALETTE(1:6, :); b(2).EdgeColor = 'k'; b(2).FaceAlpha = 0.3;
errorbar(ax2, b(2).XEndPoints, REF_MED, REF_MED - REF_MIN, REF_MAX - REF_MED, 'k', ...
    'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 6);
for k = 1:6
    t = sprintf('%d', n6(k));
    if hors(k), t = [t ' *']; end
    text(ax2, b(1).XEndPoints(k), pct(k) + 1, t, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 9, 'FontWeight', ifelse(hors(k), 'bold', 'normal'));
end
ax2.XTick = 1:6; ax2.XTickLabel = NOMS(1:6);
ylabel(ax2, '% of coloured cells (classes 1-6)');
ylim(ax2, [0, max([pct, REF_MAX]) * 1.15 + 2]);
box(ax2, 'off'); grid(ax2, 'on'); ax2.XGrid = 'off';

% test du chi-2 contre la référence (médianes normalisées) ; sans boîte à outils
attendu = N * REF_MED / sum(REF_MED);
chi2 = sum((n6 - attendu) .^ 2 ./ attendu);
p = gammainc(chi2 / 2, 5 / 2, 'upper');
title(ax2, sprintf('N = %d coloured (white %d, none %d) ; \\chi^2 vs ref p = %.2g', ...
    N, sum(c == 7), sum(c == 8), p));
l1 = bar(ax2, NaN, NaN, 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', 'k');
l2 = bar(ax2, NaN, NaN, 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', 'k', 'FaceAlpha', 0.3);
legend(ax2, [l1 l2], {'This field', 'Reference (median, min-max)'}, 'Location', 'northeast');
xlim(ax2, [0.5 6.5]);           % les barres factices de la légende décalent sinon l'axe
hold(ax2, 'off');

% résumé dans la fenêtre de commande
fprintf('\n%-8s %6s %8s %9s %14s\n', 'Colour', 'n', 'obs %', 'ref med', 'ref min-max');
for k = 1:6
    fprintf('%-8s %6d %7.1f %9.1f %8.1f-%-5.1f%s\n', NOMS{k}, n6(k), pct(k), REF_MED(k), REF_MIN(k), ...
        REF_MAX(k), ifelse(hors(k), '  *', ''));
end
fprintf('white %d, none %d ; chi2 = %.1f, p = %.2g\n', sum(c == 7), sum(c == 8), chi2, p);
end


function v = ifelse(test, a, b)
if test, v = a; else, v = b; end
end
