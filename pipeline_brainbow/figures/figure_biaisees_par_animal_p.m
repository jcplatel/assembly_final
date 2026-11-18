function fig = figure_biaisees_par_animal_p(dossier_variante, fdr, alpha, nom_variante)
%FIGURE_BIAISEES_PAR_ANIMAL_P  Comme figure_biaisees_par_animal, mais avec la
%   p hypergéométrique de chaque couleur (p_brut, non corrigée pour les 6
%   couleurs) au lieu de p_minp. Un symbole par paire assemblée x couleur
%   enrichie (une assemblée peut donc avoir 2 symboles). (17 nov., demandée
%   par l'utilisateur)
%
%   figure_biaisees_par_animal_p('E:\...\agregation\champ_globale')
%
%   Assemblées testables seulement. Écrit biased_assemblies_per_animal_p.png.
%     A. une colonne par animal :
%          rond plein   q < fdr : FDR (Benjamini-Hochberg) sur toutes les paires
%                       assemblée x couleur (comme paires_assemblee_couleur_FDR.csv)
%          rond vide    p_minp < alpha (couleur en excès de l'assemblée,
%                       corrigée pour les 6 couleurs), pas après FDR
%          carré vide   p < alpha sans correction seulement
%          trait noir   nombre de paires attendu à p < alpha par hasard (calcul
%                       exact : loi hypergéométrique, composition du champ)
%     B. % de paires à p < alpha, observé contre attendu, par animal
%   Aucun test nouveau, sauf le nombre attendu.

if nargin < 2 || isempty(fdr), fdr = 0.10; end
if nargin < 3 || isempty(alpha), alpha = 0.05; end
if nargin < 4, nom_variante = ''; end
dossier_variante = regexprep(char(dossier_variante), '[\\/]+$', '');
S = load(fullfile(dossier_variante, 'resultats_biais_clonal.mat'), 'R');
R = S.R;
if isempty(nom_variante), [~, nom_variante] = fileparts(dossier_variante); end
noms_coul = {'R', 'G', 'B', 'RG', 'RB', 'GB'};
cols = [0.86 0.16 0.16; 0.16 0.66 0.22; 0.20 0.38 0.92; ...
        0.93 0.78 0.10; 0.82 0.22 0.80; 0.10 0.74 0.80];

Rt = R(R.testable, :);
P = Rt.p_brut;                                   % n x 6, NaN = couleur non testée
ok = ~isnan(P);
Q = nan(size(P)); Q(ok) = bh_q(P(ok));            % FDR sur toutes les paires
[~, ic_dom] = ismember(Rt.couleur_dominante, noms_coul);

% nombre attendu de paires à p < alpha par hasard (exact)
Kref = Rt.Kcol_urne;
if ismember('Kcol_ref', Rt.Properties.VariableNames) && ~any(isnan(Rt.Kcol_ref(:))), Kref = Rt.Kcol_ref; end
att_paire = zeros(size(P));
for i = 1:height(Rt)
    N = sum(Kref(i, :)); n = Rt.n_colorees(i);
    for c = find(ok(i, :))
        k = 0:n;
        pk = hygecdf(k - 1, N, Kref(i, c), n, 'upper');     % p si on observe k cellules
        pk(1) = 1;
        att_paire(i, c) = sum(hygepdf(k(pk < alpha), N, Kref(i, c), n));
    end
end

animaux = unique(R.animal);
nA = numel(animaux);
fig = figure('Color', 'w', 'Position', [60 60 1400 600]);
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
n_q = 0; n_m = 0; n_s = 0;
n_t = zeros(nA, 1); n_paires = zeros(nA, 1); n_sig = zeros(nA, 1); n_att = zeros(nA, 1);

ax = nexttile(tl, [1 2]); hold(ax, 'on');
ymax = 1;
for a = 1:nA
    m = find(strcmp(Rt.animal, animaux{a}));
    n_t(a) = numel(m);
    n_paires(a) = sum(sum(ok(m, :)));
    n_att(a) = sum(sum(att_paire(m, :)));
    lst = zeros(0, 2);                           % [niveau couleur] ; 1 FDR, 2 p_minp, 3 p seul
    for i = m(:)'
        for c = find(ok(i, :) & P(i, :) < alpha)
            if Q(i, c) < fdr, niv = 1;
            elseif c == ic_dom(i) && Rt.p_minp(i) < alpha, niv = 2;
            else, niv = 3;
            end
            lst(end+1, :) = [niv c]; %#ok<AGROW>
        end
    end
    lst = sortrows(lst, [1 2]);
    n_sig(a) = size(lst, 1);
    n_q = n_q + sum(lst(:, 1) == 1); n_m = n_m + sum(lst(:, 1) == 2); n_s = n_s + sum(lst(:, 1) == 3);
    for y = 1:size(lst, 1)
        c = cols(lst(y, 2), :);
        switch lst(y, 1)
            case 1, plot(ax, a, y, 'o', 'MarkerSize', 13, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
            case 2, plot(ax, a, y, 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', c, 'LineWidth', 2.2);
            case 3, plot(ax, a, y, 's', 'MarkerSize', 11, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', c, 'LineWidth', 2.2);
        end
    end
    plot(ax, a + [-0.32 0.32], n_att(a) * [1 1], 'k-', 'LineWidth', 2);
    text(ax, a, -0.35, sprintf('n = %d', n_t(a)), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'FontSize', 8, 'Color', [0.35 0.35 0.35]);
    ymax = max([ymax, size(lst, 1), n_att(a)]);
end
set(ax, 'XTick', 1:nA, 'XTickLabel', animaux, 'XLim', [0.4 nA + 0.6], 'YLim', [-1.1 ymax + 1], ...
    'YTick', 0:max(1, ceil(ymax / 8)):ymax + 1, 'TickLabelInterpreter', 'none', 'FontSize', 10, 'TickDir', 'out');
ax.YGrid = 'on';
ylabel(ax, 'number of assembly x colour pairs');
title(ax, 'A. One symbol per enriched assembly x colour pair (p of that colour, not corrected)', 'FontSize', 10);
h = gobjects(0); t = {};
h(end+1) = plot(ax, nan, nan, 'o', 'MarkerSize', 11, 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'k');
t{end+1} = sprintf('FDR on all pairs, q < %g', fdr);
h(end+1) = plot(ax, nan, nan, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', [0.5 0.5 0.5], 'LineWidth', 2);
t{end+1} = sprintf('p_{minP} < %g (corrected for 6 colours)', alpha);
h(end+1) = plot(ax, nan, nan, 's', 'MarkerSize', 10, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', [0.5 0.5 0.5], 'LineWidth', 2);
t{end+1} = sprintf('p < %g, no correction only', alpha);
h(end+1) = plot(ax, nan, nan, 'k-', 'LineWidth', 2);
t{end+1} = sprintf('pairs expected at p < %g by chance', alpha);
for k = 1:6
    h(end+1) = plot(ax, nan, nan, 'o', 'MarkerSize', 9, 'MarkerFaceColor', cols(k, :), 'MarkerEdgeColor', 'none'); %#ok<AGROW>
    t{end+1} = noms_coul{k}; %#ok<AGROW>
end
legend(ax, h, t, 'Location', 'northoutside', 'Orientation', 'horizontal', 'NumColumns', 5, 'Box', 'off', 'FontSize', 8);
box(ax, 'off');

ax = nexttile(tl); hold(ax, 'on');
okA = n_paires > 0;
obs = 100 * n_sig ./ max(n_paires, 1); att = 100 * n_att ./ max(n_paires, 1);
for a = find(okA)'
    plot(ax, [1 2], [att(a) obs(a)], '-', 'Color', [0.7 0.7 0.7]);
end
scatter(ax, ones(sum(okA), 1), att(okA), 45, [0.45 0.45 0.45], 'filled');
scatter(ax, 2 * ones(sum(okA), 1), obs(okA), 45, [0.75 0.15 0.15], 'filled');
plot(ax, [0.8 1.2], mean(att(okA)) * [1 1], 'k-', 'LineWidth', 2);
plot(ax, [1.8 2.2], mean(obs(okA)) * [1 1], 'k-', 'LineWidth', 2);
p_w = signrank(obs(okA), att(okA), 'tail', 'right');
set(ax, 'XTick', [1 2], 'XTickLabel', {'expected', 'observed'}, 'XLim', [0.5 2.5], 'FontSize', 10, 'TickDir', 'out');
ylabel(ax, sprintf('%% of tested pairs at p < %g (no correction)', alpha));
title(ax, {'B. Each animal (one line)', sprintf('Wilcoxon one-sided on animals: p = %.3f', p_w)}, 'FontSize', 10);
box(ax, 'off');

title(tl, {sprintf(['Enriched pairs per animal (p of each colour instead of p_{minP}): %d after FDR + %d at p_{minP} < %g + ' ...
    '%d at p < %g only; %d pairs tested in %d testable assemblies, %.1f expected at p < %g'], ...
    n_q, n_m, alpha, n_s, alpha, sum(n_paires), sum(n_t), sum(n_att), alpha), strrep(nom_variante, '_', '\_')}, ...
    'FontWeight', 'bold', 'FontSize', 11, 'Interpreter', 'tex');
exportgraphics(fig, fullfile(dossier_variante, 'biased_assemblies_per_animal_p.png'), 'Resolution', 200);
end


function q = bh_q(p)
% q-values de Benjamini-Hochberg
p = p(:); n = numel(p);
[ps, o] = sort(p);
qs = min(1, cummin(ps .* n ./ (1:n)', 'reverse'));
q = zeros(n, 1); q(o) = qs;
end
