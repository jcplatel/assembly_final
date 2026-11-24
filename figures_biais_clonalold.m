%% =========================================================================
%  BIAIS CLONAL — FIGURE DES ANALYSES COMPLÉMENTAIRES (12 panneaux, en anglais)
%  -------------------------------------------------------------------------
%  Lit resultats_biais_clonal.mat (table R, écrite par agregation_stats_clonales.m)
%  et produit clonal_bias_complementary_analyses.png dans dossier_racine :
%
%    A. % d'assemblées à p < 0,05 par animal (IC de Wilson) + estimation
%       globale (rééchantillonnage des animaux) et GLMM
%    B. distribution cumulée des p minP (sans effet : sur la diagonale)
%    C. retrait de chaque animal
%    D. biais selon la compacité, par classe de taille
%    E. p minP selon la compacité (rapport d'étalement si disponible)
%    F. p minP selon le nombre de cellules colorées
%    G. part des SCE de la session où l'assemblée est active (%)
%    H. nombre de SCE où l'assemblée est active (n)
%    I. recrutement : fraction des cellules de l'assemblée active par SCE
%    J. recrutement selon la taille (tendance + écart des biaisées)
%    K. taille de l'assemblée (% des cellules du champ)
%    L. taille de l'assemblée (nombre de cellules)
%
%  Catégories (panneaux E-L) : biaisée (q < fdr_fig), non biaisée
%  (p > seuil_non_biaise_fig), indéterminée (entre les deux).
%  Tout est descriptif / exploratoire, sauf A (estimation et IC).
%  Peut être lancé seul, ou à la fin de agregation_stats_clonales.m.
%  =========================================================================

%% ---------- Paramètres --------------------------------------------------
if ~exist('dossier_racine', 'var')
    dossier_racine = 'E:\Data\Aurelie\analysis\Final\assembly\nocues\';
end
seuil_p              = 0.05;   % seuil par assemblée (p minP)
fdr_fig              = 0.20;   % FDR qui définit les assemblées « biaisées »
seuil_non_biaise_fig = 0.5;    % p minP > 0,5 -> « non biaisée »
seuil_compact        = 5;      % percentile de compacité <= 5 -> « très compacte »
graine_fig           = 1;

%% ---------- Données -----------------------------------------------------
S = load(fullfile(dossier_racine, 'resultats_biais_clonal.mat'), 'R');
T = S.R(S.R.testable, :);
T.sig     = double(T.p_minp < seuil_p);
T.mlogp   = -log10(T.p_minp);
T.biaisee = T.q < fdr_fig;
a_compac  = ismember('compacite_pct', T.Properties.VariableNames) && any(~isnan(T.compacite_pct));
if a_compac, T.compacte = T.compacite_pct <= seuil_compact; end
a_ratio   = ismember('ratio_etalement', T.Properties.VariableNames) && any(~isnan(T.ratio_etalement));

ROUGE = [200 40 62] / 255;   GRIS   = [138 138 133] / 255;
BLEU  = [46 109 180] / 255;  AMBRE  = [227 161 59] / 255;
ENCRE = [0.12 0.12 0.12];    ENCRE2 = [0.37 0.37 0.35];
CLAIR = [0.88 0.88 0.85];    GRIS_C = [201 199 190] / 255;
PTICKS = {'1', '0.1', '0.01', '0.001', '0.0001'};

cat_T = repmat({'indeterminate'}, height(T), 1);
cat_T(T.p_minp > seuil_non_biaise_fig) = {'not_biased'};
cat_T(T.biaisee) = {'biased'};
% nom interne, légende, nom court (axe x), couleur, opacité
CATS = {'not_biased',    sprintf('not biased (p > %.1f)', seuil_non_biaise_fig), 'not biased',    GRIS,  0.5; ...
        'indeterminate', 'indeterminate',                                     'indeterminate', AMBRE, 0.85; ...
        'biased',        sprintf('biased (FDR %.0f%%)', 100 * fdr_fig),        'biased',        ROUGE, 1};

%% ---------- Chiffres ----------------------------------------------------
anim_u = unique(T.animal);
k_a = cellfun(@(a) sum(T.sig(strcmp(T.animal, a))), anim_u);
n_a = cellfun(@(a) sum(strcmp(T.animal, a)), anim_u);
[~, ordre] = sort(k_a ./ n_a, 'descend');
anim_u = anim_u(ordre); k_a = k_a(ordre); n_a = n_a(ordre);
p_a = k_a ./ n_a;
[lo_a, hi_a] = wilson(k_a, n_a);

rng(graine_fig);
sig_par = cellfun(@(a) T.sig(strcmp(T.animal, a)), anim_u, 'UniformOutput', false);
nB = 10000; boot = zeros(nB, 1);
for b = 1:nB
    boot(b) = mean(vertcat(sig_par{randi(numel(anim_u), numel(anim_u), 1)}));
end
p_glob  = mean(T.sig);
ic_boot = prctile(boot, [2.5 97.5]);

glmm = [NaN NaN NaN];
try
    lg = log(seuil_p / (1 - seuil_p));
    g0 = fitglme(T, 'sig ~ 1 + (1|animal_c) + (1|champ_c)', 'Distribution', 'Binomial', ...
        'Offset', repmat(lg, height(T), 1));
    co = g0.Coefficients;
    glmm = 1 ./ (1 + exp(-(lg + [co.Estimate(1), co.Lower(1), co.Upper(1)])));
catch err_g
    fprintf('GLMM impossible (%s) : panneau A sans la ligne GLMM.\n', err_g.message);
end
p_sans = arrayfun(@(i) mean(T.sig(~strcmp(T.animal, anim_u{i}))), (1:numel(anim_u))');

%% ---------- Figure ------------------------------------------------------
fig = figure('Color', 'w', 'Position', [40 20 1500 1800]);
tl = tiledlayout(fig, 4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Clonal bias of assemblies: complementary analyses (%d testable assemblies, %d animals, %d fields)', ...
    height(T), numel(anim_u), numel(unique(T.champ_u))), 'FontWeight', 'bold', 'FontSize', 13);

% ===== A. Biased assemblies per animal ===================================
ax = nexttile(tl); hold(ax, 'on');
lab = [anim_u; {'All (bootstrap)'; 'GLMM'}];
val = 100 * [p_a; p_glob; glmm(1)];
lo  = 100 * [lo_a; ic_boot(1); glmm(2)];
hi  = 100 * [hi_a; ic_boot(2); glmm(3)];
ny = numel(lab); y = (ny:-1:1)';
for i = 1:ny
    resume = i > numel(anim_u);
    if resume, col = ENCRE; mk = 'd'; lw = 2.2; else, col = BLEU; mk = 'o'; lw = 1.4; end
    if isnan(val(i)), continue; end
    plot(ax, [lo(i) hi(i)], [y(i) y(i)], '-', 'Color', col, 'LineWidth', lw);
    plot(ax, val(i), y(i), mk, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', 'MarkerSize', 7);
end
txt = [arrayfun(@(k, n) sprintf('%d/%d', k, n), k_a, n_a, 'UniformOutput', false); ...
       {sprintf('%d/%d', sum(T.sig), height(T))}; {''}];
text(ax, repmat(101, ny, 1), y, txt, 'FontSize', 8, 'Color', ENCRE2);
xline(ax, 100 * seuil_p, '--', 'Color', ROUGE, 'LineWidth', 1.1);
yline(ax, 2.5, '-', 'Color', CLAIR);
text(ax, 100 * seuil_p + 1, 0.35, sprintf('%.0f%% expected by chance', 100 * seuil_p), 'Color', ROUGE, 'FontSize', 8);
set(ax, 'YTick', 1:ny, 'YTickLabel', flipud(lab), 'XLim', [0 100], 'YLim', [0 ny + 0.6]);
xlabel(ax, '% of assemblies with minP p < 0.05 (95% CI)');
titre_gauche(ax, 'A. Biased assemblies per animal');
habiller(ax, 'x');

% ===== B. Cumulative distribution of p-values ============================
ax = nexttile(tl); hold(ax, 'on');
ps = sort(T.p_minp); n_p = numel(ps); xg = linspace(0, 1, 401);
hband = fill(ax, [xg fliplr(xg)], 100 * [binoinv(0.025, n_p, xg) fliplr(binoinv(0.975, n_p, xg))] / n_p, ...
    CLAIR, 'EdgeColor', 'none');
hdiag = plot(ax, [0 1], [0 100], '--k', 'LineWidth', 1);
hobs  = stairs(ax, [0; ps; 1], 100 * [0; (1:n_p)' / n_p; 1], 'Color', BLEU, 'LineWidth', 1.8);
pos_txt = [0.30 14; 0.45 30]; xv_list = [0.05 0.25];
for i = 1:2
    xv = xv_list(i); yv = 100 * mean(T.p_minp < xv);
    plot(ax, xv, yv, 'o', 'MarkerFaceColor', ROUGE, 'MarkerEdgeColor', 'w', 'MarkerSize', 6);
    plot(ax, [xv pos_txt(i, 1)], [yv pos_txt(i, 2)], '-', 'Color', ROUGE, 'LineWidth', 0.8);
    text(ax, pos_txt(i, 1) + 0.01, pos_txt(i, 2), {sprintf('p < %.2f: %.0f%% of assemblies', xv, yv), ...
        sprintf('(%.0f%% expected)', 100 * xv)}, 'Color', ROUGE, 'FontSize', 8);
end
[~, p_ks] = kstest(T.p_minp, 'CDF', makedist('Uniform'));
text(ax, 0.02, 72, {sprintf('KS test vs uniform: p = %.3f', p_ks), '(assumes independent', 'assemblies)'}, ...
    'FontSize', 8, 'Color', ENCRE2, 'VerticalAlignment', 'top');
legend(ax, [hband hdiag hobs], {'95% of null draws', 'expected (no effect)', 'observed'}, ...
    'Location', 'northwest', 'Box', 'off');
set(ax, 'XLim', [0 1], 'YLim', [0 100]);
xlabel(ax, 'minP p-value'); ylabel(ax, '% of assemblies with a lower p');
titre_gauche(ax, 'B. Cumulative distribution of p-values');
habiller(ax, 'y');

% ===== C. Leave-one-animal-out ===========================================
ax = nexttile(tl); hold(ax, 'on');
yy = (numel(anim_u):-1:1)';
for i = 1:numel(anim_u)
    plot(ax, [100 * seuil_p, 100 * p_sans(i)], [yy(i) yy(i)], '-', 'Color', CLAIR, 'LineWidth', 1.2);
    plot(ax, 100 * p_sans(i), yy(i), 'o', 'MarkerFaceColor', BLEU, 'MarkerEdgeColor', 'w', 'MarkerSize', 7);
    text(ax, 100 * p_sans(i) + 0.3, yy(i), sprintf('%.1f%%', 100 * p_sans(i)), 'FontSize', 8, ...
        'BackgroundColor', 'w', 'Margin', 0.5);
end
h1 = xline(ax, 100 * p_glob, '-k', 'LineWidth', 1.1);
h2 = xline(ax, 100 * seuil_p, '--', 'Color', ROUGE, 'LineWidth', 1.1);
set(ax, 'YTick', 1:numel(anim_u), 'YTickLabel', strcat({'without '}, flipud(anim_u)), ...
    'XLim', [100 * seuil_p - 1, 100 * max(p_sans) + 2.5], 'YLim', [0.4 numel(anim_u) + 0.6]);
legend(ax, [h1 h2], {sprintf('all animals: %.1f%%', 100 * p_glob), sprintf('chance: %.0f%%', 100 * seuil_p)}, ...
    'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');
xlabel(ax, '% of assemblies with p < 0.05 (other animals)');
titre_gauche(ax, 'C. Leave-one-animal-out');
habiller(ax, 'x');

% ===== D. Bias vs compactness, at equal size =============================
ax = nexttile(tl); hold(ax, 'on');
if a_compac
    classe = discretize(T.n_colorees, [0 12 25 Inf], 'IncludedEdge', 'right');
    V = zeros(3, 2); Kt = cell(3, 2);
    for c = 1:3
        for j = 1:2
            m = classe == c & (T.compacte == (j == 1));
            if any(m), V(c, j) = 100 * mean(T.sig(m)); end
            Kt{c, j} = sprintf('%d/%d', sum(T.sig(m)), sum(m));
        end
    end
    hbd = bar(ax, 1:3, V, 'grouped', 'EdgeColor', 'none');
    hbd(1).FaceColor = BLEU; hbd(2).FaceColor = GRIS_C;
    for j = 1:2
        text(ax, hbd(j).XEndPoints, V(:, j)' + 0.8, Kt(:, j)', 'HorizontalAlignment', 'center', 'FontSize', 8);
    end
    yline(ax, 100 * seuil_p, '--', 'Color', ROUGE, 'LineWidth', 1.1);
    set(ax, 'XTick', 1:3, 'XTickLabel', {'<= 12', '13-25', '> 25'}, 'YLim', [0 max(32, max(V(:)) + 6)]);
    legend(ax, hbd, {sprintf('very compact (percentile <= %d)', seuil_compact), 'others'}, ...
        'Location', 'northwest', 'Box', 'off');
    xlabel(ax, 'colored cells in the assembly'); ylabel(ax, '% of assemblies with p < 0.05');
else
    manquant(ax, 'compacite\_pct');
end
titre_gauche(ax, 'D. Bias vs compactness, at equal size');
habiller(ax, 'y');

% ===== E. Bias vs compactness ============================================
ax = nexttile(tl); hold(ax, 'on');
if a_ratio || a_compac
    if a_ratio
        xc = T.ratio_etalement;
    else
        xc = sqrt(max(T.compacite_pct, 0));
    end
    hcat = scatter_categories(ax, xc, T.mlogp, 12 + 1.6 * T.n_colorees, cat_T, CATS);
    yline(ax, -log10(seuil_p), '--', 'Color', ENCRE2);
    if a_ratio
        xline(ax, 1, ':', 'Color', BLEU, 'LineWidth', 1);
        xlabel(ax, 'spread / random spread (1 = random; < 1 = clustered)');
    else
        xline(ax, sqrt(seuil_compact), ':', 'Color', BLEU, 'LineWidth', 1);
        tk = [0 1 5 10 25 50 100];
        set(ax, 'XTick', sqrt(tk), 'XTickLabel', string(tk), 'XLim', [-0.3 10.3]);
        xlabel(ax, 'compactness (percentile vs random groups; left = compact)');
    end
    legend(ax, flipud(hcat), flipud(CATS(:, 2)), 'Location', 'northeast', 'Box', 'off');
    text(ax, 0.98, 0.62, '(point size = no. of cells)', 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'FontSize', 8, 'Color', ENCRE2);
end
set(ax, 'YTick', 0:4, 'YTickLabel', PTICKS, 'YLim', [-0.1 4.3]);
ylabel(ax, 'minP p-value');
titre_gauche(ax, 'E. Bias vs compactness');
habiller(ax, 'y');

% ===== F. Bias vs number of colored cells ================================
ax = nexttile(tl); hold(ax, 'on');
scatter_categories(ax, T.n_colorees, T.mlogp, [], cat_T, CATS);
yline(ax, -log10(seuil_p), '--', 'Color', ENCRE2);
set(ax, 'XScale', 'log', 'XTick', [3 5 10 20 50 100], 'YTick', 0:4, 'YTickLabel', PTICKS, 'YLim', [-0.1 4.3]);
xlabel(ax, 'colored cells in the assembly (log)'); ylabel(ax, 'minP p-value');
titre_gauche(ax, 'F. Bias vs number of colored cells');
habiller(ax, 'y');

% ===== G-I. SCE descriptors by category ==================================
rng(graine_fig);
ax = nexttile(tl);
if all(ismember({'n_SCE_signif', 'n_SCE_total'}, T.Properties.VariableNames))
    strip_categories(ax, 100 * T.n_SCE_signif ./ T.n_SCE_total, cat_T, CATS, ...
        '% of the session''s SCEs', 'G. SCEs with significant activation (%)', ENCRE2);
else
    manquant(ax, 'n\_SCE\_signif / n\_SCE\_total');
end
ax = nexttile(tl);
if ismember('n_SCE_signif', T.Properties.VariableNames)
    strip_categories(ax, T.n_SCE_signif, cat_T, CATS, 'number of SCEs', ...
        'H. SCEs with significant activation (n)', ENCRE2);
else
    manquant(ax, 'n\_SCE\_signif');
end
ax = nexttile(tl);
if ismember('frac_active_moy', T.Properties.VariableNames)
    strip_categories(ax, T.frac_active_moy, cat_T, CATS, 'fraction of assembly cells active per SCE', ...
        'I. Recruitment within its SCEs', ENCRE2);
else
    manquant(ax, 'frac\_active\_moy');
end

% ===== J. Recruitment vs assembly size ===================================
% Le recrutement baisse quand la taille augmente (un SCE ne mobilise
% qu'un nombre limité de cellules) : on compare les biaisées à la tendance.
ax = nexttile(tl); hold(ax, 'on');
if ismember('frac_active_moy', T.Properties.VariableNames)
    scatter_categories(ax, T.n_total, T.frac_active_moy, [], cat_T, CATS);
    lt = log(T.n_total); yv = T.frac_active_moy; okr = ~isnan(yv) & ~isnan(lt);
    cf = polyfit(lt(okr), yv(okr), 1); xs = linspace(min(lt(okr)), max(lt(okr)), 50);
    plot(ax, exp(xs), polyval(cf, xs), '--', 'Color', ENCRE2, 'LineWidth', 1.2);
    residus = yv - polyval(cf, lt);
    est_b = strcmp(cat_T, 'biased') & okr;
    if sum(est_b) > 1 && sum(~est_b & okr) > 1
        p_res = ranksum(residus(est_b), residus(~est_b & okr));
        text(ax, 0.98, 0.97, {'dashed = size trend', sprintf('biased vs trend: p = %.2f', p_res)}, ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
            'FontSize', 8, 'Color', ENCRE2);
    end
    set(ax, 'XScale', 'log', 'XTick', [10 20 50 100]);
    xlabel(ax, 'cells in the assembly (log)'); ylabel(ax, 'fraction of assembly cells active per SCE');
else
    manquant(ax, 'frac\_active\_moy');
end
titre_gauche(ax, 'J. Recruitment vs assembly size');
habiller(ax, 'y');

% ===== K-L. Assembly size ================================================
ax = nexttile(tl);
if ismember('pct_reseau', T.Properties.VariableNames)
    strip_categories(ax, T.pct_reseau, cat_T, CATS, '% of the field''s cells', 'K. Assembly size (%)', ENCRE2);
else
    manquant(ax, 'pct\_reseau');
end
ax = nexttile(tl);
strip_categories(ax, T.n_total, cat_T, CATS, 'number of cells', 'L. Assembly size (n)', ENCRE2);

exportgraphics(fig, fullfile(dossier_racine, 'clonal_bias_complementary_analyses.png'), 'Resolution', 200);

%% ---------- Résumé dans la console --------------------------------------
fprintf('\n=== Analyses complémentaires ===\n');
fprintf('Proportion à p < %.2f : %.1f %% (%d/%d) ; IC 95 %% (animaux) : %.1f à %.1f %%\n', ...
    seuil_p, 100 * p_glob, sum(T.sig), height(T), 100 * ic_boot);
if ~isnan(glmm(1)), fprintf('GLMM : %.1f %% ; IC 95 %% : %.1f à %.1f %%\n', 100 * glmm); end
fprintf('p < 0,25 : %.1f %% des assemblées (25 %% attendus) ; test KS : p = %.3f\n', ...
    100 * mean(T.p_minp < 0.25), p_ks);
for i = 1:numel(anim_u)
    fprintf('  sans %-8s : %.1f %%\n', anim_u{i}, 100 * p_sans(i));
end
if ismember('K', S.R.Properties.VariableNames) && any(~isnan(S.R.K))
    [~, iu] = unique(S.R.session_u); Ks = S.R.K(iu);
    fprintf('Sessions au plus petit K retenu (K = %d) : %d / %d\n', min(Ks), sum(Ks == min(Ks)), numel(Ks));
end


%% ---------- Fonctions locales -------------------------------------------
function [lo, hi] = wilson(k, n)
% IC à 95 % de Wilson pour une proportion k/n
z = 1.96; p = k ./ n; den = 1 + z^2 ./ n;
c = (p + z^2 ./ (2 * n)) ./ den;
d = z * sqrt(p .* (1 - p) ./ n + z^2 ./ (4 * n.^2)) ./ den;
lo = max(0, c - d); hi = min(1, c + d);
end

function h = scatter_categories(ax, x, y, tailles, cat_T, CATS)
% Nuage de points coloré par catégorie ; tailles = [] -> taille fixe
hold(ax, 'on');
h = gobjects(size(CATS, 1), 1);
for k = 1:size(CATS, 1)
    m = strcmp(cat_T, CATS{k, 1});
    if isempty(tailles)
        if strcmp(CATS{k, 1}, 'biased'), s = 36; else, s = 20; end
    else
        s = tailles(m);
    end
    h(k) = scatter(ax, x(m), y(m), s, CATS{k, 4}, 'filled', 'MarkerFaceAlpha', CATS{k, 5}, 'MarkerEdgeColor', 'w');
end
end

function strip_categories(ax, valeurs, cat_T, CATS, ylab, titre, couleur_txt)
% Une colonne de points par catégorie, médiane en trait noir, Kruskal-Wallis
hold(ax, 'on');
v_all = []; g_all = []; etiq = cell(1, size(CATS, 1));
for k = 1:size(CATS, 1)
    v = valeurs(strcmp(cat_T, CATS{k, 1})); v = v(~isnan(v));
    if strcmp(CATS{k, 1}, 'biased'), s = 36; else, s = 18; end
    scatter(ax, (k - 1) + (rand(size(v)) - 0.5) * 0.34, v, s, CATS{k, 4}, 'filled', ...
        'MarkerFaceAlpha', CATS{k, 5}, 'MarkerEdgeColor', 'w');
    if ~isempty(v), plot(ax, (k - 1) + [-0.3 0.3], median(v) * [1 1], '-k', 'LineWidth', 2); end
    v_all = [v_all; v]; g_all = [g_all; k * ones(numel(v), 1)]; %#ok<AGROW>
    etiq{k} = sprintf('%s\\newline n = %d', CATS{k, 3}, numel(v));
end
if numel(unique(g_all)) > 1
    pkw = kruskalwallis(v_all, g_all, 'off');
    if pkw < 0.001, ptxt = 'p < 0.001'; else, ptxt = sprintf('p = %.3f', pkw); end
    text(ax, 0.98, 0.97, {['Kruskal-Wallis ' ptxt], '(bar = median)'}, 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 8, 'Color', couleur_txt);
end
set(ax, 'XTick', 0:size(CATS, 1) - 1, 'XTickLabel', etiq, 'XLim', [-0.6 size(CATS, 1) - 0.4], ...
    'TickLabelInterpreter', 'tex', 'XTickLabelRotation', 0);
ylabel(ax, ylab);
titre_gauche(ax, titre);
habiller(ax, 'y');
end

function manquant(ax, nom)
text(ax, 0.5, 0.5, [nom ' missing'], 'HorizontalAlignment', 'center');
axis(ax, 'off');
end

function titre_gauche(ax, txt)
title(ax, txt, 'FontWeight', 'bold', 'FontSize', 10.5);
try, ax.TitleHorizontalAlignment = 'left'; catch, end   % R2020b et plus
end

function habiller(ax, grille)
box(ax, 'off');
ax.TickDir = 'out';
ax.FontSize = 9;
ax.GridColor = [0.85 0.85 0.82]; ax.GridAlpha = 1;
if strcmp(grille, 'x'), ax.XGrid = 'on'; else, ax.YGrid = 'on'; end
end
