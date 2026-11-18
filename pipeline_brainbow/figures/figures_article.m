function figures_article(dossier_variante, fdr, alpha, nom_variante)
%FIGURES_ARTICLE  Deux figures pour l'article, à partir des résultats de
%   l'étape 4 d'une variante (aucun calcul nouveau, sauf l'IC des animaux).
%
%   figures_article('E:\...\agregation\champ_globale')
%   (lancée aussi par l'étape 4 pour chaque variante)
%
%   Écrit dans dossier_variante, en PNG (300 dpi) et en PDF vectoriel (pour
%   le montage) :
%   assembly_sizes_by_animal   (descripteurs SANS couleur)
%     a. assemblées par session, par animal
%     b. taille des assemblées (cellules), par animal ; pleins = testables
%     c. taille (% des cellules du champ), par animal
%     d. cellules colorées par assemblée, par animal (testabilité)
%   clonal_bias_quantification
%     a. par animal : % d'assemblées à p_minp < alpha, observé contre hasard
%        exact (taux_nul) ; ensemble : IC par rééchantillonnage des animaux ;
%        statistiques de l'étape 4 (bilan)
%     b. par animal : assemblées à p < alpha, colorées par la couleur en
%        excès (plein = biaisée après FDR)
%     c. force du biais : part de la couleur en excès dans l'assemblée
%        contre sa part dans le champ (diagonale = pas d'excès)
%     d. couleur en excès des assemblées à p < alpha, et composition du champ
%     e. par taille (cellules colorées) : observé contre hasard exact
%     f. test principal contre l'autre variante (même assemblée), si elle existe

if nargin < 2 || isempty(fdr), fdr = 0.10; end
if nargin < 3 || isempty(alpha), alpha = 0.05; end
dossier_variante = regexprep(char(dossier_variante), '[\\/]+$', '');
[parent, nom_dossier] = fileparts(dossier_variante);
if nargin < 4 || isempty(nom_variante), nom_variante = nom_dossier; end
S = load(fullfile(dossier_variante, 'resultats_biais_clonal.mat'));
R = S.R;
bilan = struct();
if isfield(S, 'bilan'), bilan = S.bilan; end

% ---------- style ----------
ROUGE = [200 40 62] / 255;  AMBRE = [227 161 59] / 255;  GRIS = [0.62 0.62 0.60];
GRIS_C = [0.85 0.85 0.83];  ENCRE = [0.15 0.15 0.15];   BLEU = [46 109 180] / 255;
noms_coul = {'R', 'G', 'B', 'RG', 'RB', 'GB'};
cols = [0.86 0.16 0.16; 0.16 0.66 0.22; 0.20 0.38 0.92; 0.93 0.78 0.10; 0.82 0.22 0.80; 0.10 0.74 0.80];
clair = 0.45 * cols + 0.55;
police = 'Arial';
rng(1);

% ---------- données ----------
R.cle_session = strcat(R.animal, '_', R.session, '_', R.champ);
anim_u = unique(R.animal);
nA = numel(anim_u);
[~, ia] = ismember(R.animal, anim_u);
T = R.testable;
sig = T & R.p_minp < alpha;
bia = T & R.q < fdr;
[~, ic_coul] = ismember(R.couleur_dominante, noms_coul);
taux = R.taux_nul; taux(isnan(taux)) = 0;

%% ===================== FIGURE 1 : tailles par animal =====================
f1 = figure('Color', 'w', 'Position', [40 40 1250 820]);
t1 = tiledlayout(f1, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t1, sprintf('Assemblies (no colour used): %d assemblies, %d sessions, %d animals', ...
    height(R), numel(unique(R.cle_session)), nA), 'FontName', police, 'FontSize', 11, 'FontWeight', 'bold');

% a. assemblées par session
ax = nexttile(t1); hold(ax, 'on');
[sess_u, i_s] = unique(R.cle_session);
n_par_sess = cellfun(@(k) sum(strcmp(R.cle_session, k)), sess_u);
a_sess = ia(i_s);
for a = 1:nA
    v = n_par_sess(a_sess == a);
    scatter(ax, a + jitter(numel(v), 0.18), v, 28, BLEU, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', 'w');
    if ~isempty(v), plot(ax, a + [-0.28 0.28], median(v) * [1 1], '-', 'Color', ENCRE, 'LineWidth', 2); end
end
habiller(ax, anim_u, police);
ylabel(ax, 'assemblies per session');
lettre(ax, 'a', 'Assemblies per session', police);

% b. taille (cellules)
ax = nexttile(t1); hold(ax, 'on');
strip_animal(ax, R.n_total, ia, T, nA, BLEU, ENCRE);
set(ax, 'YScale', 'log', 'YTick', [5 10 20 50 100 200]);
habiller(ax, anim_u, police);
ylabel(ax, 'cells in the assembly');
lettre(ax, 'b', 'Assembly size (filled: testable)', police);

% c. taille (% du champ)
ax = nexttile(t1); hold(ax, 'on');
if ismember('pct_reseau', R.Properties.VariableNames)
    strip_animal(ax, R.pct_reseau, ia, T, nA, BLEU, ENCRE);
    set(ax, 'YScale', 'log', 'YTick', [0.5 1 2 5 10 20 50]);
end
habiller(ax, anim_u, police);
ylabel(ax, '% of the field''s cells');
lettre(ax, 'c', 'Assembly size relative to the field', police);

% d. cellules colorées
ax = nexttile(t1); hold(ax, 'on');
strip_animal(ax, max(R.n_colorees, 0.8), ia, T, nA, BLEU, ENCRE);
set(ax, 'YScale', 'log', 'YTick', [1 2 5 10 20 50 100]);
habiller(ax, anim_u, police);
for a = 1:nA
    m = ia == a;
    text(ax, a, 1, sprintf('%d/%d', sum(T & m), sum(m)), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'Units', 'data', 'FontSize', 7, 'Color', ENCRE, 'FontName', police);
end
ylabel(ax, 'coloured cells in the assembly');
lettre(ax, 'd', 'Coloured cells (numbers: testable / all)', police);
exporter(f1, fullfile(dossier_variante, 'assembly_sizes_by_animal'));

%% ===================== FIGURE 2 : quantification du biais =================
f2 = figure('Color', 'w', 'Position', [40 40 1500 900]);
t2 = tiledlayout(f2, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t2, {sprintf('Clonal bias of assemblies: %d testable assemblies, %d animals', sum(T), ...
    numel(unique(R.animal(T)))), texte_variante(nom_variante)}, 'FontName', police, 'FontSize', 11, ...
    'FontWeight', 'bold', 'Interpreter', 'none');

% a. observé contre hasard, par animal
ax = nexttile(t2); hold(ax, 'on');
obs_a = nan(nA, 1); att_a = nan(nA, 1); n_a = zeros(nA, 1);
for a = 1:nA
    m = T & ia == a;
    n_a(a) = sum(m);
    if n_a(a) == 0, continue; end
    obs_a(a) = 100 * mean(sig(m));
    att_a(a) = 100 * mean(taux(m));
    plot(ax, [1 2], [att_a(a) obs_a(a)], '-', 'Color', [GRIS 0.8], 'LineWidth', 1);
    scatter(ax, [1 2], [att_a(a) obs_a(a)], 16 + 2 * n_a(a), GRIS, 'filled', 'MarkerEdgeColor', 'w');
    text(ax, 2.08, obs_a(a), anim_u{a}, 'FontSize', 7, 'Color', ENCRE, 'FontName', police);
end
% ensemble : IC en rééchantillonnant les animaux
ok_a = find(n_a > 0);
boot = zeros(5000, 1);
for b = 1:5000
    tir = ok_a(randi(numel(ok_a), numel(ok_a), 1));
    k = arrayfun(@(x) sum(sig(T & ia == x)), tir);
    n = arrayfun(@(x) sum(T & ia == x), tir);
    boot(b) = 100 * sum(k) / max(sum(n), 1);
end
obs_tot = 100 * mean(sig(T)); att_tot = 100 * mean(taux(T)); ic = prctile(boot, [2.5 97.5]);
plot(ax, [1 2], [att_tot obs_tot], '-', 'Color', ROUGE, 'LineWidth', 2.2);
plot(ax, [2 2], ic, '-', 'Color', ROUGE, 'LineWidth', 2.2);
scatter(ax, [1 2], [att_tot obs_tot], 70, ROUGE, 'd', 'filled', 'MarkerEdgeColor', 'w');
txt = {sprintf('all: %.1f%% vs %.1f%% by chance', obs_tot, att_tot), ...
    sprintf('95%% CI (animals): %.1f-%.1f%%', ic)};
if isfield(bilan, 'OR') && ~isnan(bilan.OR)
    txt{end+1} = sprintf('GLMM OR %.2f [%.2f-%.2f], p = %s', bilan.OR, bilan.OR_ic, ptexte(bilan.p_glmm));
end
if isfield(bilan, 'p_wilcoxon'), txt{end+1} = sprintf('Wilcoxon (animals) p = %s', ptexte(bilan.p_wilcoxon)); end
if isfield(bilan, 'p_wil_score') && ~isnan(bilan.p_wil_score)
    txt{end+1} = sprintf('all p-values (animals) p = %s', ptexte(bilan.p_wil_score));
end
if isfield(bilan, 'pi1') && ~isnan(bilan.pi1)
    txt{end+1} = sprintf('clonal assemblies >= %.0f%% [%.0f-%.0f]', 100 * bilan.pi1, 100 * bilan.pi1_ic);
end
text(ax, 0.03, 0.98, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 7.5, ...
    'Color', ENCRE, 'FontName', police);
set(ax, 'XTick', [1 2], 'XTickLabel', {'chance (exact)', 'observed'}, 'XLim', [0.6 2.6], ...
    'YLim', [0 max([obs_a; ic(2); 10], [], 'omitnan') * 1.15]);
style(ax, police);
ylabel(ax, sprintf('%% of testable assemblies with p_{minP} < %g', alpha));
lettre(ax, 'a', 'Observed vs chance, per animal', police);

% b. par animal, par couleur
ax = nexttile(t2); hold(ax, 'on');
Mb = zeros(nA, 12);
for a = 1:nA
    for c = 1:6
        Mb(a, c)     = sum(bia & ia == a & ic_coul == c);
        Mb(a, 6 + c) = sum(sig & ~bia & ia == a & ic_coul == c);
    end
end
hb = bar(ax, 1:nA, Mb, 'stacked', 'EdgeColor', 'none', 'BarWidth', 0.65);
for c = 1:6
    hb(c).FaceColor = cols(c, :); hb(6 + c).FaceColor = clair(c, :);
end
for a = 1:nA
    text(ax, a, sum(Mb(a, :)) + 0.15, sprintf('%d/%d', sum(Mb(a, :)), n_a(a)), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 7, 'FontName', police, 'Color', ENCRE);
end
habiller(ax, anim_u, police);
ylim(ax, [0 max(1, max(sum(Mb, 2))) + 1.2]);
ylabel(ax, sprintf('assemblies with p_{minP} < %g', alpha));
legend(ax, [hb(1:6) hb(7)], [noms_coul, {sprintf('light: not after FDR %g%%', 100 * fdr)}], ...
    'Location', 'northeast', 'Box', 'off', 'FontSize', 7, 'NumColumns', 2);
lettre(ax, 'b', 'Biased assemblies per animal (/ testable)', police);

% c. force du biais
ax = nexttile(t2); hold(ax, 'on');
a_kcol = ismember('Kcol_urne', R.Properties.VariableNames);
if a_kcol
    % référence de chaque assemblée (champ, ou sa surface pour 'territoire')
    Kref = R.Kcol_urne;
    if ismember('Kcol_ref', R.Properties.VariableNames), Kref = R.Kcol_ref; end
    part_champ = nan(height(R), 1);
    for i = 1:height(R)
        if ic_coul(i) > 0 && sum(Kref(i, :)) > 0
            part_champ(i) = Kref(i, ic_coul(i)) / sum(Kref(i, :));
        end
    end
    groupes = {T & ~sig, GRIS, 14; sig & ~bia, AMBRE, 24; bia, ROUGE, 40};
    for g = 1:3
        m = groupes{g, 1} & ~isnan(part_champ);
        scatter(ax, 100 * part_champ(m), 100 * R.part_couleur_dominante(m), groupes{g, 3}, groupes{g, 2}, ...
            'filled', 'MarkerEdgeColor', 'w', 'MarkerFaceAlpha', 0.85);
    end
    plot(ax, [0 100], [0 100], '--', 'Color', ENCRE);
    axis(ax, [0 60 0 100]);
    legend(ax, {sprintf('p >= %g', alpha), sprintf('p < %g', alpha), sprintf('biased (FDR %g%%)', 100 * fdr)}, ...
        'Location', 'southeast', 'Box', 'off', 'FontSize', 7);
    xlabel(ax, '% of the reference''s coloured cells, that colour (field or own area)');
    ylabel(ax, '% of the assembly''s coloured cells');
else
    text(ax, 0.5, 0.5, {'field composition missing', '(run step 4 again)'}, 'Units', 'normalized', ...
        'HorizontalAlignment', 'center');
end
style(ax, police);
lettre(ax, 'c', 'Strength of the bias (colour in excess)', police);

% d. couleurs
ax = nexttile(t2); hold(ax, 'on');
n_sig_c = arrayfun(@(c) sum(sig & ic_coul == c), 1:6);
n_bia_c = arrayfun(@(c) sum(bia & ic_coul == c), 1:6);
comp = nan(1, 6);
if a_kcol
    [~, i_u] = unique(R.cle_session(T));
    Kt = R.Kcol_urne(T, :);
    Kt = Kt(i_u, :);
    comp = 100 * sum(Kt, 1, 'omitnan') / sum(Kt(:), 'omitnan');
end
yy = [100 * n_sig_c / max(sum(n_sig_c), 1); comp]';
hd = bar(ax, 1:6, yy, 'grouped', 'EdgeColor', 'none');
hd(1).FaceColor = 'flat'; hd(1).CData = cols;
hd(2).FaceColor = GRIS_C;
for c = 1:6
    text(ax, hd(1).XEndPoints(c), yy(c, 1) + 1, sprintf('%d (%d)', n_sig_c(c), n_bia_c(c)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7, 'FontName', police);
end
set(ax, 'XTick', 1:6, 'XTickLabel', noms_coul);
style(ax, police);
ylabel(ax, '% ');
legend(ax, hd, {sprintf('colour in excess, assemblies p < %g', alpha), 'field composition (coloured cells)'}, ...
    'Location', 'northeast', 'Box', 'off', 'FontSize', 7);
lettre(ax, 'd', 'Which colour (numbers: p < alpha (FDR))', police);

% e. par taille
ax = nexttile(t2); hold(ax, 'on');
bornes = [0 12 25 Inf];
etiq = {'3-12', '13-25', '> 25'};
obs_t = nan(1, 3); att_t = nan(1, 3); nn = zeros(1, 3);
for k = 1:3
    m = T & R.n_colorees > bornes(k) & R.n_colorees <= bornes(k + 1);
    nn(k) = sum(m);
    if nn(k) > 0
        obs_t(k) = 100 * mean(sig(m)); att_t(k) = 100 * mean(taux(m));
    end
end
he = bar(ax, 1:3, [att_t; obs_t]', 'grouped', 'EdgeColor', 'none');
he(1).FaceColor = GRIS_C; he(2).FaceColor = ROUGE;
for k = 1:3
    text(ax, he(2).XEndPoints(k), max(obs_t(k), 0) + 0.8, sprintf('%d/%d', sum(sig & T & ...
        R.n_colorees > bornes(k) & R.n_colorees <= bornes(k + 1)), nn(k)), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 7, 'FontName', police);
end
set(ax, 'XTick', 1:3, 'XTickLabel', etiq);
style(ax, police);
xlabel(ax, 'coloured cells in the assembly');
ylabel(ax, sprintf('%% with p_{minP} < %g', alpha));
legend(ax, he, {'chance (exact)', 'observed'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 7);
lettre(ax, 'e', 'By assembly size (chance depends on size)', police);

% f. autre variante
ax = nexttile(t2); hold(ax, 'on');
autre = autre_variante(parent, nom_dossier);
if ~isempty(autre)
    S2 = load(fullfile(parent, autre, 'resultats_biais_clonal.mat'), 'R');
    R2 = S2.R;
    k1 = strcat(R.cle_session, '#', arrayfun(@(x) sprintf('%d', x), R.assemblee, 'UniformOutput', false));
    k2 = strcat(R2.animal, '_', R2.session, '_', R2.champ, '#', ...
        arrayfun(@(x) sprintf('%d', x), R2.assemblee, 'UniformOutput', false));
    [ok, j] = ismember(k1, k2);
    m = ok & T;
    x = -log10(R.p_minp(m)); y = -log10(R2.p_minp(j(m)));
    b1 = bia(m); s1 = sig(m);
    scatter(ax, x(~s1), y(~s1), 14, GRIS, 'filled', 'MarkerEdgeColor', 'w');
    scatter(ax, x(s1 & ~b1), y(s1 & ~b1), 24, AMBRE, 'filled', 'MarkerEdgeColor', 'w');
    scatter(ax, x(b1), y(b1), 40, ROUGE, 'filled', 'MarkerEdgeColor', 'w');
    lim = max([x; y; 2]) * 1.05;
    plot(ax, [0 lim], [0 lim], '--', 'Color', ENCRE);
    xline(ax, -log10(alpha), ':', 'Color', ENCRE); yline(ax, -log10(alpha), ':', 'Color', ENCRE);
    axis(ax, [0 lim 0 lim]);
    tk = 0:ceil(lim);
    etq = cellstr(compose('10^{-%d}', tk));
    set(ax, 'XTick', tk, 'XTickLabel', etq, 'YTick', tk, 'YTickLabel', etq);
    xlabel(ax, sprintf('p_{minP}, %s', court(nom_dossier)));
    ylabel(ax, sprintf('p_{minP}, %s', court(autre)));
    text(ax, 0.03, 0.97, sprintf('p < %g in both: %d of %d', alpha, sum(x > -log10(alpha) & y > -log10(alpha)), ...
        sum(s1)), 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 7.5, 'FontName', police);
else
    text(ax, 0.5, 0.5, 'no other variant', 'Units', 'normalized', 'HorizontalAlignment', 'center');
end
style(ax, police);
lettre(ax, 'f', 'Same assemblies, other variant', police);
exporter(f2, fullfile(dossier_variante, 'clonal_bias_quantification'));
end


%% ---------------------------------------------------------------------------
function strip_animal(ax, v, ia, testable, nA, BLEU, ENCRE)
for a = 1:nA
    m = ia == a & ~isnan(v);
    x = a + jitter(sum(m), 0.22);
    vt = v(m); tt = testable(m);
    scatter(ax, x(~tt), vt(~tt), 14, BLEU, 'MarkerEdgeAlpha', 0.6);
    scatter(ax, x(tt), vt(tt), 14, BLEU, 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'none');
    if any(m)
        q = prctile(vt, [25 50 75]);
        plot(ax, a + [0.3 0.3], q([1 3]), '-', 'Color', ENCRE, 'LineWidth', 1.2);
        plot(ax, a + [-0.3 0.3], q(2) * [1 1], '-', 'Color', ENCRE, 'LineWidth', 2);
    end
end
end

function x = jitter(n, largeur)
x = (rand(n, 1) - 0.5) * 2 * largeur;
end

function habiller(ax, anim_u, police)
set(ax, 'XTick', 1:numel(anim_u), 'XTickLabel', anim_u, 'XLim', [0.4 numel(anim_u) + 0.6], ...
    'XTickLabelRotation', 30, 'TickLabelInterpreter', 'none');
style(ax, police);
end

function style(ax, police)
box(ax, 'off');
ax.TickDir = 'out';
ax.FontName = police;
ax.FontSize = 8;
ax.YGrid = 'on';
ax.GridColor = [0.88 0.88 0.86]; ax.GridAlpha = 1;
end

function lettre(ax, l, titre, police)
title(ax, ['\bf' l '\rm   ' titre], 'FontName', police, 'FontSize', 9.5, 'FontWeight', 'normal');
try, ax.TitleHorizontalAlignment = 'left'; catch, end
end

function exporter(f, base)
exportgraphics(f, [base '.png'], 'Resolution', 300);
try
    exportgraphics(f, [base '.pdf'], 'ContentType', 'vector');
catch err
    warning('PDF non écrit (%s) : %s', base, err.message);
end
end

function t = ptexte(p)
if isempty(p) || isnan(p), t = 'NA'; elseif p < 0.001, t = '< 0.001'; else, t = sprintf('%.3f', p); end
end

function a = autre_variante(parent, nom)
% l'autre variante à comparer : le contrôle spatial pour le test principal, sinon le test principal
a = '';
L = dir(parent);
noms = {L([L.isdir]).name};
noms = noms(~ismember(noms, {'.', '..', nom}));
noms = noms(cellfun(@(n) isfile(fullfile(parent, n, 'resultats_biais_clonal.mat')), noms));
if isempty(noms), return; end
% test principal -> contrôle « territoire » (sinon les carrés) ; un contrôle -> test principal
if startsWith(nom, 'champ_globale')
    pref = {'champ_territoire', 'champ_blocs'};
else
    pref = {'champ_globale'};
end
for k = 1:numel(pref)
    i = find(startsWith(noms, pref{k}), 1);
    if ~isempty(i), a = noms{i}; return; end
end
a = noms{1};
end

function t = court(nom)
nom = regexprep(nom, '_PROVISOIRE$', '');
switch nom
    case 'champ_globale', t = 'main test';
    case 'champ_blocs',   t = 'spatial control';
    case 'champ_activite', t = 'activity control';
    case 'champ_territoire', t = 'territory control';
    otherwise,            t = strrep(nom, '_', ' ');
end
end

function t = texte_variante(nom)
principal = contains(nom, 'PRINCIPAL');
base = regexprep(strtrim(regexprep(nom, '\(.*\)', '')), '_PROVISOIRE$', '');
switch base
    case 'champ_globale',  t = 'field reference, colours shuffled across the field';
    case 'champ_blocs',    t = 'spatial control: colours shuffled only between nearby cells';
    case 'champ_activite', t = 'activity control: colours shuffled between cells of similar activity';
    case 'champ_territoire', t = 'territory control: reference = coloured cells within the assembly''s own area';
    otherwise,             t = base;
end
if principal, t = ['MAIN TEST: ' t]; end
end
