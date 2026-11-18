function fig = figure_paires_resultats(dossier, rayon_um, mesure)
%FIGURE_PAIRES_RESULTATS  Figure de l'étape 3b (paires), sans rien recalculer :
%   lit paires_agregation.mat (écrit par agregation_paires).
%
%   figure_paires_resultats('E:\...\paires_agregation')              % test principal
%   figure_paires_resultats('E:\...\paires_agregation', 100)         % paires <= 100 µm
%   figure_paires_resultats('E:\...\paires_agregation', 100, 'dff')  % idem, Pearson du dF/F
%   Un autre rayon ou une autre mesure que ceux du test principal : la figure
%   est écrite sous paires_<mesure>_<rayon>um.png et le titre le signale
%   (exploration, le test principal reste celui de config.paires).
%
%   Mesure et permutation principales (config.paires), l'animal est l'unité.
%   A1. paires du test principal (à moins de R µm), toutes couleurs : pour
%      chaque animal, similarité (STTC corrigé de la distance et du niveau de
%      chaque cellule) des paires de même couleur et de couleurs différentes,
%      reliées (gris) ; noir = moyenne des animaux +/- sem (même présentation
%      que l'ancienne figure G, mais à distance égale et par animal).
%   A1b. le même graphe avec la mesure BRUTE (sans résidus : ni distance ni
%       niveau de chaque cellule retirés) ; Wilcoxon sur les animaux.
%   A2. le même graphe avec le Pearson du dF/F (ou le STTC si la figure est
%       en Pearson).
%   A3. différence (même - différente) de chaque session, rangée par animal
%       (point plein : p de permutation < 0,05), trait = moyenne de l'animal.
%   B. différence même - différente selon la distance (disques cumulés ou
%      tranches) : un trait fin par animal, moyenne +/- sem ; p = Wilcoxon
%      unilatéral sur les animaux. Bande grise : le test principal. Dans le
%      titre : le modèle mixte sur toutes les distances (anneaux disjoints).
%   C. chaque couleur au rayon principal : un point par animal.
%   D. regroupement spatial : paires de même couleur observées / attendues
%      (couleurs mélangées dans le champ), selon la distance.
%   E. sensibilité au rayon principal : autres mesures (STTC à 0,1 et 0,5 s,
%      corrélation de Pearson du dF/F) et autres permutations, en % de la
%      similarité moyenne des paires.
%   Écrit paires_agregation.png.

dossier = regexprep(char(dossier), '[\\/]+$', '');
S = load(fullfile(dossier, 'paires_agregation.mat'));
A = S.A; Ani = S.Ani; opts = S.opts;
mes0 = opts.mesures{1};
perm0 = opts.permutations{1};
pd = S.principal_dist;
pd_principal = pd;
if nargin >= 2 && ~isempty(rayon_um)
    if S.cumulatif
        pd = [0, rayon_um];
    else
        warning('Distances en tranches : rayon ignoré.');
    end
end
if nargin >= 3 && ~isempty(mesure), mes0 = char(mesure); end
est_principal = isequal(pd, pd_principal) && strcmp(mes0, opts.mesures{1});
if ~any(strcmp(A.mesure, mes0) & A.dist_min_um == pd(1) & A.dist_max_um == pd(2))
    error('Pas de résultat pour %s à %g-%g µm dans %s.', mes0, pd(1), pd(2), dossier);
end
noms = {'R', 'G', 'B', 'RG', 'RB', 'GB'};
cols = [0.86 0.16 0.16; 0.16 0.66 0.22; 0.20 0.38 0.92; ...
        0.93 0.78 0.10; 0.82 0.22 0.80; 0.10 0.74 0.80];
GRIS = [0.7 0.7 0.7];
BLEU = [0.15 0.39 0.65];

% distances, dans l'ordre
Ad = A(strcmp(A.mesure, mes0) & strcmp(A.permutation, perm0) & strcmp(A.couleur, 'toutes'), :);
Ad = sortrows(Ad, {'dist_min_um', 'dist_max_um'});
dist = [Ad.dist_min_um, Ad.dist_max_um];
nd = size(dist, 1);
if S.cumulatif
    etiq = arrayfun(@(r) sprintf('<= %g', r), dist(:, 2), 'UniformOutput', false);
    etiq = strrep(etiq, '<= Inf', 'all');
    texte_d = sprintf('pairs <= %g um', pd(2));
    xlab = 'pairs closer than ... um (cumulative)';
else
    etiq = arrayfun(@(i) sprintf('%g-%g', dist(i, 1), dist(i, 2)), (1:nd)', 'UniformOutput', false);
    etiq = strrep(etiq, '-Inf', '+');
    texte_d = sprintf('pairs %g-%g um', pd(1), pd(2));
    xlab = 'distance between cells (um)';
end
i0 = find(dist(:, 1) == pd(1) & dist(:, 2) == pd(2), 1);
if S.cumulatif
    texte_t = sprintf('\\leq %g \\mum', pd(2));       % pour les titres (TeX)
else
    texte_t = sprintf('%g-%g \\mum', pd(1), pd(2));
end
sel_ani = @(m, p, c, d) Ani(strcmp(Ani.mesure, m) & strcmp(Ani.permutation, p) & strcmp(Ani.couleur, c) & ...
    Ani.dist_min_um == d(1) & Ani.dist_max_um == d(2), :);
r0 = Ad(i0, :);

fig = figure('Color', 'w', 'Position', [20 40 2500 920]);
tl = tiledlayout(fig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
if est_principal
    texte_expl = '';
else
    texte_expl = sprintf(' | EXPLORATION (main test: %s, pairs <= %g um)', nom_mesure(opts.mesures{1}), pd_principal(2));
end
title(tl, sprintf('Same-colour vs different-colour pairs at rest (%d animals, %d sessions)%s', ...
    S.n_animaux, S.n_sessions, texte_expl), 'FontWeight', 'bold', 'FontSize', 12);

%% A. même couleur contre couleurs différentes, par animal
ax = nexttile(tl);
panneau_paires_meme_diff(ax, S, mes0, perm0, pd, 'A1. ', 'residu');
% A1b : la même mesure brute (sans résidus : ni distance ni niveau des cellules retirés)
ax = nexttile(tl);
panneau_paires_meme_diff(ax, S, mes0, perm0, pd, 'A1b. ', 'brut');
% A2 : le même graphe avec l'autre mesure (Pearson du dF/F si la figure est en STTC)
if strcmp(mes0, 'dff'), mes_autre = opts.mesures{1}; else, mes_autre = 'dff'; end
if strcmp(mes_autre, mes0), mes_autre = 'sttc'; end
ax = nexttile(tl);
panneau_paires_meme_diff(ax, S, mes_autre, perm0, pd, 'A2. ', 'residu');
T0 = sel_ani(mes0, perm0, 'toutes', pd);

%% A2. différence de chaque session, rangée par animal (comme la figure du
% regroupement spatial des couleurs) : point = session (plein : p de
% permutation < 0,05), trait noir = moyenne de l'animal
ax = nexttile(tl); hold(ax, 'on');
if ~isempty(T0)
    animaux_o = sort(T0.animal);
    nA = numel(animaux_o);
    n_sig = 0; n_s = 0;
    rng(2);
    for a = 1:nA
        if isfield(S, 'Rsess')
            Rs = S.Rsess(strcmp(S.Rsess.mesure, mes0) & strcmp(S.Rsess.permutation, perm0) & ...
                strcmp(S.Rsess.couleur, 'toutes') & S.Rsess.dist_min_um == pd(1) & ...
                S.Rsess.dist_max_um == pd(2) & strcmp(S.Rsess.animal, animaux_o{a}) & S.Rsess.testable, :);
            x = a + 0.16 * (rand(height(Rs), 1) - 0.5) * 2;
            sig = Rs.p < 0.05;
            scatter(ax, x(~sig), Rs.delta(~sig), 36, [0.45 0.45 0.45], 'LineWidth', 1);
            scatter(ax, x(sig), Rs.delta(sig), 36, [0.45 0.45 0.45], 'filled');
            n_sig = n_sig + sum(sig); n_s = n_s + height(Rs);
        end
        da = T0.delta(strcmp(T0.animal, animaux_o{a}));
        plot(ax, a + [-0.25 0.25], da(1) * [1 1], 'k-', 'LineWidth', 2);
    end
    yline(ax, 0, '--', 'Color', [0.5 0.5 0.5]);
    set(ax, 'XTick', 1:nA, 'XTickLabel', animaux_o, 'XLim', [0.4 nA + 0.6], ...
        'TickLabelInterpreter', 'none', 'FontSize', 9, 'TickDir', 'out');
    xtickangle(ax, 30);
    ylabel(ax, sprintf('\\Delta %s (same - different)', nom_mesure(mes0)));
    d_a = T0.delta;
    if n_s > 0
        ligne2 = sprintf('%d/%d sessions p < 0.05', n_sig, n_s);
    else
        ligne2 = '';
    end
    ax.YAxis.Exponent = 0;
    title(ax, {sprintf('A3. %s, each session, %s', nom_mesure(mes0), texte_t), ligne2}, 'FontSize', 10);
end
box(ax, 'off');

%% B. différence selon la distance
ax = nexttile(tl); hold(ax, 'on');
Tb = Ani(strcmp(Ani.mesure, mes0) & strcmp(Ani.permutation, perm0) & strcmp(Ani.couleur, 'toutes'), :);
animaux = unique(Tb.animal);
if ~isempty(i0)
    patch(ax, i0 + [-0.4 0.4 0.4 -0.4], [-1 -1 1 1] * 1e3, [0.92 0.92 0.92], 'EdgeColor', 'none');
end
for a = 1:numel(animaux)
    m = strcmp(Tb.animal, animaux{a});
    [~, ix] = ismember([Tb.dist_min_um(m), Tb.dist_max_um(m)], dist, 'rows');
    [ix, o] = sort(ix);
    v = Tb.delta(m);
    plot(ax, ix, v(o), '-', 'Color', GRIS, 'LineWidth', 0.8);
end
errorbar(ax, 1:nd, Ad.delta_moyen, Ad.sem, '-o', 'Color', 'k', 'LineWidth', 2, 'MarkerFaceColor', 'k', ...
    'MarkerSize', 5, 'CapSize', 4);
yline(ax, 0, ':k');
yl = [min([Tb.delta; Ad.delta_moyen - Ad.sem]), max([Tb.delta; Ad.delta_moyen + Ad.sem])];
if any(~isfinite(yl)) || yl(1) == yl(2), yl = [-1 1] * 1e-3; end
h_y = yl(2) - yl(1);
ylim(ax, [yl(1) - 0.05 * h_y, yl(2) + 0.15 * h_y]);
for i = 1:nd
    text(ax, i, yl(2) + 0.13 * h_y, texte_p(Ad.p_wilcoxon(i)), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'FontSize', 8, 'Color', couleur_p(Ad.p_wilcoxon(i)));
end
set(ax, 'XTick', 1:nd, 'XTickLabel', etiq, 'XLim', [0.5 nd + 0.5], 'FontSize', 10, 'TickDir', 'out');
xlabel(ax, xlab);
ylabel(ax, sprintf('\\Delta %s (same - different)', nom_mesure(mes0)));
if est_principal
    texte_bande = 'grey band = main test';
else
    texte_bande = 'grey band = radius shown';
end
titre_b = {'B. Difference vs distance'};
if isfield(S, 'MD') && S.MD.ok && strcmp(mes0, opts.mesures{1})
    titre_b{end+1} = sprintf('all distances: %+.1f%%, p = %s', S.MD.pct, texte_p(S.MD.p));
end
ax.YAxis.Exponent = 0;
title(ax, titre_b, 'FontSize', 10);
box(ax, 'off');

%% C. chaque couleur au rayon principal
ax = nexttile(tl); hold(ax, 'on');
rng(1);
vals = [];
for c = 1:numel(noms)
    Tc = sel_ani(mes0, perm0, noms{c}, pd);
    if isempty(Tc), continue; end
    x = c + 0.12 * (rand(height(Tc), 1) - 0.5) * 2;
    scatter(ax, x, Tc.delta, 40, cols(c, :), 'filled', 'MarkerFaceAlpha', 0.85);
    plot(ax, c + [-0.25 0.25], mean(Tc.delta, 'omitnan') * [1 1], 'k-', 'LineWidth', 2);
    vals = [vals; Tc.delta]; %#ok<AGROW>
end
yline(ax, 0, ':k');
if isempty(vals), vals = [-1e-3; 1e-3]; end
yl = [min(vals), max(vals)];
if yl(1) == yl(2), yl = yl + [-1 1] * 1e-3; end
h_y = yl(2) - yl(1);
ylim(ax, [yl(1) - 0.05 * h_y, yl(2) + 0.15 * h_y]);
for c = 1:numel(noms)
    rc = A(strcmp(A.mesure, mes0) & strcmp(A.permutation, perm0) & strcmp(A.couleur, noms{c}) & ...
        A.dist_min_um == pd(1) & A.dist_max_um == pd(2), :);
    if isempty(rc), continue; end
    text(ax, c, yl(2) + 0.13 * h_y, texte_p(rc.p_wilcoxon(1)), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'FontSize', 8, 'Color', couleur_p(rc.p_wilcoxon(1)));
end
set(ax, 'XTick', 1:numel(noms), 'XTickLabel', noms, 'XLim', [0.4 numel(noms) + 0.6], 'FontSize', 10, ...
    'TickDir', 'out');
ylabel(ax, sprintf('\\Delta %s', nom_mesure(mes0)));
ax.YAxis.Exponent = 0;
title(ax, sprintf('C. Each colour, %s', texte_t), 'FontSize', 10);
box(ax, 'off');

%% D. regroupement spatial des paires de même couleur
ax = nexttile(tl); hold(ax, 'on');
Ag = A(strcmp(A.mesure, mes0) & strcmp(A.permutation, 'globale'), :);
if ~isempty(Ag)
    h = gobjects(0); t = {};
    for c = 1:numel(noms)
        Agc = sortrows(Ag(strcmp(Ag.couleur, noms{c}), :), {'dist_min_um', 'dist_max_um'});
        if isempty(Agc), continue; end
        [~, ix] = ismember([Agc.dist_min_um, Agc.dist_max_um], dist, 'rows');
        h(end+1) = plot(ax, ix, Agc.enrichissement, '-o', 'Color', cols(c, :), 'LineWidth', 1, ...
            'MarkerFaceColor', cols(c, :), 'MarkerSize', 3); %#ok<AGROW>
        t{end+1} = noms{c}; %#ok<AGROW>
    end
    Agt = sortrows(Ag(strcmp(Ag.couleur, 'toutes'), :), {'dist_min_um', 'dist_max_um'});
    [~, ix] = ismember([Agt.dist_min_um, Agt.dist_max_um], dist, 'rows');
    h(end+1) = plot(ax, ix, Agt.enrichissement, '-o', 'Color', 'k', 'LineWidth', 2.5, ...
        'MarkerFaceColor', 'k', 'MarkerSize', 5);
    t{end+1} = 'all colours';
    yline(ax, 1, ':k');
    legend(ax, h, t, 'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
    set(ax, 'XTick', 1:nd, 'XTickLabel', etiq, 'XLim', [0.5 nd + 0.5], 'FontSize', 10, 'TickDir', 'out');
    xlabel(ax, xlab);
    ylabel(ax, 'observed / expected');
    title(ax, 'D. Same-colour pairs: observed / expected', 'FontSize', 10);
else
    text(ax, 0.5, 0.5, 'No ''globale'' permutation in config.paires.permutations', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center');
    axis(ax, 'off');
end
box(ax, 'off');

%% E. sensibilité au rayon principal : une ligne par mesure (5 nov. : la
% différence et le Wilcoxon sur les animaux ne dépendent pas de la permutation)
ax = nexttile(tl); hold(ax, 'on');
[mesures, dts] = liste_mesures_paires(opts);
garde = ~strcmp(mesures, 'sce');            % co-participation aux SCE : sa propre figure (pas de %)
mesures = mesures(garde); dts = dts(garde);
% ordre : STTC par fenêtre croissante, puis Pearson du dF/F
cle_ordre = dts; cle_ordre(isnan(cle_ordre)) = Inf;
[~, ordre_m] = sort(cle_ordre);
mesures = mesures(ordre_m); dts = dts(ordre_m);
lab = {}; v = []; e = []; vp = []; ep = []; pw = []; mes_lab = {};
for im = 1:numel(mesures)
    r = A(strcmp(A.mesure, mesures{im}) & strcmp(A.permutation, perm0) & ...
        strcmp(A.couleur, 'toutes') & A.dist_min_um == pd(1) & A.dist_max_um == pd(2), :);
    if isempty(r), continue; end
    if isnan(dts(im))
        lab{end+1} = nom_mesure(mesures{im}); %#ok<AGROW>
    else
        lab{end+1} = sprintf('STTC %g ms', round(1000 * dts(im))); %#ok<AGROW>
    end
    mes_lab{end+1} = mesures{im}; %#ok<AGROW>
    v(end+1) = r.delta_moyen(1); e(end+1) = r.sem(1); %#ok<AGROW>
    vp(end+1) = r.delta_pct(1); ep(end+1) = 100 * r.sem(1) / r.similarite_moyenne(1); %#ok<AGROW>
    pw(end+1) = r.p_wilcoxon(1); %#ok<AGROW>
end
if ~isempty(lab)
    if all(isfinite(vp))
        v = vp; e = ep; unite = '%';
    else
        unite = 'raw difference';
    end
    nl = numel(lab);
    yy = nl:-1:1;
    errorbar(ax, v, yy, e, 'horizontal', 'o', 'Color', [0.45 0.45 0.45], 'MarkerFaceColor', [0.45 0.45 0.45], ...
        'LineWidth', 1.2, 'CapSize', 4);
    k0 = find(strcmp(mes_lab, mes0), 1);
    if isempty(k0), k0 = 1; end
    errorbar(ax, v(k0), yy(k0), e(k0), 'horizontal', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', ...
        'LineWidth', 2, 'CapSize', 5);
    xline(ax, 0, ':k');
    xl = xlim(ax);
    w = xl(2) - xl(1);
    xlim(ax, [xl(1), xl(2) + 0.35 * w]);
    for k = 1:nl
        text(ax, xl(2) + 0.03 * w, yy(k), sprintf('p = %s', texte_p(pw(k))), 'FontSize', 8, ...
            'Color', couleur_p(pw(k)), 'VerticalAlignment', 'middle');
    end
    set(ax, 'YTick', 1:nl, 'YTickLabel', fliplr(lab), 'YLim', [0.4 nl + 0.6], 'FontSize', 9, ...
        'TickLabelInterpreter', 'none', 'TickDir', 'out');
    xlabel(ax, ['\Delta same - different (' unite ')']);
    title(ax, sprintf('E. Sensitivity, %s', texte_t), 'FontSize', 10);
end
box(ax, 'off');

if est_principal
    nom_png = 'paires_agregation.png';
else
    nom_png = sprintf('paires_%s_%gum.png', mes0, pd(2));
end
exportgraphics(fig, fullfile(dossier, nom_png), 'Resolution', 200);
fprintf('%s, %s, %s : différence %+.4f (%s) ; Wilcoxon animaux p = %s ; sessions combinées p = %s -> %s\n', ...
    nom_mesure(mes0), texte_d, perm0, r0.delta_moyen, texte_pct(r0.delta_pct), texte_p(r0.p_wilcoxon), ...
    texte_p(r0.p_combine), nom_png);
end


function t = texte_pct(x)
% pourcentage, ou rien pour une mesure sans échelle positive (phi des SCE)
if isfinite(x), t = sprintf('%+.1f %%', x); else, t = 'no % for this measure'; end
end


function t = nom_mesure(m)
switch m
    case 'sttc', t = 'STTC';
    case 'dff',  t = 'Pearson dF/F';
    case 'sce',  t = 'SCE co-participation';
    otherwise
        t = strrep(strrep(m, 'sttc_', 'STTC '), 'ms', ' ms');
end
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


function c = couleur_p(p)
if p < 0.05
    c = [0.75 0.1 0.1];
else
    c = [0.45 0.45 0.45];
end
end
