function [Ts, Ta, Ra] = meta_analyse_par_animal(dossier_analyse, fichier_sessions)
%META_ANALYSE_PAR_ANIMAL  Bilan animal par animal d'une analyse (18 nov.,
%   demandée par l'utilisateur) : activité, couleurs, SCE, assemblées, et
%   distribution des p_minp selon la couleur en excès et la taille.
%   Aucun nouveau test sur les assemblées : relit l'étape 1 (cache), le
%   choix du K, les couleurs (colonne fichier_colorcell, sinon GT) et la
%   table R de l'étape 4 (champ_globale).
%
%   [Ts, Ta, Ra] = meta_analyse_par_animal( ...
%       'E:\...\analyses\selection34_MaxwellAuto_v2', 'E:\...\sessions_brainbow.xlsx')
%   Ts : une ligne par session ; Ta : une ligne par animal ; Ra : assemblées
%   Écrit dans <dossier_analyse>\agregation\meta_animaux\ :
%     meta_par_animal.png, pminp_couleur_taille.png, meta_sessions.csv,
%     meta_animaux.csv, meta_assemblees.csv

dossier_analyse = regexprep(char(dossier_analyse), '[\\/]+$', '');
racine1 = fileparts(fileparts(dossier_analyse));          % dossier de l'étape 1
sortie = fullfile(dossier_analyse, 'agregation', 'meta_animaux');
if ~isfolder(sortie), mkdir(sortie); end
noms_coul = {'R', 'G', 'B', 'RG', 'RB', 'GB'};
cols = [0.86 0.16 0.16; 0.16 0.66 0.22; 0.20 0.38 0.92; 0.93 0.78 0.10; 0.82 0.22 0.80; 0.10 0.74 0.80];
REF_MED = [34.0 5.5 5.0 19.5 27.5 10.0];        % médianes de référence (violon de l'utilisateur)

%% ---------- sessions ----------
T = readtable(fichier_sessions);
T = T(T.inclure == 1, :);
B = load(fullfile(dossier_analyse, 'best_k_par_session.mat'), 'best_k'); bk = B.best_k;
n = height(T);
Ts = table();
Ts.animal = T.id_animal; Ts.session = strcat(T.id_session, '_', T.id_champ);
cle = strcat(T.id_animal, '_', T.id_session, '_', T.id_champ, '_', T.id_condition);
z = nan(n, 1);
[Ts.cellules, Ts.colorees, Ts.transients, Ts.repos_min, Ts.seuil_SCE, Ts.SCE, Ts.SCE_hasard, ...
    Ts.K, Ts.assemblees, Ts.SCE_clusters] = deal(z, z, z, z, z, z, z, z, z, z);
Ts.prop = nan(n, 6);
for i = 1:n
    P = load(fullfile(racine1, cle{i}, 'pretraitement.mat'), 'NCell', 'TRace', 'n_transients_total', 'diagnostic');
    D = P.diagnostic;
    Ts.cellules(i) = P.NCell; Ts.transients(i) = P.n_transients_total; Ts.SCE(i) = numel(P.TRace);
    Ts.repos_min(i) = D.duree_repos_s / 60; Ts.seuil_SCE(i) = D.seuil_SCE; Ts.SCE_hasard(i) = D.n_SCE_hasard;
    fc = T.fichier_colorcell{i};
    if isempty(fc) || ~ischar(fc), fc = fullfile(T.dossier_suite2p{i}, 'colorcell_GT.mat'); end
    if isfile(fc)
        C = load(fc); fn = fieldnames(C); c = C.(fn{1}); c = c(:);
        col = c(c >= 1 & c <= 6);
        Ts.colorees(i) = numel(col);
        Ts.prop(i, :) = 100 * histcounts(col, 0.5:1:6.5) / max(numel(col), 1);
    end
    j = find(strcmp(bk.session, cle{i}), 1);
    if ~isempty(j) && strcmp(bk.statut{j}, 'ok') && ~isnan(bk.K(j))
        Ts.K(i) = bk.K(j);
        X = load(fullfile(racine1, cle{i}, sprintf('k%d', bk.K(j)), 'assemblees_v2.mat'), 'NClOK', 'NRaceOK');
        Ts.assemblees(i) = X.NClOK; Ts.SCE_clusters(i) = X.NRaceOK;
    else
        Ts.assemblees(i) = 0; Ts.SCE_clusters(i) = 0;
    end
end
Ts.pct_colorees = 100 * Ts.colorees ./ Ts.cellules;
Ts.transients_cell_min = Ts.transients ./ Ts.cellules ./ Ts.repos_min;
Ts.SCE_min = Ts.SCE ./ Ts.repos_min;
Ts.pct_SCE_clusters = 100 * Ts.SCE_clusters ./ Ts.SCE;

%% ---------- assemblées ----------
S = load(fullfile(dossier_analyse, 'agregation', 'champ_globale', 'resultats_biais_clonal.mat'), 'R');
R = S.R;
Ra = R(:, {'animal', 'session', 'champ', 'K', 'assemblee', 'n_total', 'n_colorees', 'part_colorees', ...
    'n_SCE_signif', 'compacite_pct', 'ratio_etalement', 'frac_active_moy', 'p_minp', 'taux_nul', ...
    'couleur_dominante', 'part_couleur_dominante', 'testable', 'exclue_criteres', 'q'});
[~, Ra.ic] = ismember(Ra.couleur_dominante, noms_coul);
% part de la couleur en excès dans le champ (référence de l'urne)
pc = nan(height(R), 1);
for i = 1:height(R)
    if Ra.ic(i) > 0 && sum(R.Kcol_urne(i, :)) > 0, pc(i) = R.Kcol_urne(i, Ra.ic(i)) / sum(R.Kcol_urne(i, :)); end
end
Ra.part_champ = pc;
Ra.enrichissement = Ra.part_couleur_dominante ./ Ra.part_champ;

%% ---------- animaux ----------
an = unique(Ts.animal); nA = numel(an);
Ta = table(an, 'VariableNames', {'animal'});
f_sum = {'cellules', 'colorees', 'SCE', 'assemblees'};
for k = 1:numel(f_sum), Ta.(f_sum{k}) = zeros(nA, 1); end
[Ta.sessions, Ta.pct_colorees, Ta.transients_cell_min, Ta.SCE_min, Ta.seuil_SCE, Ta.pct_SCE_clusters, ...
    Ta.testables, Ta.sig05, Ta.sig20, Ta.FDR, Ta.taille_med, Ta.SCE_signif_med, Ta.compacite_med, Ta.p_med] = ...
    deal(zeros(nA, 1), nan(nA, 1), nan(nA, 1), nan(nA, 1), nan(nA, 1), nan(nA, 1), zeros(nA, 1), zeros(nA, 1), ...
    zeros(nA, 1), zeros(nA, 1), nan(nA, 1), nan(nA, 1), nan(nA, 1), nan(nA, 1));
Ta.prop = nan(nA, 6);
for a = 1:nA
    m = strcmp(Ts.animal, an{a});
    Ta.sessions(a) = sum(m);
    for k = 1:numel(f_sum), Ta.(f_sum{k})(a) = sum(Ts.(f_sum{k})(m), 'omitnan'); end
    Ta.pct_colorees(a) = 100 * Ta.colorees(a) / Ta.cellules(a);
    Ta.transients_cell_min(a) = median(Ts.transients_cell_min(m));
    Ta.SCE_min(a) = median(Ts.SCE_min(m));
    Ta.seuil_SCE(a) = median(Ts.seuil_SCE(m));
    Ta.pct_SCE_clusters(a) = 100 * sum(Ts.SCE_clusters(m)) / sum(Ts.SCE(m));
    Ta.prop(a, :) = sum(Ts.prop(m, :) .* Ts.colorees(m), 1, 'omitnan') / sum(Ts.colorees(m), 'omitnan');
    r = strcmp(Ra.animal, an{a});
    t = r & Ra.testable;
    Ta.testables(a) = sum(t); Ta.sig05(a) = sum(t & Ra.p_minp < 0.05); Ta.sig20(a) = sum(t & Ra.p_minp < 0.2);
    Ta.FDR(a) = sum(t & Ra.q < 0.1);
    Ta.taille_med(a) = median(Ra.n_total(r)); Ta.SCE_signif_med(a) = median(Ra.n_SCE_signif(r));
    Ta.compacite_med(a) = median(Ra.compacite_pct(r)); Ta.p_med(a) = median(Ra.p_minp(t));
end

%% ---------- figure 1 : par animal ----------
f1 = figure('Color', 'w', 'Position', [30 30 1800 1000]);
tl = tiledlayout(f1, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Animal by animal: %d sessions, %d animals (%s)', n, nA, strrep(dossier_analyse, '\', '/')), ...
    'FontWeight', 'bold', 'Interpreter', 'none');
[~, ia] = ismember(Ts.animal, an);
x = ia + 0.12 * randn(n, 1);
pan = {'cellules', 'cells per session'; 'pct_colorees', '% coloured cells'; 'transients_cell_min', ...
    'transients / cell / min of rest'; 'SCE_min', 'SCE / min of rest'; 'pct_SCE_clusters', '% SCE in retained clusters'};
for k = 1:size(pan, 1)
    ax = nexttile(tl); hold(ax, 'on');
    v = Ts.(pan{k, 1});
    scatter(ax, x, v, 30, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.7);
    for a = 1:nA, plot(ax, a + [-0.3 0.3], median(v(ia == a), 'omitnan') * [1 1], 'k-', 'LineWidth', 2); end
    habiller(ax, an); title(ax, pan{k, 2}, 'FontWeight', 'normal');
end
% couleurs par animal
ax = nexttile(tl); hold(ax, 'on');
hb = bar(ax, 1:nA + 1, [Ta.prop; REF_MED], 'stacked', 'EdgeColor', 'none');
for c = 1:6, hb(c).FaceColor = cols(c, :); end
set(ax, 'XTick', 1:nA + 1, 'XTickLabel', [an; {'ref.'}], 'XTickLabelRotation', 45, 'TickDir', 'out');
ylim(ax, [0 100]); title(ax, 'colour proportions (coloured cells)', 'FontWeight', 'normal');
% assemblées : taille et SCE significatifs
ax = nexttile(tl); hold(ax, 'on');
[~, ir] = ismember(Ra.animal, an); xr = ir + 0.12 * randn(height(Ra), 1);
scatter(ax, xr(~Ra.testable), Ra.n_total(~Ra.testable), 25, [0.6 0.6 0.6]);
scatter(ax, xr(Ra.testable), Ra.n_total(Ra.testable), 30, [0.2 0.4 0.8], 'filled');
habiller(ax, an); title(ax, 'assembly size (cells); filled = testable', 'FontWeight', 'normal');
ax = nexttile(tl); hold(ax, 'on');
for i = find(Ra.testable)'
    c = [0.5 0.5 0.5]; if Ra.ic(i) > 0, c = cols(Ra.ic(i), :); end
    scatter(ax, xr(i), Ra.p_minp(i), 40, c, 'filled', 'MarkerEdgeColor', 'k');
end
set(ax, 'YScale', 'log'); yline(ax, [0.05 0.2], '--', {'0.05', '0.2'});
habiller(ax, an); title(ax, 'p_{minP} of testable assemblies (colour in excess)', 'FontWeight', 'normal');
exportgraphics(f1, fullfile(sortie, 'meta_par_animal.png'), 'Resolution', 150);
close(f1);

%% ---------- figure 2 : p_minp par couleur et par taille ----------
f2 = figure('Color', 'w', 'Position', [30 30 1800 950]);
tl = tiledlayout(f2, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
t = Ra.testable; nt = sum(t);
title(tl, sprintf('p_{minP} by colour in excess and by size (%d testable assemblies; grey: %d excluded by size/SCE criteria)', ...
    nt, sum(~t)), 'FontWeight', 'bold');
% a. histogramme empilé par couleur
ax = nexttile(tl); hold(ax, 'on');
bords = 0:0.05:1;
H = zeros(numel(bords) - 1, 6);
for c = 1:6, H(:, c) = histcounts(Ra.p_minp(t & Ra.ic == c), bords)'; end
hb = bar(ax, bords(1:end-1) + 0.025, H, 'stacked', 'BarWidth', 1, 'EdgeColor', 'w');
for c = 1:6, hb(c).FaceColor = cols(c, :); end
yline(ax, nt * 0.05, 'k--', 'uniform (no effect)');
xlabel(ax, 'p_{minP}'); ylabel(ax, 'testable assemblies'); legend(ax, hb, noms_coul, 'Box', 'off', 'NumColumns', 2);
title(ax, 'a. Distribution (bins of 0.05)', 'FontWeight', 'normal');
% b. cumul observé contre attendu, par seuil
ax = nexttile(tl); hold(ax, 'on');
seuils = [0.01 0.05 0.1 0.2 0.3 0.5];
Cc = zeros(numel(seuils), 6);
for k = 1:numel(seuils), for c = 1:6, Cc(k, c) = sum(t & Ra.ic == c & Ra.p_minp <= seuils(k)); end, end
hb = bar(ax, 1:numel(seuils), Cc, 'stacked', 'EdgeColor', 'none');
for c = 1:6, hb(c).FaceColor = cols(c, :); end
plot(ax, 1:numel(seuils), nt * seuils, 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'w');
for k = 1:numel(seuils)
    text(ax, k, sum(Cc(k, :)) + 0.4, sprintf('%d / %.1f', sum(Cc(k, :)), nt * seuils(k)), ...
        'HorizontalAlignment', 'center', 'FontSize', 8);
end
set(ax, 'XTick', 1:numel(seuils), 'XTickLabel', compose('<= %g', seuils));
ylabel(ax, 'testable assemblies'); title(ax, 'b. Observed (bars) vs expected without effect (line, \alpha x N)', 'FontWeight', 'normal');
% c. part des couleurs en excès à p < 0,2 contre le champ
ax = nexttile(tl); hold(ax, 'on');
m20 = t & Ra.p_minp < 0.2;
obs = arrayfun(@(c) sum(m20 & Ra.ic == c), 1:6); obs = 100 * obs / max(sum(obs), 1);
tous = arrayfun(@(c) sum(t & Ra.ic == c), 1:6); tous = 100 * tous / max(sum(tous), 1);
hb = bar(ax, 1:6, [obs; tous; REF_MED]', 'grouped', 'EdgeColor', 'none');
hb(1).FaceColor = 'flat'; hb(1).CData = cols; hb(2).FaceColor = [0.6 0.6 0.6]; hb(3).FaceColor = [0.85 0.85 0.85];
set(ax, 'XTick', 1:6, 'XTickLabel', noms_coul); ylabel(ax, '%');
legend(ax, {sprintf('colour in excess, p < 0.2 (%d)', sum(m20)), sprintf('colour in excess, all testable (%d)', nt), ...
    'field reference (median)'}, 'Box', 'off', 'Location', 'northeast');
title(ax, 'c. Which colour wins?', 'FontWeight', 'normal');
% d. p_minp contre nombre de cellules colorées
ax = nexttile(tl); hold(ax, 'on');
scatter(ax, Ra.n_colorees(~t), Ra.p_minp(~t), 25, [0.7 0.7 0.7]);
for c = 1:6
    m = t & Ra.ic == c;
    scatter(ax, Ra.n_colorees(m), Ra.p_minp(m), 45, cols(c, :), 'filled', 'MarkerEdgeColor', 'k');
end
set(ax, 'YScale', 'log'); yline(ax, [0.05 0.2], '--');
xlabel(ax, 'coloured cells in the assembly'); ylabel(ax, 'p_{minP}');
[rho, pr] = corr(Ra.n_colorees(t), Ra.p_minp(t), 'type', 'Spearman');
title(ax, sprintf('d. p_{minP} vs size (testable: Spearman \\rho = %.2f, p = %.2f)', rho, pr), 'FontWeight', 'normal');
% e. p_minp contre enrichissement
ax = nexttile(tl); hold(ax, 'on');
for c = 1:6
    m = t & Ra.ic == c;
    scatter(ax, Ra.enrichissement(m), Ra.p_minp(m), 45, cols(c, :), 'filled', 'MarkerEdgeColor', 'k');
end
set(ax, 'YScale', 'log', 'XScale', 'log'); yline(ax, [0.05 0.2], '--'); xline(ax, 1, 'k:');
xlabel(ax, 'enrichment of the colour in excess (assembly share / field share)'); ylabel(ax, 'p_{minP}');
title(ax, 'e. p_{minP} vs strength of the excess', 'FontWeight', 'normal');
% f. par animal, part à p < 0,2
ax = nexttile(tl); hold(ax, 'on');
hb = bar(ax, 1:nA, [Ta.sig05, Ta.sig20 - Ta.sig05, Ta.testables - Ta.sig20], 'stacked', 'EdgeColor', 'none');
hb(1).FaceColor = [0.75 0.15 0.15]; hb(2).FaceColor = [0.95 0.6 0.5]; hb(3).FaceColor = [0.85 0.85 0.85];
plot(ax, 1:nA, 0.2 * Ta.testables, 'k_', 'MarkerSize', 18, 'LineWidth', 2);
habiller(ax, an); ylabel(ax, 'testable assemblies');
legend(ax, {'p < 0.05', '0.05 - 0.2', 'p >= 0.2', 'expected at p < 0.2'}, 'Box', 'off');
title(ax, 'f. Per animal', 'FontWeight', 'normal');
exportgraphics(f2, fullfile(sortie, 'pminp_couleur_taille.png'), 'Resolution', 150);
close(f2);

writetable(splitvars(Ts), fullfile(sortie, 'meta_sessions.csv'));
writetable(splitvars(Ta), fullfile(sortie, 'meta_animaux.csv'));
writetable(Ra, fullfile(sortie, 'meta_assemblees.csv'));
end


function habiller(ax, an)
set(ax, 'XTick', 1:numel(an), 'XTickLabel', an, 'XLim', [0.4 numel(an) + 0.6], 'XTickLabelRotation', 45, ...
    'TickDir', 'out', 'TickLabelInterpreter', 'none');
box(ax, 'off'); grid(ax, 'on');
end
