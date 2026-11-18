function T = projection_assemblees(cfg, dossier_analyse, dossier_pc, dossier_sortie, n_perm)
%% =========================================================================
%  PROJECTION_ASSEMBLEES (18 nov. 2026, demandée par l'utilisateur)
%  Les assemblées (réactivations au repos) se projettent-elles sur
%  l'activité de course : place cells (sessions cues) et séquences de
%  course (étape 5) ? AUCUNE couleur dans les tests 1 à 3.
%
%  Assemblées : K retenu à l'étape 2 (best_k_par_session.mat), assemblees_v2.mat,
%  toutes celles du clustering (>= 6 cellules) ; une cellule peut appartenir
%  à plusieurs assemblées (v2).
%   1. Les membres des assemblées sont-ils plus souvent place cells que les
%      autres cellules ? (p par permutation des étiquettes « membre »)
%   2. Les place cells d'une même assemblée ont-elles des champs plus proches
%      que des place cells tirées au hasard ? (distance circulaire moyenne des
%      paires de même assemblée ; permutation des assemblées entre place
%      cells membres) -> une assemblée = un segment du tapis (Malvache 2016,
%      fig. S5)
%   3. Idem avec le décalage dans la séquence de course (étape 5,
%      decalage_tout, cellules de la séquence)
%   4. Couleur : pour chaque assemblée, distance des champs de ses place cells
%      colorées de la couleur dominante de l'assemblée contre les autres
%      (descriptif, à relier au biais clonal)
%  Écrit projection_assemblees_par_session.csv, projection_assemblees.png.
%  =========================================================================
if nargin < 5, n_perm = 5000; end
if ~isfolder(dossier_sortie), mkdir(dossier_sortie); end
cfg = dossiers_analyse(cfg);
S = lire_table_sessions(cfg);
B = load(fullfile(dossier_analyse, 'best_k_par_session.mat'), 'best_k'); B = B.best_k;
lignes = {};
for s = 1:height(S)
    nom = S.nom{s};
    r = struct('session', nom, 'animal', S.id_animal{s}, 'statut', '', 'n_asm', NaN, 'n_membres', NaN, ...
        'pct_pc_membres', NaN, 'pct_pc_autres', NaN, 'p_pc_membres', NaN, ...
        'n_paires_pc', NaN, 'dist_champs_meme', NaN, 'dist_champs_hasard', NaN, 'p_champs', NaN, ...
        'n_paires_seq', NaN, 'decal_meme', NaN, 'decal_hasard', NaN, 'p_sequence', NaN);
    try
        i = find(strcmp(B.session, nom), 1);
        if isempty(i) || isnan(B.K(i)), error('pas de K retenu'); end
        A = load(fullfile(S.dossier_source{s}, sprintf('k%d', B.K(i)), 'assemblees_v2.mat'), 'CellCl_final');
        C = A.CellCl_final > 0; C = C(sum(C, 2) >= 6, :);
        r.n_asm = size(C, 1);
        if r.n_asm == 0, error('aucune assemblée'); end
        membre = any(C, 1)'; r.n_membres = sum(membre);
        rng(s);
        % place cells
        fpc = fullfile(dossier_pc, nom, 'place_cells.mat');
        if isfile(fpc)
            P = load(fpc, 'est_pc', 'pc');
            if isfield(P, 'est_pc')
                pcok = P.est_pc(:);
                r.pct_pc_membres = 100 * mean(pcok(membre)); r.pct_pc_autres = 100 * mean(pcok(~membre));
                obs = mean(pcok(membre)); nul = zeros(n_perm, 1);
                for b = 1:n_perm, m = membre(randperm(numel(membre))); nul(b) = mean(pcok(m)); end
                r.p_pc_membres = (1 + sum(nul >= obs)) / (n_perm + 1);
                pos = P.pc.position(:);
                [r.n_paires_pc, r.dist_champs_meme, r.dist_champs_hasard, r.p_champs] = ...
                    test_proximite(C, pcok, pos, @(a, b) min(abs(a - b), 100 - abs(a - b)), n_perm);
            end
        end
        % séquence de course (étape 5)
        fs = fullfile(dossier_analyse, nom, 'sequences', 'sequences.mat');
        if isfile(fs)
            Q = load(fs);
            if isfield(Q, 'sequence_tout')
                sq = Q.sequence_tout(:); dec = double(Q.decalage_tout(:));
                [r.n_paires_seq, r.decal_meme, r.decal_hasard, r.p_sequence] = ...
                    test_proximite(C, sq, dec, @(a, b) abs(a - b), n_perm);
            end
        end
        r.statut = 'ok';
    catch err
        r.statut = err.message;
    end
    lignes{end+1} = r; %#ok<AGROW>
    fprintf('%d/%d %s : %s\n', s, height(S), nom, r.statut);
end
T = struct2table([lignes{:}]);
writetable(T, fullfile(dossier_sortie, 'projection_assemblees_par_session.csv'));

fig = figure('Color', 'w', 'Position', [50 50 1400 420], 'Visible', 'off');
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact');
ax = nexttile(tl); hold(ax, 'on');
ok = ~isnan(T.pct_pc_membres);
plot(ax, [1 2], [T.pct_pc_autres(ok) T.pct_pc_membres(ok)]', '-', 'Color', [.7 .7 .7]);
scatter(ax, ones(sum(ok), 1), T.pct_pc_autres(ok), 25, 'k', 'filled'); scatter(ax, 2 * ones(sum(ok), 1), T.pct_pc_membres(ok), 25, 'r', 'filled');
set(ax, 'XTick', [1 2], 'XTickLabel', {'hors assemblée', 'membres'}, 'XLim', [0.5 2.5]); ylabel(ax, '% de place cells');
title(ax, sprintf('1. Place cells parmi les membres (%d sessions, %d à p < 0,05)', sum(ok), sum(T.p_pc_membres < 0.05)));
for k = 1:2
    ax = nexttile(tl);
    if k == 1, p = T.p_champs; t = '2. Champs des place cells d''une même assemblée plus proches ?';
    else, p = T.p_sequence; t = '3. Membres d''une même assemblée voisins dans la séquence de course ?'; end
    p = p(~isnan(p)); histogram(ax, p, 0:0.1:1); yline(ax, numel(p) / 10, 'r--');
    if ~isempty(p), pf = 1 - chi2cdf(-2 * sum(log(p)), 2 * numel(p)); else, pf = NaN; end
    xlabel(ax, 'p par session'); title(ax, {t, sprintf('%d sessions, %d à p < 0,05 ; Fisher p = %.3g', numel(p), sum(p < 0.05), pf)});
end
exportgraphics(fig, fullfile(dossier_sortie, 'projection_assemblees.png'), 'Resolution', 150);
close(fig);
end


function [n_paires, d_meme, d_hasard, p] = test_proximite(C, garde, x, dist, n_perm)
% paires de cellules (gardées) d'une même assemblée : distance moyenne de x ;
% hasard : on permute les cellules gardées entre les assemblées (chaque
% assemblée garde son nombre de cellules gardées)
n_paires = 0; d_meme = NaN; d_hasard = NaN; p = NaN;
Cg = C(:, garde(:)'); xg = x(garde);
if numel(xg) < 4, return; end
D = dist(xg(:), xg(:)');
stat = @(M) moyenne_meme(M, D);
n_paires = paires(Cg);
if n_paires < 3, return; end
d_meme = stat(Cg);
nul = zeros(n_perm, 1);
for b = 1:n_perm, nul(b) = stat(Cg(:, randperm(size(Cg, 2)))); end
d_hasard = mean(nul);
p = (1 + sum(nul <= d_meme)) / (n_perm + 1);
end

function m = moyenne_meme(M, D)
s = 0; n = 0;
for a = 1:size(M, 1)
    k = find(M(a, :));
    if numel(k) < 2, continue; end
    sub = D(k, k); h = triu(true(numel(k)), 1);
    s = s + sum(sub(h)); n = n + sum(h(:));
end
m = s / max(n, 1);
end

function n = paires(M)
k = sum(M, 2); n = sum(k .* (k - 1) / 2);
end
