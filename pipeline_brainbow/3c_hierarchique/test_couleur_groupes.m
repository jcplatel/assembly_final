function R = test_couleur_groupes(M, col, groupes, opts, graine)
%TEST_COULEUR_GROUPES  Couleur des groupes d'UNE coupe de l'arbre (groupes
%   séparés), contre des couleurs mélangées (7 nov., à la place des classes
%   de taille : demandé par l'utilisateur pour ne pas multiplier les tests).
%   M       : groupes x cellules de l'arbre (logique), une coupe
%   col     : couleur de chaque cellule (1 à 6 = colorée ; autre = ignorée)
%   groupes : groupe de mélange de chaque cellule (1 = tout le champ ; blocs)
%
%   Groupe testé : au moins opts.min_colorees cellules colorées.
%   Par groupe : T = min sur les 6 couleurs du p hypergéométrique (au moins
%   k_c cellules de la couleur c parmi ses m cellules colorées) ; couleur =
%   celle du minimum.
%   1. Test de la coupe (un seul par niveau) : nombre de groupes avec
%      T < opts.alpha, contre la même chose après mélange des couleurs
%      (opts.nPerm fois) ; p unilatéral. Le mélange tient compte des 6
%      couleurs et du nombre de groupes.
%   2. Quels groupes ? p_minP de chaque groupe = part des mélanges où son T
%      est au moins aussi petit (corrige le choix de la meilleure couleur,
%      comme pour les assemblées), puis FDR de Benjamini-Hochberg sur les
%      groupes testés de cette coupe (opts.fdr_q).

colorees = col(:) >= 1 & col(:) <= 6;
idx = find(colorees);
n = numel(idx);
cc = round(col(idx));
cc = cc(:);
Mc = double(M(:, idx));
m = sum(Mc, 2);
teste = m >= opts.min_colorees;
nG = size(M, 1);

nc = accumarray(cc, 1, [6 1])';
[mm, kk] = ndgrid(0:n, 0:n);
tab = cell(1, 6);
for c = 1:6
    P = hygecdf(kk - 1, n, nc(c), mm, 'upper');   % P(X >= k)
    P(kk == 0) = 1;
    P(kk > mm) = NaN;
    tab{c} = P;
end

[T, couleur, K] = statistique(Mc, cc, m, tab);
T(~teste) = NaN;
R.m = m;
R.T = T;
R.couleur = couleur;
R.k = K;
R.teste = teste;
R.n_groupes = sum(teste);
R.obs = sum(T(teste) < opts.alpha);

fixer_flux_aleatoire(graine, 9, 1);
g = groupes(idx);
[~, ~, g] = unique(g(:));
membres = accumarray(g, (1:n)', [], @(v) {v});
it = find(teste);
Tp = zeros(opts.nPerm, numel(it));
for p = 1:opts.nPerm
    cp = cc;
    for u = 1:numel(membres)
        j = membres{u};
        if numel(j) > 1
            cp(j) = cc(j(randperm(numel(j))));
        end
    end
    t = statistique(Mc(it, :), cp, m(it), tab);
    Tp(p, :) = t';
end
R.nul = sum(Tp < opts.alpha, 2);
R.attendu = mean(R.nul);
R.p = (1 + sum(R.nul >= R.obs)) / (opts.nPerm + 1);
R.p_minP = nan(nG, 1);
R.fdr = false(nG, 1);
if ~isempty(it)
    R.p_minP(it) = (1 + sum(Tp <= T(it)' * (1 + 1e-10), 1))' / (opts.nPerm + 1);
    R.fdr(it) = bh_fdr_mask(R.p_minP(it)', opts.fdr_q)';
end
R.n_fdr = sum(R.fdr);
end


function [T, couleur, K] = statistique(Mc, cc, m, tab)
X = zeros(numel(cc), 6);
X(sub2ind(size(X), (1:numel(cc))', cc)) = 1;
K = Mc * X;
Pn = ones(size(K));
for c = 1:6
    Pn(:, c) = tab{c}(sub2ind(size(tab{c}), m + 1, K(:, c) + 1));
end
[T, couleur] = min(Pn, [], 2);
T(m < 1) = NaN;
end
