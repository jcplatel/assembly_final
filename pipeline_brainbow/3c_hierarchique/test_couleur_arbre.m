function R = test_couleur_arbre(M, col, groupes, opts, graine)
%TEST_COULEUR_ARBRE  Couleur de chaque groupe d'un arbre hiérarchique,
%   contre des couleurs mélangées.
%   M       : noeuds x cellules de l'arbre (noeuds_arbre)
%   col     : couleur de chaque cellule de l'arbre (1 à 6 = colorée ; autre
%             = non colorée, ignorée dans le test)
%   groupes : groupe de mélange de chaque cellule (tous à 1 = tout le champ ;
%             blocs spatiaux pour le contrôle)
%
%   Pour chaque noeud, avec m cellules colorées dont k_c de la couleur c :
%   p_c = P(X >= k_c), X hypergéométrique (n cellules colorées dans l'arbre,
%   n_c de la couleur c, m tirées) ; T = min des 6 p_c ; couleur = celle du
%   minimum. Le noeud « a une couleur en excès » si T < opts.alpha.
%   Comme pour les assemblées, ce p n'est pas corrigé pour les 6 couleurs ni
%   pour les nombreux noeuds (emboîtés) : c'est le hasard qui dit combien de
%   noeuds passeraient sans aucun lien couleur / activité.
%   Par classe de taille (opts.classes, en cellules colorées) :
%     obs     : nombre de noeuds avec une couleur en excès
%     nul     : la même chose, couleurs mélangées (opts.nPerm fois, à
%               l'intérieur des groupes), MÊME arbre
%     attendu : moyenne de nul ; p = part des mélanges >= obs (unilatéral)

colorees = col(:) >= 1 & col(:) <= 6;
idx = find(colorees);
n = numel(idx);
cc = round(col(idx));
cc = cc(:);
nCl = numel(opts.classes) - 1;
Mc = double(M(:, idx));
m = sum(Mc, 2);
classe = discretize(m, opts.classes);          % NaN sous la plus petite taille

% tables des p hypergéométriques (les mélanges gardent n et chaque n_c)
nc = accumarray(cc, 1, [6 1])';
[mm, kk] = ndgrid(0:n, 0:n);
tab = cell(1, 6);
for c = 1:6
    P = hygecdf(kk - 1, n, nc(c), mm, 'upper');   % P(X > k - 1) = P(X >= k)
    P(kk == 0) = 1;
    P(kk > mm) = NaN;
    tab{c} = P;
end

[R.T, R.couleur, R.k] = statistique(Mc, cc, m, tab);
R.m = m;
R.classe = classe;
R.n = n;
R.n_par_couleur = nc;
R.n_noeuds = accumarray(classe(~isnan(classe)), 1, [nCl 1])';
R.obs = compter(R.T, classe, opts.alpha, nCl);

% hasard : couleurs mélangées entre cellules colorées du même groupe
fixer_flux_aleatoire(graine, 7, 1);
g = groupes(idx);
g = g(:);
[~, ~, g] = unique(g);
membres = accumarray(g, (1:n)', [], @(v) {v});
R.nul = zeros(opts.nPerm, nCl);
for p = 1:opts.nPerm
    cp = cc;
    for u = 1:numel(membres)
        j = membres{u};
        if numel(j) > 1
            cp(j) = cc(j(randperm(numel(j))));
        end
    end
    R.nul(p, :) = compter(statistique(Mc, cp, m, tab), classe, opts.alpha, nCl);
end
R.attendu = mean(R.nul, 1);
R.p = (1 + sum(R.nul >= R.obs, 1)) / (opts.nPerm + 1);
end


function [T, couleur, K] = statistique(Mc, cc, m, tab)
X = zeros(numel(cc), 6);
X(sub2ind(size(X), (1:numel(cc))', cc)) = 1;
K = Mc * X;                                    % noeuds x couleurs : k_c
Pn = ones(size(K));
for c = 1:6
    Pn(:, c) = tab{c}(sub2ind(size(tab{c}), m + 1, K(:, c) + 1));
end
[T, couleur] = min(Pn, [], 2);
T(m < 1) = NaN;
end


function c = compter(T, classe, alpha, nCl)
ok = ~isnan(classe) & T < alpha;
c = accumarray(classe(ok), 1, [nCl 1])';
end
