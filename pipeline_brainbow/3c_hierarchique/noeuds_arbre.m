function [M, ordre, x, hauteur] = noeuds_arbre(Z)
%NOEUDS_ARBRE  Les groupes (noeuds) d'un arbre hiérarchique (sortie de linkage).
%   M       : noeuds x cellules (logique) ; M(i, j) = la cellule j est sous
%             le noeud i. Le noeud i réunit Z(i, 1) et Z(i, 2) ; il y a
%             N - 1 noeuds pour N cellules, du plus bas (petits groupes) au
%             plus haut (la racine = toutes les cellules).
%   ordre   : ordre des feuilles de gauche à droite (pour la figure)
%   x       : abscisse de chaque noeud (milieu de ses deux enfants)
%   hauteur : hauteur de chaque noeud (= Z(:, 3), distance de fusion)

N = size(Z, 1) + 1;
% ordre des feuilles : parcours depuis la racine, enfant de gauche d'abord
ordre = zeros(1, N);
pile = 2 * N - 1;
k = 0;
while ~isempty(pile)
    c = pile(end);
    pile(end) = [];
    if c <= N
        k = k + 1;
        ordre(k) = c;
    else
        pile = [pile, Z(c - N, 2), Z(c - N, 1)]; %#ok<AGROW>
    end
end
x_feuille = zeros(1, N);
x_feuille(ordre) = 1:N;

M = false(N - 1, N);
x = zeros(N - 1, 1);
for i = 1:N - 1
    xs = zeros(1, 2);
    for e = 1:2
        c = Z(i, e);
        if c <= N
            M(i, c) = true;
            xs(e) = x_feuille(c);
        else
            M(i, :) = M(i, :) | M(c - N, :);
            xs(e) = x(c - N);
        end
    end
    x(i) = mean(xs);
end
hauteur = Z(:, 3);
end
