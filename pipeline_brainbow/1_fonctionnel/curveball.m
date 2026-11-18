function R = curveball(R, n_echanges)
% CURVEBALL  Mélange une matrice binaire en gardant EXACTEMENT les sommes
% de lignes (nombre de SCE de chaque cellule) et de colonnes (nombre de
% cellules de chaque SCE). Strona et al. 2014, Nature Communications.
%
% À chaque échange, deux cellules mettent en commun les SCE que l'une a et
% pas l'autre, puis se les redistribuent au hasard (chacune garde son
% nombre de SCE).
%   R : matrice logique (cellules x SCE)
%   n_echanges : nombre d'échanges (défaut : 5 x nombre de cellules actives)

lignes_actives = find(any(R, 2));
nA = numel(lignes_actives);
if nA < 2, return; end
if nargin < 2, n_echanges = 5 * nA; end

for t = 1:n_echanges
    paire = lignes_actives(randperm(nA, 2));
    a = R(paire(1), :);
    b = R(paire(2), :);
    seul_a = find(a & ~b);
    seul_b = find(b & ~a);
    na = numel(seul_a);
    if na == 0 || isempty(seul_b), continue; end
    pool = [seul_a, seul_b];
    pool = pool(randperm(numel(pool)));
    R(paire(1), pool) = false;
    R(paire(2), pool) = false;
    R(paire(1), pool(1:na))     = true;
    R(paire(2), pool(na+1:end)) = true;
end
end
