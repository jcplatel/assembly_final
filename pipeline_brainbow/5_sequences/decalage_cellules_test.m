function [corr_max, decalage, p, dans_sequence, ordre] = decalage_cellules_test(DF, debuts, courses, ...
    trames_course, max_lag, n_faux, q_fdr)
%DECALAGE_CELLULES_TEST  Place de chaque cellule par rapport au début de
%   course, avec un test par cellule (20 nov. 2026, choix de l'utilisateur).
%   Même idée que decalage_cellules (code de l'utilisateur) : corrélation
%   croisée normalisée ('coeff') entre les débuts de course (1 au début de
%   chaque course, 0 ailleurs) et l'activité de chaque cellule ; on garde le
%   maximum (décalages -max_lag..+max_lag) et son décalage.
%   Deux changements :
%     1. signal continu : la corrélation est faite sur tout l'enregistrement
%        (les courses ne sont plus collées bout à bout), pour pouvoir
%        chercher jusqu'à 40 s après le départ alors que les courses
%        durent ~8 s ;
%     2. seuil : à la place d'Otsu, test par cellule. Faux débuts : autant
%        de débuts que de vraies courses, tirés au hasard parmi les trames
%        où l'animal court (trames_course), n_faux fois ; même calcul du
%        maximum. p = (1 + nombre de faux >= vrai) / (n_faux + 1). Puis FDR
%        de Benjamini-Hochberg (q_fdr) sur les cellules de la session.
%        Une cellule seulement plus active quand l'animal court ne passe
%        pas : il faut qu'elle soit calée sur le départ.
%   Cellules de la séquence : significatives ET décalage < 0 (la cellule
%   s'active APRÈS le début de course), comme avant.
%   ordre : ces cellules rangées du décalage le plus proche de 0 au plus
%   tardif (comme decalage_cellules).
%   Le générateur aléatoire global doit être fixé avant l'appel
%   (fixer_flux_aleatoire) pour des tirages reproductibles.
%
%   DF : cellules x trames (z-score) ; debuts : trame de début de chaque
%   course ; courses : numéros des courses utilisées ; trames_course :
%   trames où l'animal court (candidats pour les faux débuts).

[NCell, Nz] = size(DF);
L = max_lag;
DFp = [zeros(NCell, L), DF, zeros(NCell, L)];     % zéros aux bords, comme xcorr
norme = sqrt(sum(DF .^ 2, 2));                    % normalisation 'coeff' (partie cellule)
d = debuts(courses);
d = d(:)';
n = numel(d);
lags = -L:L;
corr_max = nan(1, NCell);
decalage = nan(1, NCell);
p = nan(1, NCell);
dans_sequence = false(1, NCell);
ordre = zeros(1, 0);
if n == 0, return; end

% vraie corrélation : pic et décalage
C = somme_blocs(DFp, {d}, L) ./ (norme * sqrt(n));
[corr_max, k] = max(C, [], 2);
corr_max = corr_max';
decalage = lags(k);

% faux débuts : même nombre, tirés parmi les trames de course (par paquets
% de 25 tirages, pour la vitesse)
trames_course = trames_course(:)';
nf = min(n, numel(trames_course));
nb_sup = zeros(NCell, 1);
paquet = 25;
for t0 = 1:paquet:n_faux
    T = min(paquet, n_faux - t0 + 1);
    faux = cell(1, T);
    for t = 1:T
        faux{t} = trames_course(randperm(numel(trames_course), nf));
    end
    Cf = somme_blocs(DFp, faux, L) ./ (norme * sqrt(nf));
    Mf = max(reshape(Cf, NCell, 2 * L + 1, T), [], 2);     % NCell x 1 x T
    nb_sup = nb_sup + sum(reshape(Mf, NCell, T) >= corr_max' - 1e-12, 2);
end
p = ((1 + nb_sup) / (n_faux + 1))';

ok = isfinite(corr_max) & isfinite(p);
sig = false(1, NCell);
if any(ok)
    sig(ok) = bh_fdr_mask(p(ok), q_fdr);
end
dans_sequence = ok & sig & decalage < 0;
idx1 = find(dans_sequence);
[~, idx2] = sort(decalage(idx1), 2, 'descend');
ordre = idx1(idx2);
end


function C = somme_blocs(DFp, D, L)
% Pour chaque jeu de débuts D{t} : C_t(:, j) = somme sur les débuts de
% DF(:, d - lags(j)) avec lags = -L..L, c'est-à-dire la corrélation croisée
% de xcorr(debut_binaire, DF) (non normalisée). Dans DFp (décalé de L),
% DF(:, d - m) = DFp(:, d - m + L) : la colonne d + L - m. Calcul par
% produit avec une matrice creuse (DFp * B), les jeux côte à côte :
% C = [C_1, C_2, ...], chacun NCell x (2L + 1), lags croissants.
nL = 2 * L + 1;
T = numel(D);
lignes = cell(T, 1); colonnes = cell(T, 1);
for t = 1:T
    d = D{t}(:);
    li = d + L - (-L:L);                         % numel(d) x nL : colonne de DFp pour chaque lag
    co = repmat((t - 1) * nL + (1:nL), numel(d), 1);
    lignes{t} = li(:); colonnes{t} = co(:);
end
B = sparse(vertcat(lignes{:}), vertcat(colonnes{:}), 1, size(DFp, 2), T * nL);
C = full(DFp * B);
end
