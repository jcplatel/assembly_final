function [corr_max, decalage, seuil, dans_sequence, ordre] = decalage_cellules(DF, debuts, fins, courses, max_lag)
%DECALAGE_CELLULES  Place de chaque cellule par rapport au début de course
%   (code de l'utilisateur, sequence_2026_01_27, sans changement).
%   Sur les courses choisies seulement (trames debut:fin de chaque course,
%   mises bout à bout) : corrélation croisée entre les débuts de course
%   (1 au début de chaque course, 0 ailleurs) et l'activité de chaque
%   cellule, décalage jusqu'à max_lag trames ; on garde le maximum et son
%   décalage.
%   Cellules de la séquence : maximum au-dessus du seuil d'Otsu (graythresh)
%   ET décalage < 0 (la cellule s'active APRÈS le début de course).
%   ordre : ces cellules rangées du décalage le plus proche de 0 au plus
%   tardif (sort(..., 'descend') comme dans le code d'origine).
%
%   DF : cellules x trames ; debuts, fins : trames de chaque course ;
%   courses : numéros des courses à utiliser.

[NCell, Nz] = size(DF);
W = cell2mat(arrayfun(@(n) debuts(n):fins(n), courses(:)', 'UniformOutput', false));
debut_binaire = zeros(1, Nz);
debut_binaire(debuts) = 1;
x = debut_binaire(W);
corr_max = nan(1, NCell);
decalage = nan(1, NCell);
for i = 1:NCell
    [c, lags] = xcorr(x, DF(i, W), max_lag, 'coeff');
    [corr_max(i), k] = max(c);
    decalage(i) = lags(k);
end
ok = isfinite(corr_max);
seuil = graythresh(corr_max(ok));
dans_sequence = ok & corr_max > seuil & decalage < 0;
idx1 = find(dans_sequence);
[~, idx2] = sort(decalage(idx1), 2, 'descend');
ordre = idx1(idx2);
end
