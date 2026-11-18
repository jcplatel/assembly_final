function MAct = activite_cosynchrone(Raster, synchronous_frames)
%ACTIVITE_COSYNCHRONE  MAct(i) = nombre de cellules actives au moins une fois
%   entre les trames i et i + synchronous_frames (même calcul que la boucle
%   de SCE_detection4, en beaucoup plus rapide).
Nz = size(Raster, 2);
if Nz <= synchronous_frames
    MAct = zeros(1, 0);
    return
end
fenetre_max = movmax(Raster, [0 synchronous_frames], 2);
MAct = sum(fenetre_max(:, 1:Nz - synchronous_frames), 1);
end
