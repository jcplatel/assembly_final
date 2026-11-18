function [seuils, N_obs, N_hasard] = courbe_SCE_hasard(Raster, repos, opts, decalage_max, nNull, graine)
%COURBE_SCE_HASARD  SCE détectés au-dessus de chaque seuil : données et hasard.
%
%   seuils   : seuil minimum (opts.sce_n_cells_threshold) jusqu'au plus grand pic
%   N_obs    : nombre de SCE détectés dans les données au-dessus de chaque seuil
%   N_hasard : nombre moyen de SCE détectés dans nNull surrogates où le train
%              de chaque cellule est décalé d'un bloc de k trames, k tiré entre
%              -decalage_max et +decalage_max (Inf = décalage circulaire
%              complet). L'activité de chaque cellule est gardée, la synchronie
%              plus fine que le décalage est détruite.
%
%   « Au-dessus du seuil s » = pic de MAct STRICTEMENT supérieur à s, comme
%   findpeaks(..., 'MinPeakHeight', s) dans SCE_detection5.
%
%   Les décalages se font À L'INTÉRIEUR des trames de repos mises bout à bout
%   (repos = masque logique) : les transients ne sont détectés qu'au repos.
%   Décaler sur toute la session étalerait l'activité sur les périodes de
%   course, où le raster est vide, et sous-estimerait le hasard.

synchronous_frames = opts.synchronous_frames;
MinPeakDistancesce = opts.MinPeakDistancesce;
seuil_min          = opts.sce_n_cells_threshold;
if isfield(opts, 'max_cells_SCE'), max_cells = opts.max_cells_SCE; else, max_cells = 150; end

R = Raster(:, logical(repos));          % trames de repos mises bout à bout
[NCell, Nr] = size(R);
pics_obs = pics_SCE(activite_cosynchrone(R, synchronous_frames), seuil_min, MinPeakDistancesce, max_cells);
seuils = seuil_min:max([seuil_min; pics_obs(:)]);
N_obs = arrayfun(@(s) sum(pics_obs > s), seuils);

N_hasard = zeros(size(seuils));
if Nr <= 2 * synchronous_frames
    return
end
for n = 1:nNull
    fixer_flux_aleatoire(graine, n, nNull);
    Rn = zeros(NCell, Nr);
    for c = 1:NCell
        if isinf(decalage_max)
            k = randi(Nr);
        else
            k = randi([-decalage_max decalage_max]);
        end
        Rn(c, :) = circshift(R(c, :), k);
    end
    p = pics_SCE(activite_cosynchrone(Rn, synchronous_frames), seuil_min, MinPeakDistancesce, max_cells);
    N_hasard = N_hasard + arrayfun(@(s) sum(p > s), seuils) / nNull;
end
end


function pks = pics_SCE(MAct, seuil_min, MinPeakDistancesce, max_cells)
% Pics de MAct au-dessus de seuil_min, avec la même règle que SCE_detection5.
% Avec MinPeakDistance, findpeaks garde les pics du plus haut au plus bas :
% les pics au-dessus d'un seuil s plus élevé sont donc exactement ceux qu'on
% obtiendrait en relançant findpeaks avec s.
if numel(MAct) <= MinPeakDistancesce + 1 || max(MAct) <= seuil_min
    pks = zeros(1, 0);
    return
end
pks = findpeaks(double(MAct), 'MinPeakHeight', seuil_min, 'MinPeakDistance', MinPeakDistancesce);
pks = pks(pks <= max_cells);
end
