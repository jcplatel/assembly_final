function T = tailles_roi(Coomasqx, min_px, max_facteur)
%TAILLES_ROI  ROI trop petites ou trop grandes pour être un corps cellulaire
%   (7 nov., demandé par l'utilisateur : « enlève les minuscules (1, 2 ou 3
%   pixels si ça existe) et les géantes »).
%   Trop petite : moins de min_px pixels (4 par défaut : retire 1 à 3 px).
%   Trop grande : plus de max_facteur fois la surface médiane des ROI de la
%   session (3 par défaut ; relatif, pour ne pas dépendre du zoom).
%   Sortie T : npix (pixels de chaque ROI), mediane, petites, grandes,
%   retirees (les deux), min_px, max_facteur.

if nargin < 2 || isempty(min_px), min_px = 4; end
if nargin < 3 || isempty(max_facteur), max_facteur = 3; end
npix = cellfun(@numel, Coomasqx(:))';
med = median(npix);
petites = find(npix < min_px);
grandes = find(npix > max_facteur * med);
T = struct('npix', npix, 'mediane', med, 'petites', petites, 'grandes', grandes, ...
    'retirees', sort([petites, grandes]), 'min_px', min_px, 'max_facteur', max_facteur);
end
