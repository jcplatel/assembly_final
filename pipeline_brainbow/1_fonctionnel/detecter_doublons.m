function D = detecter_doublons(Tr1b, WinRest, Coomasqx, Coomasqy, r_min)
%DETECTER_DOUBLONS  ROI qui sont probablement la même cellule comptée deux
%   fois (7 nov., demandé par l'utilisateur après 444175_221121_plane1 : ROI
%   4 et 224, centres à 0,2 µm, traces corrélées à 0,94).
%   Doublon = traces au repos corrélées à plus de r_min (0,7 par défaut)
%   ET masques qui se touchent (au moins un pixel commun) ou centres à moins
%   d'un diamètre de cellule (diamètre = celui d'un disque de même surface
%   que le masque, moyenne des deux ROI).
%   Sur 444175_221121_plane1 : 5 paires à r > 0,7 (centres à 0,2-9,9 µm,
%   masques qui se touchent) ; paire suivante, r = 0,67 ; une paire
%   ordinaire est à r = 0,03.
%   Dans chaque paire, on retire la ROI au plus petit masque (le morceau) ;
%   les paires sont prises de la plus corrélée à la moins corrélée, et une
%   paire dont une ROI est déjà retirée est sautée (chaînes A-B-C).
%   Entrées : Tr1b (cellules x trames), WinRest (trames de repos),
%   Coomasqx / Coomasqy (pixels des masques, charger_geometrie).
%   Sortie D : retirees (indices des ROI retirées), paires (table : roi_1,
%   roi_2, r, distance_px, se_touchent, retiree), r_min.

if nargin < 5 || isempty(r_min), r_min = 0.7; end
NCell = size(Tr1b, 1);
D = struct('retirees', zeros(1, 0), 'paires', table(), 'r_min', r_min);
if NCell < 2 || numel(Coomasqx) ~= NCell
    return
end
X = double(Tr1b(:, WinRest));
X = X - mean(X, 2);
s = sqrt(sum(X .^ 2, 2));
s(s == 0) = Inf;
X = X ./ s;
R = X * X';
R(1:NCell + 1:end) = 0;
[i, j] = find(triu(R > r_min, 1));
if isempty(i)
    return
end
r = R(sub2ind([NCell NCell], i, j));
npix = cellfun(@numel, Coomasqx(:));
cx = cellfun(@mean, Coomasqx(:));
cy = cellfun(@mean, Coomasqy(:));
d = hypot(cx(i) - cx(j), cy(i) - cy(j));
diam = 2 * sqrt(npix / pi);
proche = d < (diam(i) + diam(j)) / 2;
touche = false(numel(i), 1);
for k = 1:numel(i)
    pi_ = round(Coomasqy{i(k)}(:)) + 1e5 * round(Coomasqx{i(k)}(:));
    pj_ = round(Coomasqy{j(k)}(:)) + 1e5 * round(Coomasqx{j(k)}(:));
    touche(k) = any(ismember(pi_, pj_));
end
garde = proche | touche;
i = i(garde); j = j(garde); r = r(garde); d = d(garde); touche = touche(garde);
[r, o] = sort(r, 'descend');
i = i(o); j = j(o); d = d(o); touche = touche(o);
retiree = zeros(numel(i), 1);
deja = false(NCell, 1);
for k = 1:numel(i)
    if deja(i(k)) || deja(j(k)), continue; end
    if npix(j(k)) < npix(i(k)) || (npix(j(k)) == npix(i(k)) && j(k) > i(k))
        retiree(k) = j(k);
    else
        retiree(k) = i(k);
    end
    deja(retiree(k)) = true;
end
D.retirees = sort(retiree(retiree > 0))';
D.paires = table(i, j, r, d, touche, retiree, 'VariableNames', ...
    {'roi_1', 'roi_2', 'r', 'distance_px', 'se_touchent', 'retiree'});
end
