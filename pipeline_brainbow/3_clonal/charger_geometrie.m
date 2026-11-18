function [outline_gcampx, outline_gcampy, gcamp_centroid, calcium, roi_cells, Coomasqx, Coomasqy] = charger_geometrie(path)
% Contours et centres des cellules suite2p (Fall.mat), et image calcique.
% Repris de brainbowassemblies2026_05_07, mais Fall.mat est chargé dans une
% struct : il n'écrase plus F et ops, et ne masque plus la fonction iscell.
%   path : dossier suite2p, terminé par \ (ou /)
%   Coomasqx, Coomasqy : pixels du masque de chaque cellule (sorties facultatives)
% suite2p compte les pixels à partir de 0 : +1 pour MATLAB (c'est suffisant
% si calcium.tif est l'image moyenne pleine taille ; si c'est max_proj ou
% Vcorr, rognées par ops.xrange / ops.yrange, il faut aussi ce décalage).

S2P = load(strcat(path, 'Fall.mat'), 'stat', 'iscell');
roi_cells = find(S2P.iscell(:, 1) == 1);
NCell = numel(roi_cells);

outline_gcampx = cell(1, NCell);
outline_gcampy = cell(1, NCell);
Coomasqx = cell(1, NCell);
Coomasqy = cell(1, NCell);
for m = 1:NCell
    n = roi_cells(m);
    if iscell(S2P.stat), st = S2P.stat{n}; else, st = S2P.stat(n); end
    outline_gcampx{m} = double(st.xext + 1);
    outline_gcampy{m} = double(st.yext + 1);
    Coomasqx{m} = double(st.xpix + 1);
    Coomasqy{m} = double(st.ypix + 1);
end

% Centre de chaque cellule = moyenne des pixels de son masque
gcamp_centroid = zeros(NCell, 2);
for n = 1:NCell
    gcamp_centroid(n, :) = [mean(Coomasqx{n}), mean(Coomasqy{n})];
end

calcium = [];
if isfile(strcat(path, 'calcium.tif'))
    calcium = imread(strcat(path, 'calcium.tif'));
end
end
