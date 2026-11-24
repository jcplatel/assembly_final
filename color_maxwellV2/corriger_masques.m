function [info, stat] = corriger_masques(dossier, taille_cible, enregistrer)
%% =========================================================================
%  CORRIGER_MASQUES (18 nov., demandée par l'utilisateur)
%  Les masques suite2p des animaux de 2025 (446178, 446198) sont trop gros :
%  médiane ~115 pixels contre ~70-80 pour les animaux de 2022-2023, alors que
%  le rayon du cœur (stat.radius) est presque le même. L'excès est un halo de
%  pixels de faible poids (stat.lam) autour de la cellule.
%  Correction : pour chaque ROI, on garde les pixels de plus fort poids, en
%  nombre r x npix, avec r = taille_cible / médiane des npix des cellules
%  (iscell = 1) de la session ; puis on garde le morceau qui contient le pixel
%  le plus fort et on bouche les trous. Les tailles relatives entre cellules
%  sont donc gardées. Si r >= 1, rien ne change.
%  taille_cible : défaut 73 pixels (médiane des animaux nocues de 2022-2023 :
%  411582, 443669, 443695, 444113, 444119, 444152, 444175).
%
%  Fall.mat n'est JAMAIS modifié. enregistrer = true : écrit
%  masques_corriges.mat dans le dossier de la session (stat corrigé pour
%  toutes les ROI, dans le même ordre que Fall.mat, et info). Les fonctions
%  de couleurs le lisent à la place de Fall.mat (charger_masques).
%  =========================================================================
if nargin < 2 || isempty(taille_cible), taille_cible = 73; end
if nargin < 3, enregistrer = false; end
d = regexprep(char(dossier), '[\\/]+$', '');
F = load(fullfile(d, 'Fall.mat'), 'stat', 'iscell', 'ops');
sz = [F.ops.Ly F.ops.Lx];
stat = F.stat; est_cell = iscell(stat);
n = numel(stat);
npix_avant = zeros(n, 1); npix_apres = zeros(n, 1);
for k = 1:n
    if est_cell, st = stat{k}; else, st = stat(k); end
    npix_avant(k) = numel(st.xpix);
end
cel = F.iscell(:, 1) == 1;
r = taille_cible / median(npix_avant(cel));
for k = 1:n
    if est_cell, st = stat{k}; else, st = stat(k); end
    if r < 1
        st = rogner(st, r, sz);
    end
    npix_apres(k) = numel(st.xpix);
    if est_cell, stat{k} = st; else, stat(k) = st; end
end
info = struct('date', datestr(now, 'yyyy-mm-dd HH:MM:SS'), 'taille_cible', taille_cible, ...
    'rapport', min(r, 1), 'npix_avant', npix_avant, 'npix_apres', npix_apres, ...
    'mediane_avant', median(npix_avant(cel)), 'mediane_apres', median(npix_apres(cel)), ...
    'methode', 'pixels de plus fort lam, r x npix, morceau du pixel le plus fort, trous bouchés'); %#ok<DATST,TNOW1>
if enregistrer
    save(fullfile(d, 'masques_corriges.mat'), 'stat', 'info');
end
end


function st = rogner(st, r, sz)
x = double(st.xpix(:)); y = double(st.ypix(:)); lam = double(st.lam(:));
garder = max(5, round(r * numel(x)));
if garder >= numel(x), return; end
[~, o] = sort(lam, 'descend');
bw = false(sz);
idx = sub2ind(sz, y(o(1:garder)) + 1, x(o(1:garder)) + 1);
bw(idx) = true;
bw = imfill(bw, 'holes');
L = bwlabel(bw, 8);
centre = L(sub2ind(sz, y(o(1)) + 1, x(o(1)) + 1));
bw = L == centre;
% on ne garde que des pixels de la ROI d'origine (imfill peut en ajouter)
dans = false(sz); dans(sub2ind(sz, y + 1, x + 1)) = true;
bw = bw & dans;
k = find(bw(sub2ind(sz, y + 1, x + 1)));
st.xpix = st.xpix(k); st.ypix = st.ypix(k); st.lam = st.lam(k);
if isfield(st, 'npix'), st.npix = numel(k); end
end
