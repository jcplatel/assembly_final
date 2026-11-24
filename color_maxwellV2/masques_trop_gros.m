function [rejet, info] = masques_trop_gros(dossier, facteur, allongement_max)
%% =========================================================================
%  MASQUES_TROP_GROS (17 nov., demandée par l'utilisateur)
%  Masques suite2p (iscell = 1, dans l'ordre de colorcell) à passer en noir (8) :
%   - trop gros : plus de facteur x la taille médiane des masques de la
%     session (défaut 2,5) : plusieurs cellules collées ou du neuropile. La
%     taille est comparée à la médiane de la session, car les masques des
%     animaux de 2025 sont plus gros (réglages suite2p).
%   - trop allongés : grand axe / petit axe > allongement_max (défaut 2,5),
%     mesuré sur le plus grand morceau du masque (regionprops).
%  Les masques en virgule ou en crochet (solidité faible) ne sont PAS retirés
%  (choix de l'utilisateur : ce peut être un prolongement ; correction à la
%  main).
%  =========================================================================
if nargin < 2 || isempty(facteur), facteur = 2.5; end
if nargin < 3 || isempty(allongement_max), allongement_max = 2.5; end
d = regexprep(char(dossier), '[\\/]+$', '');
F = load(fullfile(d, 'Fall.mat'), 'ops');
[F.stat, F.iscell] = charger_masques(d);   % 18 nov. : masques corrigés s'ils existent
sz = [F.ops.Ly F.ops.Lx];
roi = find(F.iscell(:, 1) == 1);
npix = zeros(numel(roi), 1); allong = zeros(numel(roi), 1);
for j = 1:numel(roi)
    if iscell(F.stat), st = F.stat{roi(j)}; else, st = F.stat(roi(j)); end
    npix(j) = numel(st.xpix);
    bw = false(sz);
    bw(sub2ind(sz, double(st.ypix(:)) + 1, double(st.xpix(:)) + 1)) = true;
    p = regionprops(bw, 'Area', 'MajorAxisLength', 'MinorAxisLength');
    [~, g] = max([p.Area]);
    allong(j) = p(g).MajorAxisLength / max(p(g).MinorAxisLength, 1);
end
gros = npix > facteur * median(npix);
long = allong > allongement_max;
rejet = gros | long;
info = struct('facteur', facteur, 'allongement_max', allongement_max, 'npix', npix, ...
    'npix_median', median(npix), 'allongement', allong, 'id_gros', find(gros), 'id_allonges', find(long), ...
    'id', find(rejet));
end
