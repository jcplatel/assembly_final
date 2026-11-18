function A = alignement_composite(fichier_image, dossier_s2p, decalage_max)
%ALIGNEMENT_COMPOSITE  L'image Brainbow est-elle bien posée sur les masques
%   suite2p ? (11 nov., après 443669_220615_plane0 : image décalée d'environ
%   150 pixels sans qu'aucun contrôle ne le voie)
%   Une ROI bien placée tombe sur une cellule, plus brillante que le reste de
%   l'image. On déplace l'image (de -decalage_max à +decalage_max pixels,
%   défaut 80) et on cherche le décalage qui met le plus de lumière dans les
%   masques (corrélation croisée par FFT). Bien aligné : meilleur décalage
%   proche de 0.
%   Image : canaux 1-3 = rouge, vert, bleu (somme, chaque canal ramené à son
%   99,5e percentile) ; canal 4 = calcium (contrôle : doit être à 0).
%   Sorties (A) :
%     decalage_bb, contraste_bb_0, contraste_bb_max : meilleur décalage [x y]
%       (pixels) de l'image Brainbow ; lumière moyenne dans les masques /
%       moyenne de l'image, sans décalage et au meilleur décalage ;
%     decalage_ca, contraste_ca_0 : pareil pour le canal calcium ;
%     bord_vide_px : largeur du plus grand bord où un canal Brainbow vaut 0
%       (trace d'un recalage par grande translation) ;
%     n_canaux, taille.
if nargin < 3 || isempty(decalage_max), decalage_max = 80; end
A = struct('decalage_bb', [NaN NaN], 'contraste_bb_0', NaN, 'contraste_bb_max', NaN, ...
    'decalage_ca', [NaN NaN], 'contraste_ca_0', NaN, 'bord_vide_px', NaN, 'n_canaux', 0, 'taille', [NaN NaN]);
info = imfinfo(fichier_image);
C = {};
for pg = 1:numel(info)
    I = double(imread(fichier_image, pg));
    for c = 1:size(I, 3), C{end+1} = I(:, :, c); end %#ok<AGROW>
end
A.n_canaux = numel(C);
if A.n_canaux < 3, return; end
[Ly, Lx] = size(C{1});
A.taille = [Ly Lx];

S2P = load(fullfile(dossier_s2p, 'Fall.mat'), 'stat', 'iscell');
M = zeros(Ly, Lx);
for n = find(S2P.iscell(:, 1) == 1)'
    if iscell(S2P.stat), st = S2P.stat{n}; else, st = S2P.stat(n); end
    x = double(st.xpix(:)) + 1; y = double(st.ypix(:)) + 1;
    ok = x >= 1 & x <= Lx & y >= 1 & y <= Ly;
    M(sub2ind([Ly Lx], y(ok), x(ok))) = 1;
end
if ~any(M(:)), return; end

norm = @(I) I / max(prctile(I(:), 99.5), eps);
B = norm(C{1}) + norm(C{2}) + norm(C{3});
[A.decalage_bb, A.contraste_bb_0, A.contraste_bb_max] = meilleur_decalage(B, M, decalage_max);
if A.n_canaux >= 4
    [A.decalage_ca, A.contraste_ca_0] = meilleur_decalage(norm(C{4}), M, decalage_max);
end
% bords vides : lignes ou colonnes entières à 0 dans un canal Brainbow, depuis un bord
bord = 0;
for c = 1:3
    Z = C{c} == 0;
    lz = all(Z, 2); cz = all(Z, 1);
    bord = max([bord, longueur_debut(lz), longueur_debut(flipud(lz)), longueur_debut(cz'), longueur_debut(flipud(cz'))]);
end
A.bord_vide_px = bord;
end


function [d, c0, cmax] = meilleur_decalage(I, M, dmax)
% lumière moyenne de I sous les masques M déplacés de d = [dx dy]
% (masque en (y, x) lu dans I en (y + dy, x + dx)), divisée par la moyenne de I
[Ly, Lx] = size(I);
P = [2 * Ly, 2 * Lx];
FI = fft2(I, P(1), P(2));
FU = fft2(ones(Ly, Lx), P(1), P(2));
FM = conj(fft2(M, P(1), P(2)));
somme = real(ifft2(FI .* FM));            % somme de I sous les masques déplacés
nb = real(ifft2(FU .* FM));               % pixels des masques encore dans l'image
moy = mean(I(:));
best = -Inf; d = [0 0]; c0 = NaN;
for dy = -dmax:dmax
    iy = mod(dy, P(1)) + 1;
    for dx = -dmax:dmax
        ix = mod(dx, P(2)) + 1;
        if nb(iy, ix) < 0.5 * sum(M(:)), continue; end
        v = somme(iy, ix) / nb(iy, ix) / moy;
        if dx == 0 && dy == 0, c0 = v; end
        if v > best, best = v; d = [dx dy]; end
    end
end
cmax = best;
end


function n = longueur_debut(v)
% nombre de valeurs vraies au début de v
n = find(~v(:), 1) - 1;
if isempty(n), n = numel(v); end
end
