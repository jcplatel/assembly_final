function [D, r, g, b, classe, marge, info] = proportions_auto(R2, G2, B2, images_corrigees, dossier, rejet)
%% =========================================================================
%  PROPORTIONS_AUTO (17 nov., demandée par l'utilisateur)
%  Étape 7 de la méthode MaxwellAuto, commune à colorcell_maxwellauto et à
%  l'étape 7 de interface_couleurs :
%    - fond retiré à chaque canal = médiane de l'image corrigée de la fuite
%      HORS des ROI de suite2p (toutes, élargies de 2 pixels)
%    - poids du bleu 0,7 (le bleu, canal le plus faible, prenait sinon des
%      cellules magenta et cyan)
%    - proportions r, v, b, triangle de Maxwell, color_classification sans
%      classe blanche, frontières pur / mélange [R V B] = [0,4 0,3 0,4]
%  Entrées : valeurs des ROI après fuite et seuils (NaN = rejetées),
%            images corrigées {rouge, vert, bleu}, dossier de la session
%            (pour Fall.mat), rejet (logique, facultatif).
%  Sorties : D (triangle, NaN pour les rejetées), r g b, classe (avant le
%            passage en noir des marges faibles), marge, info (fond, poids...).
%  =========================================================================

poids = [1, 1, 0.7];                    % [R V B]
seuils = [0.4, 0.3, 0.4];               % frontières pur / mélange [R V B]
R2 = R2(:); G2 = G2(:); B2 = B2(:);
if nargin < 6 || isempty(rejet), rejet = ~isfinite(R2) | ~isfinite(G2) | ~isfinite(B2); end

d = [regexprep(char(dossier), '[\\/]+$', '') filesep];
taille = size(images_corrigees{1});
hors_roi = ~masque_toutes_roi(d, taille);
fond = cellfun(@(I) median(I(hors_roi)), images_corrigees(1:3));

net = max(0, [R2, G2, B2] - fond) .* poids;
somme = sum(net, 2);
P3 = net ./ (somme + 1e-6);
r = P3(:, 1); g = P3(:, 2); b = P3(:, 3);
D = [r + b * 0.5, b * (sqrt(3) / 2)];
[classe, marge] = color_classification(D, numel(r), r, g, b, 0, Inf, seuils);
classe = classe(:);
sans = rejet(:) | ~(somme > 0);
classe(sans) = 8;
D(sans, :) = NaN;
info = struct('methode', 'fond = médiane de l''image corrigée hors ROI ; poids [1 1 0.7] ; seuils [0.4 0.3 0.4]', ...
    'fond', fond, 'poids', poids, 'seuils', seuils, 'somme', somme);
end


function M = masque_toutes_roi(d, taille)
% pixels de toutes les ROI de suite2p (cellules ou non), élargis de 2 pixels
F = load([d 'Fall.mat'], 'stat');
M = false(taille);
for n = 1:numel(F.stat)
    if iscell(F.stat), st = F.stat{n}; else, st = F.stat(n); end
    x = double(st.xpix(:)) + 1; y = double(st.ypix(:)) + 1;
    ok = x >= 1 & x <= taille(2) & y >= 1 & y <= taille(1);
    M(sub2ind(taille, y(ok), x(ok))) = true;
end
M = imdilate(M, strel('disk', 2));
end
