function [X, valide, info] = variables_couleurs(R, G, B, Ca)
%% =========================================================================
%  VARIABLES_COULEURS (17 nov., demandée par l'utilisateur)
%  Prétraitement et variables de colorcell_maxwell4.mlx / prediction_color.m
%  (code de l'utilisateur, copié ici à l'identique) :
%    1. rejet : max(R, V, B) < 0,20 -> NaN       (removal_low_intensities2)
%    2. correction de la fuite du GCaMP, sans limite sur k
%                                                 (bleedthrough_correction)
%    3. rejet : max(R, V, B) corrigés < 0,15 -> NaN
%    4. variables : R, V, B, Ca, I = sqrt(R²+V²+B²), R/(V+0,01),
%       R/(B+0,01), V/(B+0,01)                    (prediction_color)
%  Entrées : valeurs des ROI de ROI_extraction2 (vecteurs colonnes).
%  Sorties : X (n x 8), valide (cellules non rejetées), info (k, b, rejets).
%  =========================================================================

R = R(:); G = G(:); B = B(:); Ca = Ca(:);
bas1 = max([R, G, B], [], 2) < 0.20;
R(bas1) = NaN; G(bas1) = NaN; B(bas1) = NaN;
[Gc, kg, bg] = fuite(G, Ca);
[Rc, kr, br] = fuite(R, Ca);
[Bc, kb, bb] = fuite(B, Ca);
bas2 = max([Rc, Gc, Bc], [], 2) < 0.15;
Rc(bas2) = NaN; Gc(bas2) = NaN; Bc(bas2) = NaN;
I = sqrt(Rc.^2 + Gc.^2 + Bc.^2);
X = [Rc, Gc, Bc, Ca, I, Rc ./ (Gc + 0.01), Rc ./ (Bc + 0.01), Gc ./ (Bc + 0.01)];
valide = all(isfinite(X), 2);
info = struct('k_fuite', [kg kr kb], 'b_fuite', [bg br bb], 'ordre', 'vert rouge bleu', ...
    'rejet_1', find(bas1), 'rejet_2', find(bas2 & ~bas1), ...
    'noms_variables', {{'R', 'G', 'B', 'Ca', 'I', 'R/G', 'R/B', 'G/B'}});
end


function [y, k, b] = fuite(roi_color, roi_ca)
% bleedthrough_correction de l'utilisateur (minimum par tranche, fit robuste)
idx = (roi_ca > 0) & (roi_color > 0);
x_val = roi_ca(idx); y_val = roi_color(idx);
k = 0; b = 0;
if numel(x_val) < 10, y = roi_color; return; end
edges = linspace(0, max(x_val), 13); x_min = []; y_min = [];
for i = 1:12
    dans = (x_val >= edges(i)) & (x_val < edges(i + 1));
    if sum(dans) > 3
        [mn, im] = min(y_val(dans)); xi = x_val(dans);
        x_min = [x_min; xi(im)]; y_min = [y_min; mn]; %#ok<AGROW>
    end
end
garde = y_min <= median(y_min);
if sum(garde) >= 2
    bk = robustfit(x_min(garde), y_min(garde)); b = bk(1); k = bk(2);
end
if k < 0, k = 0; b = min(y_min); end
y = roi_color - (k * roi_ca + b);
y(y < 0) = 0;
end
