function [colorcell, details] = colorcell_auto(dossier, varargin)
%% =========================================================================
%  COLORCELL_AUTO (17 nov., demandée par l'utilisateur)
%  Les étapes 1 à 7 de interface_couleurs / color_predictor_model, sans
%  fenêtre, avec ce que les fenêtres proposent quand on valide sans rien
%  changer :
%    1-2  lecture, retrait du fond, normalisation, recalage sur le GCaMP
%    3    valeurs des ROI (ROI_extraction2)
%    4    rejet des cellules faibles : sqrt(R²+V²+B²) < seuil (0,25), comme
%         interactive_threshold_selection sans clic ; ou la liste
%         id_low_intensity.mat du dossier si 'id_low', 'fichier'
%    5    correction de la fuite : k et b estimés par bleedthrough_correction3
%         (même calcul, copié ci-dessous), sans retouche
%    6    rien de plus (sans clic, l'étape 6 reprend la liste de l'étape 4)
%    7    maxwell_proportions + color_classification
%  Rien n'est écrit : voir enregistrer_colorcell_auto.
%
%  Options (nom, valeur) :
%    'seuil'       0.25        seuil des étapes 4 et 6
%    'id_low'      'seuil'     'seuil' ou 'fichier' (id_low_intensity.mat)
%    'ep'          0           bande incertaine (color_classification)
%    'seuil_blanc' 0.20        Inf = pas de classe blanche
%    'methode'     'roi95'     échelle de maxwell_proportions ('ancien')
%    'limits_k'    [0.4 0.2 0.1]  limites de k [vert rouge bleu]
%  =========================================================================

p = inputParser;
p.addParameter('seuil', 0.25);
p.addParameter('id_low', 'seuil');
p.addParameter('ep', 0);
p.addParameter('seuil_blanc', 0.20);
p.addParameter('methode', 'roi95');
p.addParameter('limits_k', [0.4, 0.2, 0.1]);
p.addParameter('pre', []);          % sortie 'details' d'un appel précédent : saute les étapes 1-3
p.parse(varargin{:});
o = p.Results;
pre = o.pre; o = rmfield(o, 'pre');

d = [regexprep(char(dossier), '[\\/]+$', '') filesep];

if isempty(pre)
    % 1-2
    [red, red890, green, blue, calcium] = brainbow_image_read(d);
    In = brainbow_image_processing(d, red, red890, green, blue, calcium);
    [ar, ag, ab, cn, ~, info_alignement] = image_alignment_multimodal(In);
    % 3
    [R, G, B, Ca] = ROI_extraction2(d, ar, ag, ab, In.calcium);
else
    ar = pre.images_alignees{1}; ag = pre.images_alignees{2}; ab = pre.images_alignees{3};
    cn = pre.images_alignees{4}; info_alignement = pre.info_alignement;
    R = pre.roi_red; G = pre.roi_green; B = pre.roi_blue; Ca = pre.roi_ca;
end
n = numel(R);

% 4
if strcmpi(o.id_low, 'fichier') && isfile([d 'id_low_intensity.mat'])
    X = load([d 'id_low_intensity.mat'], 'id_low_intensity');
    bas = false(n, 1);
    v = X.id_low_intensity(:);
    if islogical(v) && numel(v) == n
        bas = v;
    else
        v = v(v >= 1 & v <= n);
        bas(v) = true;
    end
    source_bas = 'id_low_intensity.mat';
else
    bas = sqrt(R.^2 + G.^2 + B.^2) < o.seuil;
    source_bas = sprintf('seuil %.2f', o.seuil);
end
R1 = R; G1 = G; B1 = B;
R1(bas) = NaN; G1(bas) = NaN; B1(bas) = NaN;

% 5 : mêmes formules que bleedthrough_correction3, k et b estimés
[k_g, b_g] = estimer_fuite(G1, Ca, o.limits_k(1));
[k_r, b_r] = estimer_fuite(R1, Ca, o.limits_k(2));
[k_b, b_b] = estimer_fuite(B1, Ca, o.limits_k(3));
G2 = max(0, G1 - (k_g * Ca + b_g));
R2 = max(0, R1 - (k_r * Ca + b_r));
B2 = max(0, B1 - (k_b * Ca + b_b));
img_g = max(0, ag - (k_g * cn + b_g));
img_r = max(0, ar - (k_r * cn + b_r));
img_b = max(0, ab - (k_b * cn + b_b));
% max(0, NaN) = 0 : les cellules rejetées redeviennent NaN, comme dans le pipeline
R2(bas) = NaN; G2(bas) = NaN; B2(bas) = NaN;

% 7
[D, r, g, b, info_maxwell] = maxwell_proportions(R2, G2, B2, o.methode);
[colorcell, marge] = color_classification(D, n, r, g, b, o.ep, o.seuil_blanc);
colorcell = colorcell(:);

details = struct();
details.date = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<DATST,TNOW1>
details.methode = 'colorcell_auto (étapes 1-7 sans fenêtre, color_maxwellV2)';
details.options = o;
details.colorcell = colorcell;
details.id_low_intensity = find(bas);
details.source_id_low = source_bas;
details.roi_red = R; details.roi_green = G; details.roi_blue = B; details.roi_ca = Ca;
details.roi_red_corrected3 = R2; details.roi_green_corrected3 = G2; details.roi_blue_corrected3 = B2;
details.r = r; details.g = g; details.b = b; details.marge = marge;
details.Donnees_Clustering = D;
details.info_maxwell = info_maxwell;
details.info_alignement = info_alignement;
details.params_fuite = struct('k', [k_g, k_r, k_b], 'b', [b_g, b_r, b_b], ...
    'limites_k', o.limits_k, 'ordre', 'vert rouge bleu');
details.images_corrigees = {img_r, img_g, img_b};     % non enregistrées dans le .mat
details.images_alignees = {ar, ag, ab, cn};            % idem
end


function [k_est, b_est] = estimer_fuite(roi_col, roi_ca, k_limit)
% copie de estimate_params (bleedthrough_correction3) : élastique inférieur
idx_v = (roi_ca > 0) & (roi_col > 0) & isfinite(roi_ca) & isfinite(roi_col);
x_v = roi_ca(idx_v); y_v = roi_col(idx_v);
if length(x_v) < 10
    k_est = 0; b_est = 0; return
end
x_aug = [x_v; max(x_v)]; y_aug = [y_v; max(y_v) * 10];
try
    k_h = convhull(x_aug, y_aug);
    x_h = x_aug(k_h); y_h = y_aug(k_h);
    [~, i_min_x] = min(x_h); i_max_x = max(x_v);
    x_l = []; y_l = [];
    for i = 1:length(x_h)
        if x_h(i) >= x_h(i_min_x) && x_h(i) <= i_max_x && y_h(i) <= max(y_v)
            x_l = [x_l; x_h(i)]; y_l = [y_l; y_h(i)]; %#ok<AGROW>
        end
    end
    if length(x_l) < 2
        x_l = [min(x_v); max(x_v)]; y_l = [min(y_v); min(y_v)];
    end
catch
    x_l = [min(x_v); max(x_v)]; y_l = [min(y_v); min(y_v)];
end
if length(x_l) >= 2
    pp = robustfit(x_l, y_l); b_est = pp(1); k_est = pp(2);
else
    k_est = 0; b_est = min(y_v);
end
if k_est < 0, k_est = 0; end
b_fit_orig = b_est; k_brid = false;
if k_est > k_limit, k_est = k_limit; k_brid = true; end
y_virt = y_v - k_est * x_v; b_est = prctile(y_virt, 5);
if k_brid
    b_max = min([b_fit_orig, max(y_v) * 0.05]);
    if b_est > b_max, b_est = b_max; end
end
if b_est < 0, b_est = 0; end
if b_est > min(y_v), b_est = min(y_v); end
end
