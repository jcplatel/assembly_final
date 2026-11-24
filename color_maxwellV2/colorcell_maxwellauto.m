function [colorcell, details] = colorcell_maxwellauto(dossier, enregistrer, marge_min, suffixe)
%% =========================================================================
%  COLORCELL_MAXWELLAUTO (17 nov., demandée par l'utilisateur)
%  Couleurs d'une session sans intervention :
%    1-3  comme interface_couleurs : brainbow_image_read,
%         brainbow_image_processing, image_alignment_multimodal (recalage
%         sur le GCaMP), ROI_extraction2
%    4    rejet des cellules faibles : max(R, V, B) < 0,20
%    5    correction de la fuite du GCaMP : k et b estimés comme
%         bleedthrough_correction3 (estimer_fuite3), k borné par
%         [vert rouge bleu] = [0,4 0,2 0,1]
%    6    rejet : max(R, V, B) corrigés < 0,15
%    7    fond retiré = médiane de l'image corrigée hors des ROI ; poids du
%         bleu 0,7 ; proportions + color_classification (ep = 0, pas de
%         classe blanche, frontière du bleu pur 0,4) ; noire (8) si la
%         couleur est proche d'une frontière du triangle (marge < marge_min,
%         défaut 0,05)
%  Validation (17 nov., GT 444119 220919/220923, 444175 221125) : bonne
%  couleur pour 83-91 % des cellules colorées gardées ; 35-63 % des noires
%  du GT retrouvées ; 79-90 % des cellules colorées du GT gardées.
%  ATTENTION : les cellules hors de la couche pyramidale ou absentes du plan
%  ex vivo ne sont pas toutes retirées : vérifier chaque session à l'œil.
%
%  enregistrer = true : écrit dans le dossier de la session
%    colorcell_MaxwellAuto.mat   (variable colorcell)
%    couleurs_details_Auto.mat   (valeurs, proportions, marges, rejets...)
%  Rien n'est écrasé : si colorcell_MaxwellAuto.mat existe, erreur.
%  18 nov. : suffixe (facultatif) ajouté aux deux noms de fichier, par ex.
%  '_masquescorr' -> colorcell_MaxwellAuto_masquescorr.mat. Les masques sont
%  lus par charger_masques : masques_corriges.mat s'il existe, sinon Fall.mat.
%  =========================================================================

if nargin < 2, enregistrer = false; end
if nargin < 3 || isempty(marge_min), marge_min = 0.05; end
if nargin < 4, suffixe = ''; end
d = [regexprep(char(dossier), '[\\/]+$', '') filesep];
f_couleur = [d 'colorcell_MaxwellAuto' suffixe '.mat'];
f_details = [d 'couleurs_details_Auto' suffixe '.mat'];
if enregistrer && isfile(f_couleur)
    error('%s existe déjà (rien n''est écrasé).', f_couleur);
end
limites_k = [0.4, 0.2, 0.1];            % [vert rouge bleu], comme color_predictor_model

% 1-3
[red, red890, green, blue, calcium] = brainbow_image_read(d);
In = brainbow_image_processing(d, red, red890, green, blue, calcium);
[ar, ag, ab, cn, ~, info_alignement] = image_alignment_multimodal(In);
% 17 nov. : chaque canal recalé sur les cellules si besoin (ajuster_canaux_masques)
[ar, ag, ab, info_alignement.ajustement_masques] = ajuster_canaux_masques(ar, ag, ab, d);
[R, G, B, Ca] = ROI_extraction2(d, ar, ag, ab, In.calcium);
R = R(:); G = G(:); B = B(:); Ca = Ca(:);

% 4
bas1 = max([R, G, B], [], 2) < 0.20;
R1 = R; G1 = G; B1 = B;
R1(bas1) = NaN; G1(bas1) = NaN; B1(bas1) = NaN;
% 5
[kg, bg] = estimer_fuite3(G1, Ca, limites_k(1));
[kr, br] = estimer_fuite3(R1, Ca, limites_k(2));
[kb, bb] = estimer_fuite3(B1, Ca, limites_k(3));
G2 = max(0, G1 - (kg * Ca + bg)); R2 = max(0, R1 - (kr * Ca + br)); B2 = max(0, B1 - (kb * Ca + bb));
R2(bas1) = NaN; G2(bas1) = NaN; B2(bas1) = NaN;
% 6
bas2 = ~bas1 & max([R2, G2, B2], [], 2) < 0.15;
R2(bas2) = NaN; G2(bas2) = NaN; B2(bas2) = NaN;
rejet = bas1 | bas2;
% 7 : fond de l'image hors ROI, poids du bleu 0,7, frontières [0,4 0,3 0,4] (proportions_auto)
img_c = {max(0, ar - (kr * cn + br)), max(0, ag - (kg * cn + bg)), max(0, ab - (kb * cn + bb))};
[D, r, g, b, classe, marge, info_maxwell] = proportions_auto(R2, G2, B2, img_c, d, rejet);
fond_image = info_maxwell.fond;
colorcell = classe;
colorcell(isnan(marge) | marge < marge_min) = 8;
% 17 nov. : masques trop gros (> 2,5 x la médiane de la session) ou trop allongés (axes > 2,5) en noir
[gros, info_gros] = masques_trop_gros(d);
colorcell(gros) = 8;

details = struct();
details.date = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<DATST,TNOW1>
details.methode = ['colorcell_maxwellauto : rejet max(R,V,B) 0,20 / 0,15, fuite bleedthrough_correction3 ' ...
    '(k et b estimés, bornés), fond de l''image hors ROI, poids du bleu 0,7, frontières [0,4 0,3 0,4], ' ...
    'sans blanc (color_maxwellV2)'];
details.marge_min = marge_min;
details.colorcell = colorcell;
details.classe_triangle = classe;          % avant le passage en noir des rejetées et des marges faibles
details.marge = marge;
details.id_low_intensity = find(rejet);
details.rejet_etape4 = find(bas1); details.rejet_etape6 = find(bas2);
details.roi_red = R; details.roi_green = G; details.roi_blue = B; details.roi_ca = Ca;
details.roi_red_corrected3 = R2; details.roi_green_corrected3 = G2; details.roi_blue_corrected3 = B2;
details.r = r; details.g = g; details.b = b;
details.Donnees_Clustering = D;            % pour le bouton « Triangle » de interface_couleurs
details.info_maxwell = info_maxwell;
details.params_fuite = struct('k', [kg kr kb], 'b', [bg br bb], 'limites_k', limites_k, 'ordre', 'vert rouge bleu');
details.info_alignement = info_alignement;

details.fond_image = fond_image;           % [R V B], retiré avant les proportions
details.masques_trop_gros = info_gros;     % 17 nov. : passés en noir
details.masques_corriges = isfile([d 'masques_corriges.mat']);   % 18 nov.

if enregistrer
    save(f_couleur, 'colorcell');
    save(f_details, '-struct', 'details');
end
end
