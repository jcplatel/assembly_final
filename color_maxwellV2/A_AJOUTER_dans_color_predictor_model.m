%% =========================================================================
%  LIGNES À CHANGER / AJOUTER dans color_predictor_model.mlx (11 nov.)
%  Les 4 fonctions modifiées gardent leurs sorties d'avant : ton script
%  marche sans rien changer. Ces lignes servent seulement à récupérer et
%  enregistrer les informations en plus. Rien ne change dans colorcell.
%  =========================================================================

% 1. Recalage : ajouter la 6e sortie
[aligned_red, aligned_green, aligned_blue, calcium_norm, aligned_red890, info_alignement] = ...
    image_alignment_multimodal(images_norm);

% 2. Extraction : ajouter la 11e sortie
[roi_red, roi_green, roi_blue, roi_ca, Coomasqx, Coomasqy, Coox, Cooy, outline_gcampx, outline_gcampy, extra_roi] = ...
    ROI_extraction2(pathload, aligned_red, aligned_green, aligned_blue, images_norm.calcium);

% 3. Fuite : ajouter la 7e sortie
[roi_green_corrected2, roi_red_corrected2, roi_blue_corrected2, ...
 aligned_green_corrected, aligned_red_corrected, aligned_blue_corrected, params_fuite] = ...
    bleedthrough_correction3(roi_green_corrected, roi_red_corrected, roi_blue_corrected, ...
                             roi_ca, limits_k, aligned_green, aligned_red, aligned_blue, calcium_norm);

% 4. Proportions et classification (corrections du 11 nov.)
%    maxwell_proportions : nouvelle échelle par défaut ('roi95') ; mettre 'ancien'
%    en 4e argument pour retrouver le calcul d'avant.
%    color_classification : classe blanc (7) si la plus faible des 3 proportions
%    dépasse 0.20 (7e argument ; Inf = pas de blanc, comme avant).
[Donnees_Clustering, r, g, b, info_maxwell] = maxwell_proportions(roi_red_corrected3, roi_green_corrected3, roi_blue_corrected3);
[colorcell_classification, marge_couleur] = color_classification(Donnees_Clustering, NCell, r, g, b, ep);

%    Pour comparer avec l'ancien calcul sur une session, avant de tout refaire :
[D_anc, r_anc, g_anc, b_anc] = maxwell_proportions(roi_red_corrected3, roi_green_corrected3, roi_blue_corrected3, 'ancien');
colorcell_ancien = color_classification(D_anc, NCell, r_anc, g_anc, b_anc, ep, Inf);
fprintf('Classes changées par les corrections : %d / %d cellules (%.0f %%) ; blanches : %d\n', ...
    sum(colorcell_ancien ~= colorcell_classification), NCell, ...
    100 * mean(colorcell_ancien ~= colorcell_classification), sum(colorcell_classification == 7));
confusion_matrix(colorcell_ancien, colorcell_classification)

% 5. À la fin de la boucle, juste après les deux « save ... colorcell_Maxwellnew » :
%    tout ce qui a servi à choisir les couleurs de cette session
details = struct();
details.date                 = datestr(now, 'yyyy-mm-dd HH:MM:SS');
details.session              = char(identifier);
details.colorcell            = colorcell;                 % couleur finale (après tes corrections)
details.colorcell_automatique = colorcell_classification; % avant les corrections à la main
details.cellules_corrigees   = find(colorcell(:) ~= colorcell_classification(:));
details.id_low_intensity     = id_low_intensity;
details.roi_brut             = roi_original;              % [rouge vert bleu] avant toute correction
details.roi_corrige          = roi_original_corrected;    % après fuite et seuils
details.roi_calcium          = roi_ca;
details.proportions_rgb      = [r, g, b];                 % proportions de Maxwell
details.triangle             = Donnees_Clustering;
details.maxwell              = info_maxwell;              % fond, échelle de chaque canal, intensité totale
details.marge_couleur        = marge_couleur;             % distance à la frontière la plus proche
details.fond_local           = extra_roi;                 % anneau autour de chaque cellule
details.fuite                = params_fuite;              % k et b proposés et validés
details.alignement           = info_alignement;           % décalages et corrélations
details.seuil_initial        = seuil_initial;
details.ep                   = ep;
save(strcat(namefull, 'couleurs_details.mat'), '-struct', 'details');
save(strcat(pathload, 'couleurs_details.mat'), '-struct', 'details');
fprintf('Recalage : rouge 890 / rouge 1040 r = %.2f (avant %.2f) ; calcium / rouge r = %.2f (avant %.2f) ; décalages %s et %s px\n', ...
    info_alignement.corr_rouges_apres, info_alignement.corr_rouges_avant, info_alignement.corr_calcium_apres, ...
    info_alignement.corr_calcium_avant, mat2str(round(info_alignement.decalage_1040_px, 1)), ...
    mat2str(round(info_alignement.decalage_890_px, 1)));
