function bilan = refaire_images_recalees(dossiers)
%% =========================================================================
%  REFAIRE_IMAGES_RECALEES (17 nov., demandée par l'utilisateur)
%  Étapes 1 et 2 de color_predictor_model, sans fenêtre : lecture des images
%  brutes (red, red890, green, blue, calcium.tif), retrait du fond et
%  normalisation (brainbow_image_processing), recalage sur le GCaMP
%  (image_alignment_multimodal). Pas de seuil ni de correction de fuite.
%
%  Écrit dans chaque dossier de session (rien n'est écrasé : dossier sauté
%  si aligned_red_V2.tif existe déjà) :
%    aligned_red_V2.tif, aligned_green_V2.tif, aligned_blue_V2.tif,
%    calcium_norm_V2.tif           (16 bits, valeurs 0-1 x 65535)
%    compositecalcium_V2.tif       (4 pages : rouge, vert, bleu, GCaMP)
%    recalage_V2.mat               (décalages, corrélations, date)
%  Les anciens aligned_*.tif / compositecalcium.tif (recalage de 2023-2025)
%  ne sont pas touchés.
%
%  Lancer :  bilan = refaire_images_recalees({'E:\...\220919_plane0', ...})
%  =========================================================================

noms = {'aligned_red_V2.tif', 'aligned_green_V2.tif', 'aligned_blue_V2.tif', 'calcium_norm_V2.tif'};
n = numel(dossiers);
bilan = table(dossiers(:), repmat({''}, n, 1), nan(n, 2), nan(n, 2), nan(n, 1), nan(n, 1), ...
    'VariableNames', {'dossier', 'etat', 'decalage_1040_px', 'decalage_890_px', 'corr_calcium_apres', 'corr_rouges_apres'});
for i = 1:n
    d = regexprep(char(dossiers{i}), '[\\/]+$', '');
    try
        if isfile(fullfile(d, noms{1}))
            bilan.etat{i} = 'déjà fait (sauté)';
            continue
        end
        [red, red890, green, blue, calcium] = brainbow_image_read([d filesep]);
        In = brainbow_image_processing([d filesep], red, red890, green, blue, calcium);
        [ar, ag, ab, cn, ~, info] = image_alignment_multimodal(In);
        [ar, ag, ab, info.ajustement_masques] = ajuster_canaux_masques(ar, ag, ab, d);   % 17 nov.
        I = {ar, ag, ab, cn};
        for k = 1:4
            imwrite(im2uint16(min(max(I{k}, 0), 1)), fullfile(d, noms{k}), 'tif');
        end
        fc = fullfile(d, 'compositecalcium_V2.tif');
        for k = 1:4
            if k == 1, mode = 'overwrite'; else, mode = 'append'; end
            imwrite(im2uint16(min(max(I{k}, 0), 1)), fc, 'tif', 'WriteMode', mode);
        end
        info.date = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<DATST,TNOW1>
        info.fonctions = 'brainbow_image_read, brainbow_image_processing, image_alignment_multimodal (color_maxwellV2)';
        save(fullfile(d, 'recalage_V2.mat'), '-struct', 'info');
        bilan.etat{i} = 'ok';
        bilan.decalage_1040_px(i, :) = info.decalage_1040_px;
        bilan.decalage_890_px(i, :) = info.decalage_890_px;
        bilan.corr_calcium_apres(i) = info.corr_calcium_apres;
        bilan.corr_rouges_apres(i) = info.corr_rouges_apres;
    catch err
        bilan.etat{i} = ['ERREUR : ' err.message];
    end
    fprintf('%d/%d %s : %s\n', i, n, d, bilan.etat{i});
end
end
