function [colorcell_new, certitude_max, T_predictions_new, T_low_conf_new] = predict_colors_LDA(MdlLDA, roi_red, roi_green, roi_blue, roi_ca, cell_IDs, seuil_faible)

    if nargin < 7 || isempty(seuil_faible)
        seuil_faible = 0.60;
    end

    classNames = {'red','green','blue','yellow','magenta','cyan','white','black'};
    NCell = length(roi_red);

    % =========================================================================
    % 1. IDENTIFICATION DES CELLULES NOIRES (NaN)
    % =========================================================================
    % Si un des canaux est NaN, la cellule est considérée comme noire/rejetée
    idx_nan = isnan(roi_red) | isnan(roi_green) | isnan(roi_blue);
    idx_valides = ~idx_nan;

    % Pré-allocation avec les valeurs par défaut pour les cellules noires
    pred_labels_new = repmat({'black'}, NCell, 1);
    colorcell_new = repmat(8, NCell, 1); % 8 = black
    certitude_max = nan(NCell, 1);       % NaN pour la certitude (car forcé manuellement)

    % =========================================================================
    % 2. PRÉPARATION ET PRÉDICTION (UNIQUEMENT SUR LES CELLULES VALIDES)
    % =========================================================================
    if any(idx_valides)
        % On n'extrait les variables que pour les cellules non-noires
        R_v = roi_red(idx_valides);
        G_v = roi_green(idx_valides);
        B_v = roi_blue(idx_valides);
        Ca_v = roi_ca(idx_valides);

        % Création des features
        I_v = sqrt(R_v.^2 + G_v.^2 + B_v.^2);
        ratio_RG = R_v ./ (G_v + 0.01);
        ratio_RB = R_v ./ (B_v + 0.01);
        ratio_GB = G_v ./ (B_v + 0.01);

        X_valid = [R_v, G_v, B_v, Ca_v, I_v, ratio_RG, ratio_RB, ratio_GB];

        % Prédiction par la LDA
        [y_pred, scores] = predict(MdlLDA, X_valid);
        [cert_max, ~] = max(scores, [], 2);
        
        % On réinjecte les prédictions valides dans le vecteur complet
        pred_labels_new(idx_valides) = cellstr(y_pred);
        certitude_max(idx_valides) = cert_max;

        % Conversion des labels en numéros (1 à 7)
        labels_v = pred_labels_new(idx_valides);
        color_v = zeros(length(labels_v), 1);
        for k = 1:length(labels_v)
            idx = find(strcmp(classNames, labels_v{k}));
            if ~isempty(idx)
                color_v(k) = idx;
            end
        end
        colorcell_new(idx_valides) = color_v;
    end

    % =========================================================================
    % 3. CRÉATION DES TABLES ET RÉSULTATS
    % =========================================================================
    if isrow(cell_IDs)
        cell_IDs = cell_IDs';
    end

    T_predictions_new = table( ...
        cell_IDs, ...
        pred_labels_new, ...
        colorcell_new, ...
        round(certitude_max * 100, 1), ...
        'VariableNames', {'Cellule', 'PredictionTexte', 'Colorcell', 'CertitudePct'});

    % On identifie les incertaines en excluant les noires (qui ont certitude_max = NaN)
    idx_low_conf = (certitude_max < seuil_faible) & idx_valides;
    T_low_conf_new = T_predictions_new(idx_low_conf, :);

    % =========================================================================
    % 4. AFFICHAGE CONSOLE
    % =========================================================================
    nb_noires = sum(idx_nan);
    nb_valides_ok = sum(certitude_max >= seuil_faible);
    nb_incertaines = sum(idx_low_conf);

    fprintf('\n====== RÉSULTATS DE LA PRÉDICTION HYBRIDE ======\n');
    fprintf('Cellules totales   : %d\n', NCell);
    fprintf('  - Noires (forcées): %d\n', nb_noires);
    fprintf('  - Sûres (>= %.0f%%) : %d\n', seuil_faible*100, nb_valides_ok);
    fprintf('  - Incertaines     : %d\n\n', nb_incertaines);
end