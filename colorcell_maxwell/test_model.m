function [certitude_max, scores, accuracy_pct,accuracy_pct_sans_noires, C_matrix,truth_labels] = test_model(MdlFinal, X, Ycat,cell_IDs,tableau,seuil_faible)

fprintf('\n\n====== ANALYSE DES PROBABILITÉS ET VÉRIFICATION GROUND TRUTH ======\n');

class_order = {'rouge','vert','bleu','jaune','magenta','cyan','blanc','noir'};
[y_pred, scores, ~] = predict(MdlFinal, X);
[certitude_max, ~] = max(scores, [], 2);

if iscategorical(Ycat)
    truth_labels = cellstr(Ycat);
elseif ischar(Ycat)
    truth_labels = cellstr(Ycat);
else
    truth_labels = Ycat;
end

pred_labels = cellstr(y_pred);

% --- 1. ACCURACY GLOBALE (Toutes les cellules) ---
nb_total = length(truth_labels);
nb_correct = sum(strcmp(truth_labels, pred_labels));
accuracy_pct = (nb_correct / nb_total) * 100;

% --- 2. ACCURACY SANS LES NOIRES ---
% On identifie les indices où la Ground Truth n'est pas "noir"
idx_sans_noires = find(~strcmp(truth_labels, 'noir'));
nb_total_sans_noires = length(idx_sans_noires);

if nb_total_sans_noires > 0
    nb_correct_sans_noires = sum(strcmp(truth_labels(idx_sans_noires), pred_labels(idx_sans_noires)));
    accuracy_pct_sans_noires = (nb_correct_sans_noires / nb_total_sans_noires) * 100;
else
    accuracy_pct_sans_noires = NaN; % Sécurité si jamais il n'y a QUE des noires
end

fprintf('Précision globale (Accuracy) : %.2f %% (%d / %d cellules correctes)\n', ...
    accuracy_pct, nb_correct, nb_total);
fprintf('Précision sur les couleurs (sans les noires) : %.2f %% (%d / %d cellules correctes)\n\n', ...
    accuracy_pct_sans_noires, nb_correct_sans_noires, nb_total_sans_noires);

% --- 3. MATRICE DE CONFUSION ---
[C_matrix, ~] = confusionmat(truth_labels, pred_labels, 'Order', class_order);

figure('Name', 'Matrice de Consensus', 'NumberTitle', 'off');
cm = confusionchart(C_matrix, class_order);
% On affiche les deux accuracy dans le titre
cm.Title = sprintf('Matrice de Confusion\nGlobale: %.1f%% | Sans noires: %.1f%%', ...
    accuracy_pct, accuracy_pct_sans_noires);
cm.RowSummary = 'row-normalized';
cm.ColumnSummary = 'column-normalized';
sortClasses(cm, class_order);

% =========================================================================
% 1. Cellules divergentes
% =========================================================================
erreurs_idx = find(~strcmp(truth_labels, pred_labels));

fprintf('L''algorithme est en désaccord avec la Ground Truth sur %d cellules.\n\n', length(erreurs_idx));

% CORRECTION : J'ai supprimé la ligne qui réinitialisait "cell_IDs" ici !

if isrow(cell_IDs)
    cell_IDs = cell_IDs';
end

divergent_table = table( ...
    cell_IDs(erreurs_idx), ...
    truth_labels(erreurs_idx), ...
    pred_labels(erreurs_idx), ...
    round(certitude_max(erreurs_idx) * 100, 1), ...
    'VariableNames', {'Cellule','GroundTruth','Prediction','CertitudePct'});

% =========================================================================
% 2. Cellules à faible probabilité
% =========================================================================
% seuil_faible = 0.60;
low_conf_idx = find(certitude_max < seuil_faible);

lowconf_table = table( ...
    cell_IDs(low_conf_idx), ...
    truth_labels(low_conf_idx), ...
    pred_labels(low_conf_idx), ...
    round(certitude_max(low_conf_idx) * 100, 1), ...
    'VariableNames', {'Cellule','GroundTruth','Prediction','CertitudePct'});

% =========================================================================
% 3. Affichage popup window
% =========================================================================
if tableau==1
    fig = uifigure('Name', 'Résultats de classification', 'Position', [100 100 900 700]);

uilabel(fig, ...
    'Text', sprintf('Cellules divergentes : %d', height(divergent_table)), ...
    'Position', [20 660 300 22], ...
    'FontWeight', 'bold', ...
    'FontSize', 14);

uitable(fig, ...
    'Data', divergent_table, ...
    'Position', [20 380 860 270], ...
    'ColumnWidth', {'auto','auto','auto','auto'});

uilabel(fig, ...
    'Text', sprintf('Cellules à faible probabilité (< %.0f%%) : %d', seuil_faible*100, height(lowconf_table)), ...
    'Position', [20 340 420 22], ...
    'FontWeight', 'bold', ...
    'FontSize', 14);

uitable(fig, ...
    'Data', lowconf_table, ...
    'Position', [20 40 860 290], ...
    'ColumnWidth', {'auto','auto','auto','auto'});
end
% =========================================================================
% 4. Console : erreurs très suspectes
% =========================================================================
seuil_suspect = 0.80;
suspects_idx_relatif = erreurs_idx(certitude_max(erreurs_idx) > seuil_suspect);

if ~isempty(suspects_idx_relatif)
    fprintf('\nATTENTION : cellules où l''algorithme est TRÈS SÛR (>80%%) mais contredit la GT.\n');
    disp(cell_IDs(suspects_idx_relatif)');
else
    fprintf('\nAucune erreur de classification avec une certitude extrême (>80%%).\n');
end

% =========================================================================
% 5. RECALCUL APRÈS REJET DES INCERTAINES
% =========================================================================
fprintf('\n====== ACCURACY APRÈS REJET DES CELLULES INCERTAINES (< %.0f%%) ======\n', seuil_faible*100);

% On crée un masque pour garder uniquement les prédictions sûres
idx_valides = certitude_max >= seuil_faible;
nb_rejetes = sum(~idx_valides);

% On filtre nos tableaux
truth_labels_filt = truth_labels(idx_valides);
pred_labels_filt = pred_labels(idx_valides);

if ~isempty(truth_labels_filt)
    % Recalcul global filtré
    nb_total_filt = length(truth_labels_filt);
    nb_correct_filt = sum(strcmp(truth_labels_filt, pred_labels_filt));
    accuracy_pct_filt = (nb_correct_filt / nb_total_filt) * 100;
    
    % Recalcul sans noires filtré
    idx_sans_noires_filt = find(~strcmp(truth_labels_filt, 'noir'));
    nb_total_sans_noires_filt = length(idx_sans_noires_filt);
    
    if nb_total_sans_noires_filt > 0
        nb_correct_sans_noires_filt = sum(strcmp(truth_labels_filt(idx_sans_noires_filt), pred_labels_filt(idx_sans_noires_filt)));
        accuracy_pct_sans_noires_filt = (nb_correct_sans_noires_filt / nb_total_sans_noires_filt) * 100;
    else
        accuracy_pct_sans_noires_filt = NaN;
    end
    
    fprintf('Cellules conservées : %d (%d rejetées sur %d)\n', nb_total_filt, nb_rejetes, nb_total);
    fprintf('Précision globale filtrée : %.2f %% (%d / %d)\n', accuracy_pct_filt, nb_correct_filt, nb_total_filt);
    fprintf('Précision sur couleurs filtrées : %.2f %% (%d / %d)\n\n', accuracy_pct_sans_noires_filt, nb_correct_sans_noires_filt, nb_total_sans_noires_filt);
    
    % Nouvelle matrice de confusion filtrée
    [C_matrix_filt, ~] = confusionmat(truth_labels_filt, pred_labels_filt, 'Order', class_order);
    figure('Name', 'Matrice de Consensus (Filtrée)', 'NumberTitle', 'off');
    cm_filt = confusionchart(C_matrix_filt, class_order);
    cm_filt.Title = sprintf('Matrice de Confusion (Certitude >= %.0f%%)\nGlobale: %.1f%% | Sans noires: %.1f%%', ...
        seuil_faible*100, accuracy_pct_filt, accuracy_pct_sans_noires_filt);
    cm_filt.RowSummary = 'row-normalized';
    cm_filt.ColumnSummary = 'column-normalized';
    sortClasses(cm_filt, class_order);
else
    fprintf('Aucune cellule n''a atteint le seuil de certitude.\n\n');
end

end