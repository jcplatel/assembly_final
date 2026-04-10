% 1. Charger le fichier Excel
% 'VariableNamingRule', 'preserve' permet de garder les espaces dans les noms de colonnes
T = readtable("E:\Data\Aurelie\analysis\March2026\nocues\test_stability\analysisALL.xlsx", 'VariableNamingRule', 'preserve');

% Nom exact de la colonne pour le recall dans votre fichier
col_recall = 'SClok'; 
% 2. Définir les poids (Ajustés selon votre proposition pour un total de 1.00)
w_sil = 0.30;         % sil all clusters (séparation mathématique)
w_recall = 0.20;      % recall (SClok)
w_cell = 0.20;        % cell_Purity
w_sce = 0.20;         % SCE_Purity
w_assemblies = 0.10;  % Assemblies (quantité)

% 3. Normalisation Globale Min-Max (entre 0 et 1)
% La fonction normalize ramène chaque métrique sur la même échelle pour les comparer
norm_sil = normalize(T.("sil all clusters"), 'range');
norm_recall = normalize(T.(col_recall), 'range');
norm_cell = normalize(T.cell_Purity, 'range');
norm_sce = normalize(T.SCE_Purity, 'range');
norm_assemblies = normalize(T.Assemblies, 'range');

% Remplacer les éventuels NaN générés par la normalisation (si variance nulle) par 0
norm_sil(isnan(norm_sil)) = 0;
norm_recall(isnan(norm_recall)) = 0;
norm_cell(isnan(norm_cell)) = 0;
norm_sce(isnan(norm_sce)) = 0;
norm_assemblies(isnan(norm_assemblies)) = 0;

% 4. Calcul du score global pour toutes les lignes
T.ScoreGlobal = (w_sil .* norm_sil) + ...
                (w_recall .* norm_recall) + ...
                (w_cell .* norm_cell) + ...
                (w_sce .* norm_sce) + ...
                (w_assemblies .* norm_assemblies);

% 5. Boucle d'analyse par IDENTIFIER unique
[unique_ids, ~, idx_id] = unique(T.identifier);

% Préparation d'une table pour stocker uniquement les meilleurs résultats
BestAnalyses = table();

% Création d'une figure unique qui sera mise à jour à chaque itération
fig = figure('Name', 'Analyse par Identifier', 'Color', 'w');

for i = 1:length(unique_ids)
    current_id = unique_ids(i);
    
    % Extraire les données spécifiques à cet identifier
    T_id = T(idx_id == i, :);
    
    % Trier les données par nombre de K (NCl_ini) pour avoir une belle courbe
    T_id = sortrows(T_id, 'NCl_ini');
    
    % Trouver le meilleur score et le meilleur K pour cet identifier
    [best_score, max_idx] = max(T_id.ScoreGlobal);
    best_K = T_id.NCl_ini(max_idx);
    
    % Récupérer le nombre final d'assemblées pour ce meilleur score
    best_assemblies = T_id.Assemblies(max_idx);
    
    % Ajouter la meilleure ligne au tableau final des résultats
    BestAnalyses = [BestAnalyses; T_id(max_idx, :)];
    
    % --- Visualisation ---
    clf(fig); % Nettoyer la figure pour le nouveau tracé
    
    % Tracer l'évolution du score en fonction de K
    plot(T_id.NCl_ini, T_id.ScoreGlobal, '-o', 'LineWidth', 2, ...
        'MarkerSize', 6, 'MarkerFaceColor', '#0072BD', 'Color', '#0072BD');
    hold on;
    
    % Marquer le meilleur K avec une étoile rouge
    plot(best_K, best_score, 'p', 'MarkerSize', 14, ...
        'MarkerFaceColor', '#D95319', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    
    % Mise en forme du graphique
    grid on;
    xlabel('Nombre de clusters (NCl\_ini)', 'FontWeight', 'bold');
    ylabel('Score Global Pondéré', 'FontWeight', 'bold');
    
    % Intégration des informations dans le titre et le sous-titre
    title(sprintf('Identifier : %s', string(current_id)), 'Interpreter', 'none');
    subtitle(sprintf('Meilleur K = %d (Score = %.3f) | Nombre d''assemblées = %d', ...
        best_K, best_score, best_assemblies), 'FontSize', 11, 'Color', '#D95319');
    
    legend('Score Global', 'Meilleur K', 'Location', 'best');
    
    % Pause pour permettre l'inspection visuelle
    disp(['Affichage de l''identifier : ', char(string(current_id)), ' (', num2str(i), '/', num2str(length(unique_ids)), ')']);
    disp('Appuyez sur une touche pour passer au suivant... (ou Ctrl+C pour arrêter)');
    pause; 
end

% 6. Sauvegarde des meilleurs résultats
writetable(BestAnalyses, 'Meilleures_Analyses_Par_Identifier.xlsx');
disp('Terminé ! Les meilleures analyses de chaque identifier ont été sauvegardées.');
close(fig);