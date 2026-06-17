function confusion_matrix(colorcell_method1, colorcell_method2, Donnees_Clustering, R, G, B, write_number)


labels_str = {'Red', 'Green', 'Blue', 'Yellow', 'Magenta', 'Cyan', 'White', 'Black'};
n_classes = length(labels_str);
all_classes = 1:8;

% --- 3. Calcul de la Matrice de Confusion brute ---
% C(i,j) = nombre de neurones classés 'i' par méthode 1 et 'j' par méthode 2
C = confusionmat(colorcell_method1, colorcell_method2, 'Order', all_classes);

% --- 4. Calcul de l'Accuracy sur neurones actifs (sans Black) ---
idx_black = 8;
idxs_to_keep = setdiff(1:8, idx_black);
C_active = C(idxs_to_keep, idxs_to_keep);

% Calcul de la nouvelle accuracy (uniquement sur classes 1 à 7)
accuracy_active = trace(C_active) / sum(C_active(:)) * 100;
fprintf('Accuracy sur les neurones actifs (sans Black) : %.2f%%\n', accuracy_active);

% --- 5. Création et manipulation du ConfusionChart ---
figure;

% On affiche le chart d'abord de manière brute pour garder toutes les classes visibles
cm = confusionchart(C, labels_str);
sortClasses(cm, labels_str);

% Pour exclure les noirs des pourcentages marginaux (RowSummary et ColumnSummary)
% On désactive les summary automatiques (qui prendraient les noirs en compte)
cm.RowSummary = 'row-normalized';
cm.ColumnSummary = 'column-normalized';

% On remplace les données affichées dans les cases par une matrice personnalisée
% où les pourcentages sont calculés par rapport à C_active.
% On garde les valeurs brutes dans le tableau principal, mais on ajoute
% manuellement les totaux marginaux normalisés sans les noirs.

% Note : MATLAB ne permet pas de modifier la façon dont RowSummary et ColumnSummary
% calculent leurs totaux internes en gardant la ligne affichée. 
% L'astuce consiste à fournir au graphique une matrice où la classe noire 
% est "désactivée" dans les totaux.

% Nous appliquons la normalisation manuellement en désactivant la normalisation interne
cm.Normalization = 'absolute'; 

cm.Title = sprintf('Comparaison des Méthodes de Classification\nAccuracy (hors noirs) : %d%%', ceil(accuracy_active));
cm.XLabel = 'Biapy';
cm.YLabel = 'in vivo';

% Optionnel : Si vous voulez VRAIMENT voir les pourcentages en bout de ligne 
% basés uniquement sur les actifs dans le tableau (ce qui écrase les fonctions natives),
% vous pouvez l'afficher dans la console :

fprintf('\n--- Répartition des prédictions (excluant les totaux noirs) ---\n');
for i = 1:7
    total_ligne_actif = sum(C(i, 1:7)); % Somme de la ligne i, colonnes 1 à 7
    total_colonne_actif = sum(C(1:7, i)); % Somme de la colonne i, lignes 1 à 7
    
    if total_ligne_actif > 0
        rappel = C(i,i) / total_ligne_actif * 100;
    else
        rappel = 0;
    end
    
    if total_colonne_actif > 0
        precision = C(i,i) / total_colonne_actif * 100;
    else
        precision = 0;
    end
    
    fprintf('Classe %s : Rappel (vrai %s) = %.1f%% | Précision (prédit %s) = %.1f%%\n', ...
        labels_str{i}, labels_str{i}, rappel, labels_str{i}, precision);
end

end