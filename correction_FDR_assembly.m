%% --- SCRIPT FINAL : LECTURE ET CORRECTION FDR ---
clear all; close all; clc;

PathSave = 'C:\Votre\Chemin\'; % Remplacez par votre vrai chemin
path_mat_file = strcat(PathSave, 'analysis_assemblies_pvalues.mat');

% 1. Charger la table contenant les milliers d'assemblées
load(path_mat_file, 'Table_All_Assemblies');
disp(['Fichier chargé. Nombre total d''assemblées : ' num2str(height(Table_All_Assemblies))]);

% 2. Extraire uniquement les colonnes des p-values sous forme de matrice (double)
% (On suppose que les p-values sont dans les colonnes 4 à 9)
p_values_matrice = table2array(Table_All_Assemblies(:, 4:9));

% 3. Aplatir la matrice en un seul grand vecteur (pour appliquer la FDR sur TOUT)
p_values_vector = p_values_matrice(:);

% 4. Appliquer la correction de Benjamini-Hochberg (FDR)
% Si vous avez la Bioinformatics toolbox, décommentez la ligne mafdr :
% q_values_vector = mafdr(p_values_vector, 'BHFDR', true);

% SINON, voici l'algorithme exact de Benjamini-Hochberg :
m = length(p_values_vector); % Nombre total de tests
[p_sorted, sort_idx] = sort(p_values_vector); % Trier les p-values
q_sorted = p_sorted .* m ./ (1:m)'; % Calcul q-value = p * m / rang

% On s'assure que la suite des q-values est croissante (règle BH)
for i = m-1:-1:1
    q_sorted(i) = min(q_sorted(i), q_sorted(i+1));
end

% On remet les q-values dans l'ordre d'origine
q_values_vector = zeros(m, 1);
q_values_vector(sort_idx) = q_sorted;
% (les q-values sont plafonnées à 1 maximum)
q_values_vector(q_values_vector > 1) = 1; 

% 5. Remodeler le vecteur en matrice (même forme que p_values_matrice)
q_values_matrice = reshape(q_values_vector, size(p_values_matrice));

% 6. Intégrer ces q-values dans la table et trouver la couleur gagnante
seuil_FDR = 0.05;
nom_couleur = {'Red', 'Green', 'Blue', 'RG', 'RB', 'GB'};

% Ajouter les colonnes Q-values à la table existante
Table_All_Assemblies.Q_val_Rouge = q_values_matrice(:, 1);
Table_All_Assemblies.Q_val_Vert  = q_values_matrice(:, 2);
Table_All_Assemblies.Q_val_Bleu  = q_values_matrice(:, 3);
Table_All_Assemblies.Q_val_RG    = q_values_matrice(:, 4);
Table_All_Assemblies.Q_val_RB    = q_values_matrice(:, 5);
Table_All_Assemblies.Q_val_GB    = q_values_matrice(:, 6);

% Déterminer le résultat final (la couleur significative après FDR)
Resultat_Final = cell(height(Table_All_Assemblies), 1);

for i = 1:height(Table_All_Assemblies)
    % Trouver la meilleure q-value pour cette assemblée
    [min_q, idx_q] = min(q_values_matrice(i, :));
    
    if min_q < seuil_FDR
        Resultat_Final{i} = nom_couleur{idx_q};
    else
        Resultat_Final{i} = 'Not significant';
    end
end

% Ajouter le résultat à la table
Table_All_Assemblies.Couleur_Significative_FDR = Resultat_Final;

% 7. Statistiques rapides
nb_signif = sum(~strcmp(Resultat_Final, 'Not significant'));
pourcentage = (nb_signif / height(Table_All_Assemblies)) * 100;
disp(['Nombre d''assemblées significatives après FDR : ' num2str(nb_signif) ' (' num2str(pourcentage, '%.1f') '%)']);

% 8. Sauvegarder la table finale dans Excel pour la lire facilement
path_excel_final = strcat(PathSave, 'FINAL_assemblies_FDR_corrected.xlsx');
writetable(Table_All_Assemblies, path_excel_final);
disp(['Tableau final sauvegardé dans : ' path_excel_final]);