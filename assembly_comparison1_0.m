% =========================================================================
% SCRIPT 1 : Comparaison A vs B pour chaque K (Courbe de similarité)
% =========================================================================
dossier_A = "E:\Data\Aurelie\analysis\March2026\nocues\test_stability\411582_230320_plane0_26_04_07_19_03_26"; % À MODIFIER
dossier_B = "E:\Data\Aurelie\analysis\March2026\nocues\test_stability\411582_230320_plane0_26_04_07_19_17_15"; % À MODIFIER
% "E:\Data\Aurelie\analysis\March2026\nocues\test_stability\411582_230320_plane0_26_04_07_19_31_04"
Ks = 4:20;
similarite_moyenne_K = zeros(length(Ks), 1);

for idx_k = 1:length(Ks)
    K = Ks(idx_k);
    
    % Trouver le fichier pour l'Analyse 1
    dossier_KA = dir(fullfile(dossier_A, ['*' num2str(K) '*']));
    if isempty(dossier_KA), continue; end
    fichier_A = dir(fullfile(dossier_A, dossier_KA(1).name, '*.mat'));
    if isempty(fichier_A), continue; end
    data_A = load(fullfile(fichier_A(1).folder, fichier_A(1).name), 'assemblyortho');
    assembly_A = data_A.assemblyortho;
    
    % Trouver le fichier pour l'Analyse 2
    dossier_KB = dir(fullfile(dossier_B, ['*' num2str(K) '*']));
    if isempty(dossier_KB), continue; end
    fichier_B = dir(fullfile(dossier_B, dossier_KB(1).name, '*.mat'));
    if isempty(fichier_B), continue; end
    data_B = load(fullfile(fichier_B(1).folder, fichier_B(1).name), 'assemblyortho');
    assembly_B = data_B.assemblyortho;
    
    % Matrice de Jaccard pour ce K entre Analyse A et Analyse B
    % (Pénalise la différence de taille des clusters) [web:25]
    matrice_jaccard = zeros(length(assembly_A), length(assembly_B));
    
    for i = 1:length(assembly_A)
        for j = 1:length(assembly_B)
            set_A = assembly_A{i};
            set_B = assembly_B{j};
            intersection = length(intersect(set_A, set_B));
            union_AB = length(union(set_A, set_B));
            matrice_jaccard(i, j) = intersection / max(union_AB, 1);
        end
    end
    
    % Pour chaque cluster de A, on prend son meilleur match dans B
    meilleurs_matchs = max(matrice_jaccard, [], 2);
    
    % On sauvegarde la similarité moyenne pour ce K
    similarite_moyenne_K(idx_k) = mean(meilleurs_matchs);
end

% Affichage de la courbe d'accord entre les deux analyses
figure('Name', 'Accord entre Analyse 1 et 2 selon K', 'Color', 'w');
plot(Ks, similarite_moyenne_K, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
xlabel('Nombre de clusters (K)');
ylabel('Similarité Jaccard moyenne (A vs B)');
title('Stabilité des clusters entre les deux analyses');
grid on;