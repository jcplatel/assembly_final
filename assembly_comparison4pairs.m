% =========================================================================
% SCRIPT 1 : Comparaison d'analyses pour chaque K (Similarité & Clusters finaux)
% =========================================================================

% Définir les dossiers à comparer (À MODIFIER avec vos chemins)
dossiers = {
    "E:\Data\Aurelie\analysis\March2026\nocues\test_stability\411582_230320_plane0_26_04_09_17_34_56", ...
    "E:\Data\Aurelie\analysis\March2026\nocues\test_stability\411582_230320_plane0_26_04_09_17_46_27", ...
    "E:\Data\Aurelie\analysis\March2026\nocues\test_stability\411582_230320_plane0_26_04_09_17_57_17", ...
};
Ks = 3:20;
nb_analyses = length(dossiers);
nb_paires = nchoosek(nb_analyses, 2); 

% Matrices pour stocker les résultats (remplies de NaN pour éviter les chutes à 0 si un K manque)
similarites_paires_K = NaN(length(Ks), nb_paires);
similarite_globale_K = NaN(length(Ks), 1);
nb_clusters_finaux_K = NaN(length(Ks), nb_analyses); 

for idx_k = 1:length(Ks)
    K = Ks(idx_k);
    assemblies = cell(nb_analyses, 1);
    data_missing = false;
    
    % Charger dynamiquement les données pour les analyses
    for d = 1:nb_analyses
        % On liste tous les éléments du dossier et on ne garde que les sous-dossiers
        tous = dir(dossiers{d});
        tous = tous([tous.isdir]);
        tous = tous(~ismember({tous.name}, {'.', '..'})); % Ignorer les dossiers système
        
        idx_match = [];
        % RECHERCHE CIBLÉE :
        % Cherche le format exact "K3" (ou "k3"). L'expression (?!\d) empêche
        % de confondre "K3" avec "K30" par exemple.
        motif_exact = ['K' num2str(K) '(?!\d)'];
        
        for idx_dossier = 1:length(tous)
            nom = tous(idx_dossier).name;
            if ~isempty(regexpi(nom, motif_exact, 'once'))
                idx_match = idx_dossier;
                break;
            end
        end
        
        if isempty(idx_match)
            data_missing = true; break;
        end
        
        dossier_K = tous(idx_match);
        fichier = dir(fullfile(dossiers{d}, dossier_K.name, '*.mat'));
        
        if isempty(fichier)
            data_missing = true; break;
        end
        
        data = load(fullfile(fichier(1).folder, fichier(1).name), 'assemblyortho');
        assemblies{d} = data.assemblyortho;
        
        % --- Sauvegarde du nombre de clusters finaux ---
        nb_clusters_finaux_K(idx_k, d) = length(data.assemblyortho);
    end
    
    % Si une analyse manque pour ce K, on laisse des NaN et on passe au K suivant
    if data_missing
        continue;
    end
    
    % Calculer la similarité de Jaccard pour chaque paire d'analyses
    idx_paire = 1;
    for i = 1:nb_analyses-1
        for j = (i+1):nb_analyses
            assembly_1 = assemblies{i};
            assembly_2 = assemblies{j};
            
            matrice_jaccard = zeros(length(assembly_1), length(assembly_2));
            
            for c1 = 1:length(assembly_1)
                for c2 = 1:length(assembly_2)
                    set_1 = assembly_1{c1};
                    set_2 = assembly_2{c2};
                    intersection = length(intersect(set_1, set_2));
                    union_12 = length(union(set_1, set_2));
                    matrice_jaccard(c1, c2) = intersection / max(union_12, 1);
                end
            end
            
            % Pour chaque cluster de l'analyse 1, on prend son meilleur match dans la 2
            meilleurs_matchs = max(matrice_jaccard, [], 2);
            
            % Sauvegarder la similarité moyenne de cette paire spécifique
            similarites_paires_K(idx_k, idx_paire) = mean(meilleurs_matchs);
            idx_paire = idx_paire + 1;
        end
    end
    
    % Sauvegarder la moyenne globale sur toutes les paires pour ce K ('omitnan' évite les bugs)
    similarite_globale_K(idx_k) = mean(similarites_paires_K(idx_k, :), 'omitnan');
end

% =========================================================================
% Affichage des courbes (Similarité)
% =========================================================================
figure('Name', 'Accord entre Analyses selon K', 'Color', 'w');
hold on;
% 1. Tracer les courbes individuelles des paires (en gris clair)
plot_paires = plot(Ks, similarites_paires_K, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
% 2. Tracer la moyenne globale de toutes les comparaisons (en bleu fort)
plot_moyenne = plot(Ks, similarite_globale_K, '-o', 'LineWidth', 2.5, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'b', 'Color', 'b');
xlabel('Nombre de clusters (K)', 'FontWeight', 'bold');
ylabel('Similarité Jaccard moyenne', 'FontWeight', 'bold');
title('Stabilité des clusters (comparaison croisée sur les analyses)');
grid on;
legend([plot_paires(1), plot_moyenne], {sprintf('Paires individuelles (%d)', nb_paires), 'Moyenne globale'}, ...
    'Location', 'best');
ylim ([0.50 1]);
hold off;

% =========================================================================
% Affichage des courbes (Nombre de clusters finaux)
% =========================================================================
figure('Name', 'Nombre de clusters finaux selon K', 'Color', 'w');
hold on;

% Tracer le nombre de clusters finaux pour chaque analyse
plot_clusters = plot(Ks, nb_clusters_finaux_K, '-x', 'LineWidth', 1.5);

% Tracer une ligne pointillée noire représentant y=x (K final = K initial)
plot_ref = plot(Ks, Ks, 'k--', 'LineWidth', 1);

xlabel('Nombre de clusters initial (K)', 'FontWeight', 'bold');
ylabel('Nombre de clusters finaux obtenus', 'FontWeight', 'bold');
title('Nombre de clusters finaux selon le K de départ');
grid on;

% Création dynamique de la légende
legend_labels = cell(nb_analyses, 1);
for d = 1:nb_analyses
    legend_labels{d} = sprintf('Analyse %d', d);
end
legend([plot_clusters; plot_ref], [legend_labels; {'Référence (Final = Initial)'}], 'Location', 'best');

hold off;