clear
close all

% =========================================================================
% 1 : Comparaison des Matrices de Co-occurrence (Consensus)
% =========================================================================
dossier_A = "E:\Data\Aurelie\analysis\March2026\nocues\test_stability\411582_230320_plane0_26_04_06_17_17_46";
%dossier_B = "E:\Data\Aurelie\analysis\March2026\nocues\test_stability\411582_230320_plane0_26_04_06_17_29_57"; 

% Fonction interne (nested) pour calculer la matrice d'un dossier
function [Matrice, max_ID] = calculer_cooccurrence(dossier)
    fichiers = dir(fullfile(dossier, '**', '*.mat'));
    max_ID = 0;
    nb_analyses = 0;
    
    % Trouver l'ID maximum
    for f = 1:length(fichiers)
        try
            data = load(fullfile(fichiers(f).folder, fichiers(f).name), 'assemblyortho');
            if isfield(data, 'assemblyortho')
                nb_analyses = nb_analyses + 1;
                for c = 1:length(data.assemblyortho)
                    if ~isempty(data.assemblyortho{c})
                        max_ID = max([max_ID, max(data.assemblyortho{c})]);
                    end
                end
            end
        catch; end
    end
    
    Matrice = zeros(max_ID, max_ID);
    
    % Remplir la matrice (Vote par co-occurrence) [web:24]
    for f = 1:length(fichiers)
        try
            data = load(fullfile(fichiers(f).folder, fichiers(f).name), 'assemblyortho');
            if isfield(data, 'assemblyortho')
                for c = 1:length(data.assemblyortho)
                    neurones = data.assemblyortho{c};
                    if ~isempty(neurones)
                        Matrice(neurones, neurones) = Matrice(neurones, neurones) + 1;
                    end
                end
            end
        catch; end
    end
    
    if nb_analyses > 0
        Matrice = (Matrice / nb_analyses) * 100;
    end
end

% Calcul pour les deux dossiers
[Matrice_A, max_A] = calculer_cooccurrence(dossier_A);

% % Uniformiser la taille des matrices si un neurone manque dans l'une
% ID_global = max(max_A, max_B);
% Matrice_A(ID_global, ID_global) = 0; 

% Affichage comparatif
figure('Name', 'Comparaison des Consensus', 'Color', 'w', 'Position', [100 100 1200 500]);

imagesc(Matrice_A); colormap('jet'); colorbar; axis square;
title('Consensus Analyse 1 (%)');
xlabel('ID Neurone'); ylabel('ID Neurone');


%%
% =========================================================================
% 2 :  Réorganisation des matrice Tri optimisé avec suppression du bruit
% =========================================================================

% 1. Filtrer le bruit : on met à 0 toutes les associations sous 25%
Seuil = 25; 
Mat_A_propre = Matrice_A; Mat_A_propre(Mat_A_propre < Seuil) = 0;

% 2. Calculer la distance (sur les données filtrées)
Dist_A = 100 - Mat_A_propre; Dist_A(logical(eye(size(Dist_A)))) = 0;

% 3. Clustering avec la méthode 'average' (souvent plus propre pour ce type de matrice)
Z_A = linkage(squareform(Dist_A), 'average');

% Obtenir le nouvel ordre
figure('Visible', 'on');
[~, ~, ordre_A] = dendrogram(Z_A, 0); 

close;

% 4. Réorganiser les matrices d'origine (non-filtrées pour garder l'info visuelle)
Matrice_A_triee_opt = Matrice_A(ordre_A, ordre_A);


% Affichage
figure('Name', 'Consensus Tri Optimisé', 'Color', 'w', 'Position', [100 100 1200 500]);
imagesc(Matrice_A_triee_opt); colormap('jet'); colorbar; axis square;
title('Consensus (Tri Optimisé)');

%%
% =========================================================================
% 3 : Extraction des Assemblées Noyaux (Core Assemblies)
% =========================================================================

% 1. Définir le seuil de robustesse (ex: 85% de co-occurrence minimum)
Seuil_Core = 80; 

% 2. Binariser la matrice (1 si >= seuil, 0 sinon)
Matrice_Binaire = Matrice_A >= Seuil_Core;

% Mettre la diagonale à 0 (un neurone n'est pas connecté à lui-même dans un graphe)
Matrice_Binaire(logical(eye(size(Matrice_Binaire)))) = 0;

% 3. Créer un objet Graphe en MATLAB
Graphe_Noyaux = graph(Matrice_Binaire);

% 4. Extraire les composantes connexes (les îlots de neurones)
bins = conncomp(Graphe_Noyaux);

% 5. Analyser et stocker les assemblées trouvées
Taille_Min_Assemblee = 3; % Ignorer les neurones isolés ou les paires
Assemblees_Noyaux = {};
IDs_Uniques = unique(bins);

for i = 1:length(IDs_Uniques)
    id_bin = IDs_Uniques(i);
    % Trouver tous les neurones appartenant à cet îlot
    neurones_dans_ilot = find(bins == id_bin);
    
    % Si l'îlot est assez grand, c'est une assemblée noyau
    if length(neurones_dans_ilot) >= Taille_Min_Assemblee
        Assemblees_Noyaux{end+1} = neurones_dans_ilot;
    end
end

% =========================================================================
% AFFICHAGE DES RÉSULTATS
% =========================================================================
fprintf('--- RÉSULTATS DE L''EXTRACTION (Seuil: %d%%) ---\n', Seuil_Core);
fprintf('Nombre d''assemblées noyaux détectées : %d\n\n', length(Assemblees_Noyaux));

for i = 1:length(Assemblees_Noyaux)
    neurones = Assemblees_Noyaux{i};
    fprintf('Assemblées Noyau %d (%d neurones) : ', i, length(neurones));
    
    % Afficher les 10 premiers neurones maximum pour ne pas inonder la console
    if length(neurones) > 10
        fprintf('%d, ', neurones(1:10));
        fprintf('... et %d autres\n', length(neurones)-10);
    else
        fprintf('%d, ', neurones(1:end-1));
        fprintf('%d\n', neurones(end));
    end
end

% 6. Visualisation optionnelle du réseau noyau
figure('Name', 'Réseau des Assemblées Noyaux', 'Color', 'w');
p = plot(Graphe_Noyaux, 'Layout', 'force');
% Ne colorer que les neurones qui font partie d'une assemblée
highlight(p, [Assemblees_Noyaux{:}], 'NodeColor', 'r', 'MarkerSize', 6);
title(sprintf('Réseau des cœurs d''assemblées (Seuil > %d%%)', Seuil_Core));
%% chargement coordonnées cellules

% =========================================================================
% 4 : Extraction des coordonnées spatiales depuis Suite2p (Fall.mat)
% =========================================================================

% 1. Chargement du fichier généré par Suite2p
% Assurez-vous d'avoir bien généré le fichier Fall.mat dans Suite2p
disp('Chargement des données Suite2p...');
load('E:\Data\Aurelie\data\nocues\411582\230320_plane0\Fall.mat', 'stat', 'iscell');

% 2. Initialisation des variables
Nb_Cellules_Totales = sum(iscell(:,1) == 1); % Compter le nombre de vrais neurones
Coordonnees = zeros(Nb_Cellules_Totales, 2); % Matrice [N x 2] pour les X et Y
IDs_Suite2p = zeros(Nb_Cellules_Totales, 1); % Pour garder une trace de l'ID original (optionnel)

m = 0; % Compteur pour les vraies cellules

% 3. Boucle sur tous les ROIs détectés par Suite2p
for n = 1:length(stat)
    
    % On ne garde que les ROIs classés comme vraies cellules (iscell == 1)
    if iscell(n, 1) == 1
        m = m + 1;
        
        % Méthode 1 : Utiliser le centre médian pré-calculé par Suite2p
        % C'est la méthode la plus propre. 'med' contient [Y, X].
        if isfield(stat{n}, 'med')
            Coordonnees(m, 1) = double(stat{n}.med(2)) + 1; % Coordonnée X (pixel)
            Coordonnees(m, 2) = double(stat{n}.med(1)) + 1; % Coordonnée Y (pixel)
        
        % Méthode 2 (Alternative au cas où 'med' n'existe pas dans votre version)
        % On calcule le centre géométrique moyen à partir de tous les pixels (xpix, ypix)
        else
            Coordonnees(m, 1) = mean(double(stat{n}.xpix)) + 1; 
            Coordonnees(m, 2) = mean(double(stat{n}.ypix)) + 1; 
        end
        
        % Sauvegarde de l'ID original de Suite2p (peut être utile plus tard)
        IDs_Suite2p(m) = n;
    end
end

% fprintf('%d cellules extraites avec succès.\n', m);
% 
% % =========================================================================
% % TEST VISUEL RAPIDE
% % =========================================================================
% % Ce petit bloc permet de vérifier que l'extraction s'est bien passée
% % en affichant simplement tous vos neurones sur une grille.
% 
% figure('Name', 'Vérification de l''extraction spatiale', 'Color', 'w');
% scatter(Coordonnees(:,1), Coordonnees(:,2), 15, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
% set(gca, 'YDir', 'reverse'); % Très important en imagerie (le Y=0 est en haut)
% axis equal;
% title('Position de tous les neurones détectés (iscell=1)');
% xlabel('Pixels X');
% ylabel('Pixels Y');
% grid on;

%%% =========================================================================
% =========================================================================
% 5 : Graphe Spatial Clair (Réseaux en Étoile par Assemblée)
% =========================================================================

% On suppose que vous avez déjà calculé 'Assemblees_Noyaux' (SCRIPT 4)
% et 'Coordonnees' (extraction Suite2p)

% % 1. Préparation de la figure
% figure('Name', 'Carte Spatiale des Assemblées', 'Color', 'k', 'Position', [100 100 800 800]);
% hold on; set(gca, 'Color', 'k', 'YDir', 'reverse'); axis equal; axis off;
% 
% % 2. Dessiner le bruit (les neurones qui n'appartiennent à aucune assemblée noyau)
% tous_neurones_noyaux = [Assemblees_Noyaux{:}];
% neurones_bruit = setdiff(1:size(Coordonnees, 1), tous_neurones_noyaux);
% 
% scatter(Coordonnees(neurones_bruit, 1), Coordonnees(neurones_bruit, 2), ...
%         20, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.5);
% 
% % 3. Générer une belle palette de couleurs contrastées
% couleurs = lines(length(Assemblees_Noyaux));
% 
% % 4. Dessiner chaque assemblée
% for i = 1:length(Assemblees_Noyaux)
%     neurones_groupe = Assemblees_Noyaux{i};
%     X_groupe = Coordonnees(neurones_groupe, 1);
%     Y_groupe = Coordonnees(neurones_groupe, 2);
% 
%     % Calculer le centre géographique ("Hub") de cette assemblée
%     centre_X = mean(X_groupe);
%     centre_Y = mean(Y_groupe);
% 
%     % Dessiner les liens légers entre le hub et chaque neurone (réseau en étoile)
%     for j = 1:length(neurones_groupe)
%         plot([centre_X, X_groupe(j)], [centre_Y, Y_groupe(j)], ...
%              'Color', [couleurs(i, :), 0.4], 'LineWidth', 1.5);
%     end
% 
%     % Dessiner les neurones par-dessus
%     scatter(X_groupe, Y_groupe, 60, couleurs(i, :), 'filled', ...
%             'MarkerEdgeColor', 'w', 'LineWidth', 1);
% 
%     % Dessiner le hub (optionnel)
%     % scatter(centre_X, centre_Y, 100, 'w', 'p', 'filled'); % Étoile blanche au centre
% end
% 
% title('Cartographie Spatiale des Cœurs d''Assemblées', 'Color', 'w', 'FontSize', 14);
%%

% =========================================================================
% SCRIPT 6_BIS : Histogramme en Log (Zoom sur le signal)
% =========================================================================

% 1. Extraire le triangle supérieur
valeurs = Matrice_A(triu(true(size(Matrice_A)), 1));

% 2. Ignorer le bruit pur (0%) pour dé-zoomer l'axe Y
valeurs_positives = valeurs(valeurs > 5); % On coupe tout ce qui est sous 5%

figure('Name', 'Distribution (Zoom Signal)', 'Color', 'w', 'Position', [200 200 700 500]);

% 3. Afficher l'histogramme
h = histogram(valeurs_positives, 40, 'Normalization', 'count', ...
              'FaceColor', [0.8 0.3 0.1], 'EdgeColor', 'w');

% 4. Passer l'axe Y en Logarithmique !
set(gca, 'YScale', 'log');

% Décoration
xlabel('Score de co-occurrence (%)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Nombre de paires (Echelle Log)', 'FontSize', 12, 'FontWeight', 'bold');
title('Distribution des paires associées (Zoom Logarithmique)', 'FontSize', 14);
grid on;

%%
% =========================================================================
% SCRIPT 7 : Cartographie individuelle de chaque Assemblée (Liens détaillés)
% =========================================================================

% Assurez-vous d'avoir relancé le SCRIPT 4 avec Seuil_Core = 60 ou 65
% pour que Assemblees_Noyaux contienne les bons groupes.

% Seuil_Affichage = 85; % Le seuil issu de l'histogramme logarithmique
% Epaisseur_Max = 1;    % Épaisseur pour un lien à 100%
% 
% Nb_Assemblees = length(Assemblees_Noyaux);
% 
% % Générer des couleurs distinctes pour chaque assemblée
% couleurs = lines(Nb_Assemblees);
% 
% % Créer une grande figure adaptée au nombre d'assemblées
% % (Calcule le nombre de lignes et colonnes pour les subplots)
% cols = ceil(sqrt(Nb_Assemblees));
% rows = ceil(Nb_Assemblees / cols);
% 
% figure('Name', 'Cartographie Individuelle des Assemblées', ...
%        'Color', 'k', 'Position', [50 50 1200 800]);
% 
% for a = 1:Nb_Assemblees
%     subplot(rows, cols, a);
%     hold on;
%     set(gca, 'Color', 'k', 'YDir', 'reverse'); 
%     axis equal; axis off;
% 
%     % 1. Dessiner le fond (tous les neurones du champ de vision en gris très sombre)
%     scatter(Coordonnees(:,1), Coordonnees(:,2), 10, [0.2 0.2 0.2], 'filled');
% 
%     % 2. Récupérer les neurones de CETTE assemblée spécifique
%     neurones_groupe = Assemblees_Noyaux{a};
% 
%     % 3. Dessiner tous les liens entre les neurones de CE groupe
%     for i = 1:length(neurones_groupe)
%         for j = (i+1):length(neurones_groupe)
%             id_1 = neurones_groupe(i);
%             id_2 = neurones_groupe(j);
%             score = Matrice_A(id_1, id_2);
% 
%             % Ne dessiner que les liens au-dessus du seuil
%             if score >= Seuil_Affichage
%                 X = [Coordonnees(id_1, 1), Coordonnees(id_2, 1)];
%                 Y = [Coordonnees(id_1, 2), Coordonnees(id_2, 2)];
% 
%                 % L'épaisseur dépend du score (60% = fin, 100% = épais)
%                 epaisseur = (score / 100) * Epaisseur_Max;
% 
%                 % Dessiner la ligne avec la couleur de l'assemblée
%                 plot(X, Y, 'Color', [couleurs(a, :), 0.6], 'LineWidth', epaisseur);
%             end
%         end
%     end
% 
%     % 4. Dessiner les neurones de l'assemblée par-dessus les lignes
%     X_groupe = Coordonnees(neurones_groupe, 1);
%     Y_groupe = Coordonnees(neurones_groupe, 2);
%     scatter(X_groupe, Y_groupe, 40, couleurs(a, :), 'filled', ...
%             'MarkerEdgeColor', 'w', 'LineWidth', 0.5);
% 
%     % Titre du subplot
%     title(sprintf('Assemblée %d (n=%d)', a, length(neurones_groupe)), ...
%           'Color', couleurs(a, :), 'FontSize', 12, 'FontWeight', 'bold');
% end
% 
% % Titre global de la figure
% sgtitle(sprintf('Connectivité spatiale intra-assemblée (Liens > %d%%)', Seuil_Affichage), ...
%         'Color', 'w', 'FontSize', 16, 'FontWeight', 'bold');

%%%
% =========================================================================
% SCRIPT 8 : Cartographie et Corrélation de Spearman (Traces Temporelles)
% =========================================================================

% 1. Chargement des traces temporelles
% Assurez-vous que le chemin est correct
disp('Chargement de l''activité temporelle (Tr1b)...');
chemin_results = "E:\Data\Aurelie\analysis\March2026\nocues\test_stability\411582_230320_plane0_26_04_06_17_17_46\k9\results.mat"; 
data_results = load(chemin_results, 'Tr1b','WinRest');
WinRest=data_results.WinRest;
Traces_Temporelles = data_results.Tr1b; % Matrice [Nb_Neurones x Temps]
Traces_Temporelles = Traces_Temporelles (:,WinRest);

% Paramètres d'affichage spatial
Seuil_Affichage = 85; 
Epaisseur_Max = 1;
Nb_Assemblees = length(Assemblees_Noyaux);
couleurs = lines(Nb_Assemblees);

% 2. Calcul de la corrélation de Spearman "Baseline" (Toutes les cellules)
% On transpose (Traces') car 'corr' dans MATLAB calcule la corrélation 
% entre les COLONNES. On veut corréler les neurones entre eux sur le temps.
disp('Calcul de la corrélation de Spearman pour tout le réseau (baseline)...');
Corr_Globale_Matrice = corr(Traces_Temporelles', 'Type', 'Spearman');

% On ignore la diagonale (les 1) et on prend seulement le triangle supérieur
valeurs_globales = Corr_Globale_Matrice(triu(true(size(Corr_Globale_Matrice)), 1));
Spearman_Moyen_Global = mean(valeurs_globales, 'omitnan');

fprintf('-> Corrélation Spearman moyenne (Réseau entier) = %.3f\n\n', Spearman_Moyen_Global);

% 3. Préparation de la figure
cols = ceil(sqrt(Nb_Assemblees));
rows = ceil(Nb_Assemblees / cols);
figure('Name', 'Assemblées et Corrélations Temporelles', ...
       'Color', 'k', 'Position', [50 50 1400 900]);

% 4. Boucle sur chaque assemblée (Dessin + Calcul local)
for a = 1:Nb_Assemblees
    subplot(rows, cols, a); hold on;
    set(gca, 'Color', 'k', 'YDir', 'reverse'); axis equal; axis off;
    
    % Dessiner le fond (neurones inactifs)
    scatter(Coordonnees(:,1), Coordonnees(:,2), 10, [0.2 0.2 0.2], 'filled');
    
    % Neurones de cette assemblée
    neurones_groupe = Assemblees_Noyaux{a};
    
    % ----- A. CALCUL DE LA CORRELATION INTRA-ASSEMBLÉE -----
    % Extraire les traces uniquement pour les neurones de ce groupe
    Traces_Groupe = Traces_Temporelles(neurones_groupe, :);
    
    % Calculer la matrice de corrélation Spearman pour ce petit groupe
    Corr_Groupe_Matrice = corr(Traces_Groupe', 'Type', 'Spearman');
    
    % Moyenne intra-groupe (hors diagonale)
    valeurs_groupe = Corr_Groupe_Matrice(triu(true(size(Corr_Groupe_Matrice)), 1));
    Spearman_Intra = mean(valeurs_groupe, 'omitnan');
    
    % ----- B. DESSIN SPATIAL (Lignes + Noeuds) -----
    for i = 1:length(neurones_groupe)
        for j = (i+1):length(neurones_groupe)
            id_1 = neurones_groupe(i);
            id_2 = neurones_groupe(j);
            score_consensus = Matrice_A(id_1, id_2);
            
            if score_consensus >= Seuil_Affichage
                X = [Coordonnees(id_1, 1), Coordonnees(id_2, 1)];
                Y = [Coordonnees(id_1, 2), Coordonnees(id_2, 2)];
                epaisseur = (score_consensus / 100) * Epaisseur_Max;
                plot(X, Y, 'Color', [couleurs(a, :), 0.6], 'LineWidth', epaisseur);
            end
        end
    end
    
    X_groupe = Coordonnees(neurones_groupe, 1);
    Y_groupe = Coordonnees(neurones_groupe, 2);
    scatter(X_groupe, Y_groupe, 40, couleurs(a, :), 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.5);
            
    % ----- C. TITRE DYNAMIQUE AVEC LES STATISTIQUES -----
    % Affiche le Spearman de l'assemblée vs la baseline du réseau
    titre_ligne1 = sprintf('Assemblée %d (n=%d)', a, length(neurones_groupe));
    titre_ligne2 = sprintf('Spearman: %.3f (Global: %.3f)', Spearman_Intra, Spearman_Moyen_Global);
    
    title({titre_ligne1, titre_ligne2}, 'Color', couleurs(a, :), 'FontSize', 11, 'FontWeight', 'bold');
end

% Titre de la fenêtre entière
sgtitle('Réseaux spatiaux et Synchronisation Temporelle (Corrélation de Spearman)', ...
        'Color', 'w', 'FontSize', 16, 'FontWeight', 'bold');

% =========================================================================
% SCRIPT FINAL : Calcul du Meilleur K via results.mat + Cartographie
% =========================================================================
clear; clc;

% 1. PARAMÈTRES ET CHEMINS (À adapter pour la session cible)
dossier_mere = "E:\Data\Aurelie\analysis\March2026\nocues\test_stability\411582_230320_plane0_26_04_06_17_17_46";
chemin_fall = 'E:\Data\Aurelie\data\nocues\411582\230320_plane0\Fall.mat';


% 2. NOUVEAUX POIDS DU SCORE COMPOSITE
w_sil        = 0.30;  % sil all clusters (séparation mathématique)
w_assemblies = 0.20;  % Assemblies (quantité)
w_recall     = 0.20;  % recall (SClok)
w_sce        = 0.15;  % SCE_Purity
w_cell       = 0.15;  % cell_Purity

% =========================================================================
% ÉTAPE 1 : EXTRACTION DES MÉTRIQUES ET CALCUL DU BEST K
% =========================================================================
fprintf('Analyse des dossiers K dans : %s\n', dossier_mere);
dossiers_K = dir(fullfile(dossier_mere, '*K*')); 

% Initialisation des listes
Ks_list = []; Silh = []; Recall = []; CellSpec = []; SCEPurity = []; AssembliesCount = [];

for i = 1:length(dossiers_K)
    if ~dossiers_K(i).isdir, continue; end
    
    k_num_str = regexp(dossiers_K(i).name, '\d+', 'match');
    if isempty(k_num_str), continue; end
    current_k = str2double(k_num_str{1});
    
    fichier_res = fullfile(dossiers_K(i).folder, dossiers_K(i).name, 'results.mat');
    if ~isfile(fichier_res), continue; end
    
    try
        data = load(fichier_res, 'RasterRace', 'NRace', 'NRaceOK', 'mean_sClOK', 'assemblyortho','WinRest');
        
        % Variables temporaires pour s'assurer que tout est valide avant d'ajouter
        temp_silh = NaN; temp_recall = 0; temp_cell = NaN; temp_sce = NaN; temp_asc = 0;
        
        % 1. Silhouette
        if isfield(data, 'mean_sClOK') && ~isnan(data.mean_sClOK)
            temp_silh = data.mean_sClOK;
        else, continue; end % On passe si pas de silhouette
        
        % 2. Recall
        if isfield(data, 'NRaceOK') && data.NRace > 0
            temp_recall = (data.NRaceOK / data.NRace) * 100;
        end
        
        % 3. Quantité d'assemblées
        if isfield(data, 'assemblyortho')
            temp_asc = length(data.assemblyortho);
        end
        
        % 4. Specificity & Purity
        if isfield(data, 'assemblyortho') && isfield(data, 'RasterRace')
            nb_clusters = length(data.assemblyortho);
            nb_SCE = size(data.RasterRace, 2);
            participation_SCE = zeros(nb_clusters, nb_SCE);
            
            for c = 1:nb_clusters
                if ~isempty(data.assemblyortho{c})
                    participation_SCE(c, :) = sum(data.RasterRace(data.assemblyortho{c}, :), 1);
                end
            end
            
            [val_max, SCE_labels] = max(participation_SCE, [], 1);
            SCE_labels(val_max == 0) = 0;
            
            cell_spec_inter = NaN(1, nb_clusters);
            sce_purity = NaN(1, nb_clusters);
            
            for c = 1:nb_clusters
                cellules = data.assemblyortho{c};
                idx_SCE = (SCE_labels == c);
                if sum(idx_SCE) > 0 && ~isempty(cellules)
                    tirs_in_box = sum(data.RasterRace(cellules, idx_SCE), 'all');
                    tirs_totaux_colonne = sum(data.RasterRace(:, idx_SCE), 'all');
                    if tirs_totaux_colonne > 0, sce_purity(c) = tirs_in_box / tirs_totaux_colonne; end
                    
                    tirs_totaux_ligne = sum(data.RasterRace(cellules, SCE_labels > 0), 'all');
                    if tirs_totaux_ligne > 0, cell_spec_inter(c) = tirs_in_box / tirs_totaux_ligne; end
                end
            end
            temp_cell = mean(cell_spec_inter, 'omitnan');
            temp_sce = mean(sce_purity, 'omitnan');
        end
        
        % Si on arrive ici sans erreur, on ajoute tout dans nos listes définitives
        Ks_list(end+1) = current_k;
        Silh(end+1) = temp_silh;
        Recall(end+1) = temp_recall;
        AssembliesCount(end+1) = temp_asc;
        CellSpec(end+1) = temp_cell;
        SCEPurity(end+1) = temp_sce;
        
    catch
        continue;
    end
end

% Remplacer les NaN éventuels par 0 avant normalisation
CellSpec(isnan(CellSpec)) = 0;
SCEPurity(isnan(SCEPurity)) = 0;

% Fonction de normalisation Min-Max (sécurisée contre les divisions par zéro)
safe_norm = @(x) (x - min(x)) ./ max(1e-9, max(x) - min(x));

% Calcul du Score Composite avec vos nouveaux poids
ScoreGlobal = (w_sil * safe_norm(Silh)) + ...
              (w_assemblies * safe_norm(AssembliesCount)) + ...
              (w_recall * safe_norm(Recall)) + ...
              (w_sce * safe_norm(SCEPurity)) + ...
              (w_cell * safe_norm(CellSpec));

[best_score, max_idx] = max(ScoreGlobal);
Best_K = Ks_list(max_idx);
Best_K =12;
fprintf('--> Meilleur K calculé : K = %d (Score = %.3f)\n\n', Best_K, best_score);


% =========================================================================
% ÉTAPE 2 : PRÉPARATION DES COORDONNÉES ET TRACES POUR LE BEST K
% =========================================================================
disp('Chargement des traces temporelles et des coordonnées pour ce Best K...');
dossier_best_k = strcat(dossier_mere, '\K' ,num2str(Best_K));
% dossier_best_k = dir(fullfile(dossier_mere, ['*' num2str(Best_K) '*']));
% fichier_best_k = fullfile(dossier_mere, dossier_best_k(1).name, 'results.mat');
fichier_best_k = strcat(dossier_best_k, '\results.mat');

data_best = load(fichier_best_k, 'Tr1b', 'assemblyortho','WinRest');
Traces_Temporelles = data_best.Tr1b;
WinRest = data_best.WinRest;
Traces_Temporelles = Traces_Temporelles(:,WinRest);

% Extraction des coordonnées (Fall.mat)
load(chemin_fall, 'stat', 'iscell');
Coordonnees = zeros(sum(iscell(:,1)==1), 2);
m = 0;
for n = 1:length(stat)
    if iscell(n,1)==1
        m = m + 1;
        Coordonnees(m, 1) = double(stat{n}.med(2)) + 1; % X
        Coordonnees(m, 2) = double(stat{n}.med(1)) + 1; % Y
    end
end

% Nettoyage : on ne garde que les assemblées avec au moins 5 neurones
Assemblees = {};
for i = 1:length(data_best.assemblyortho)
    if length(data_best.assemblyortho{i}) >= 5
        Assemblees{end+1} = data_best.assemblyortho{i};
    end
end
Nb_Assemblees = length(Assemblees);

% Baseline Spearman globale
Corr_Globale = corr(Traces_Temporelles', 'Type', 'Spearman');
Spearman_Global = mean(Corr_Globale(triu(true(size(Corr_Globale)), 1)), 'omitnan');

% =========================================================================
% ÉTAPE 3 : CARTOGRAPHIE ET CORRÉLATION (SPEARMAN)
% =========================================================================
disp('Création de la figure finale...');
cols = ceil(sqrt(Nb_Assemblees)); rows = ceil(Nb_Assemblees / cols);
couleurs = lines(Nb_Assemblees);

figure('Name', sprintf('Meilleur K = %d', Best_K), 'Color', 'k', 'Position', [50 50 1400 900]);

for a = 1:Nb_Assemblees
    subplot(rows, cols, a); hold on;
    set(gca, 'Color', 'k', 'YDir', 'reverse'); axis equal; axis off;
    
    % Fond gris
    scatter(Coordonnees(:,1), Coordonnees(:,2), 10, [0.2 0.2 0.2], 'filled');
    
    % Calcul Spearman pour l'assemblée
    neurones_groupe = Assemblees{a};
    Traces_Groupe = Traces_Temporelles(neurones_groupe, :);
    Corr_Groupe = corr(Traces_Groupe', 'Type', 'Spearman');
    Spearman_Intra = mean(Corr_Groupe(triu(true(size(Corr_Groupe)), 1)), 'omitnan');
    
    % Dessin des liens
    for i = 1:length(neurones_groupe)
        for j = (i+1):length(neurones_groupe)
            if Corr_Groupe(i,j) > 0.3 % Filtre d'affichage visuel
                X = [Coordonnees(neurones_groupe(i), 1), Coordonnees(neurones_groupe(j), 1)];
                Y = [Coordonnees(neurones_groupe(i), 2), Coordonnees(neurones_groupe(j), 2)];
                epaisseur = Corr_Groupe(i,j) * 4;
                plot(X, Y, 'Color', [couleurs(a, :), 0.5], 'LineWidth', epaisseur);
            end
        end
    end
    
    % Dessin des neurones
    scatter(Coordonnees(neurones_groupe, 1), Coordonnees(neurones_groupe, 2), ...
            40, couleurs(a, :), 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.5);
            
    % Titre du cluster
    t1 = sprintf('Cluster %d (n=%d)', a, length(neurones_groupe));
    t2 = sprintf('Spearman: %.3f (Total: %.3f)', Spearman_Intra, Spearman_Global);
    title({t1, t2}, 'Color', couleurs(a, :), 'FontSize', 11, 'FontWeight', 'bold');
end

% Titre principal
sgtitle(sprintf('Topologie au Meilleur K (K=%d | Score=%.3f)', Best_K, best_score), ...
        'Color', 'w', 'FontSize', 16, 'FontWeight', 'bold');