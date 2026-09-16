% =========================================================================
% SCRIPT analyse_assembly_results5parfor.mlx:
%  EXTRACT_BEST_K_FROM_RESULTS + GRAPHIQUES + SUCCESSRATE
% Objectif : Trouver les meilleurs K, extraire le Success Rate, et faire les graphes
% =========================================================================

clear; clc; close all

%% 1. PARAMÈTRES ET CHEMINS
dossier_mere = "E:\Data\Aurelie\analysis\Final\assembly\nocues\";

nbtop = 3; %nb K à sortir
k_start = 4;

%% 2. POIDS DES CRITÈRES POUR LE SCORE GLOBAL
weight_silhouette = 0;
weight_recall = 0.5; %proportion de SCE clusterisés
weight_cell_specificity = 0.3;  %proportion des activations d'une cellule qui tombent dans les SCE de son propre cluster, parmi tous les SCE valides 
weight_sce_purity = 0.20; %proportion d'un SCE donné qui est composée des cellules du cluster censé l'expliquer, parmi toutes les cellules du raster 



%% 3. RECHERCHE DE TOUS LES FICHIERS RESULTS.MAT
% fprintf('Recherche des fichiers results.mat...\n');
tous_les_fichiers = dir(fullfile(dossier_mere, '**', 'results.mat'));

% if isempty(tous_les_fichiers)
%     error('Aucun fichier results.mat trouvé.');
% end
% sessions_autorisees = {'411582_230320_plane0', '411582_230320_plane1', '443669_220614_plane0', '443669_220614_plane1', ...
%     '443695_220728_plane0', '443695_220809_plane1','444113_220927_plane1', '444113_220928_plane1', ...
%     '444119_220919_plane0', '444119_220919_plane1','444152_221129_plane0', '444152_221129_plane1', ...
%     '444175_221121_plane0', '444175_221121_plane1', '444175_221125_plane0', '444175_221125_plane1', ...
%     '446178_11062025_plane0', '446178_11062025_plane1', '446198_12062025_plane0','446178_14062025_plane0',...
% '446178_14062025_plane1'...
% };
% 
% % Filtrage vectoriel direct et réécriture de la variable
% tous_les_fichiers = tous_les_fichiers(contains({tous_les_fichiers.folder}, sessions_autorisees));


% fprintf('%d fichiers trouvés. Analyse en cours...\n\n', length(tous_les_fichiers));

%% 4. EXTRACTION ET CALCULS DES MÉTRIQUES (PARALLÉLISÉ)
% fprintf('Démarrage du pool parallèle...\n');
% Vérifie si un pool est ouvert, sinon en ouvre un
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool;
end

N_fichiers = length(tous_les_fichiers);

% Pré-allocation des variables pour le parfor
temp_sessionName     = cell(N_fichiers, 1);
temp_safeSessionName = cell(N_fichiers, 1);
temp_k               = NaN(N_fichiers, 1);
temp_silhouette      = NaN(N_fichiers, 1);
temp_sce_recall      = NaN(N_fichiers, 1);
temp_cell_spec_inter = NaN(N_fichiers, 1);
temp_cell_spec_global= NaN(N_fichiers, 1);
temp_sce_purity      = NaN(N_fichiers, 1);
temp_success_rate    = NaN(N_fichiers, 1);
temp_paths           = cell(N_fichiers, 1);

% Configuration de la barre de progression pour parfor
h_wait = waitbar(0, 'Extraction des données (Parallel)...', 'Name', 'Progression');
D = parallel.pool.DataQueue;
afterEach(D, @(~) updateWaitbar(h_wait, N_fichiers));

fprintf('Début de l''extraction en parallèle...\n');

parfor i = 1:N_fichiers
    % --- Initialisation locale ---
    matFilePath = fullfile(tous_les_fichiers(i).folder, tous_les_fichiers(i).name);
    parts = strsplit(tous_les_fichiers(i).folder, filesep);
    
    if length(parts) < 2
        send(D, i); continue;
    end
    
    kFolderName = parts{end};
    sessName = parts{end-1};
    
    k_num_str = regexp(kFolderName, '\d+', 'match');
    if isempty(k_num_str)
        send(D, i); continue; 
    end
    
    c_k = str2double(k_num_str{1});
    if c_k < k_start
        send(D, i); continue;
    end

    % --- Chargement des données (.mat) ---
    try
        data = load(matFilePath, 'RasterRace', 'NRace', 'NRaceOK', 'mean_sClOK', 'assemblyortho');
    catch
        send(D, i); continue;
    end
    
    % 1. Silhouette
    if isfield(data, 'mean_sClOK') && isnumeric(data.mean_sClOK) && ~isnan(data.mean_sClOK)
        c_silhouette = data.mean_sClOK;
    else
        send(D, i); continue;
    end
    
    % 2. Recall
    if isfield(data, 'NRaceOK') && isfield(data, 'NRace') && data.NRace > 0
        c_sce_recall = (data.NRaceOK / data.NRace) * 100;
    else
        c_sce_recall = 0;
    end
    
    % 3. Specificity & Purity
    c_cell_spec_inter = NaN;
    c_cell_spec_global = NaN;
    c_sce_purity = NaN;
    
    if isfield(data, 'assemblyortho') && ~isempty(data.assemblyortho) && isfield(data, 'RasterRace')
        nb_clusters = length(data.assemblyortho);
        nb_SCE = size(data.RasterRace, 2);
        
        if nb_clusters > 0 && nb_SCE > 0
            participation_SCE = zeros(nb_clusters, nb_SCE);
            for c = 1:nb_clusters
                cellules_c = data.assemblyortho{c};
                if ~isempty(cellules_c)
                    participation_SCE(c, :) = sum(data.RasterRace(cellules_c, :), 1);
                end
            end
            
            [val_max, SCE_labels] = max(participation_SCE, [], 1);
            SCE_labels(val_max == 0) = 0;
            
            cell_spec_inter_array = NaN(1, nb_clusters);
            cell_spec_global_array = NaN(1, nb_clusters);
            sce_purity_array = NaN(1, nb_clusters);
            
            for c = 1:nb_clusters
                cellules_du_cluster = data.assemblyortho{c};
                idx_SCE_du_cluster = (SCE_labels == c);
                
                if sum(idx_SCE_du_cluster) > 0 && ~isempty(cellules_du_cluster)
                    tirs_in_box = sum(data.RasterRace(cellules_du_cluster, idx_SCE_du_cluster), 'all');
                    
                    tirs_totaux_colonne = sum(data.RasterRace(:, idx_SCE_du_cluster), 'all');
                    if tirs_totaux_colonne > 0
                        sce_purity_array(c) = tirs_in_box / tirs_totaux_colonne;
                    end
                    
                    idx_SCE_valides = (SCE_labels > 0);
                    tirs_totaux_ligne_inter = sum(data.RasterRace(cellules_du_cluster, idx_SCE_valides), 'all');
                    if tirs_totaux_ligne_inter > 0
                        cell_spec_inter_array(c) = tirs_in_box / tirs_totaux_ligne_inter;
                    end
                    
                    tirs_totaux_ligne_global = sum(data.RasterRace(cellules_du_cluster, :), 'all');
                    if tirs_totaux_ligne_global > 0
                        cell_spec_global_array(c) = tirs_in_box / tirs_totaux_ligne_global;
                    end
                end
            end
            c_cell_spec_inter = mean(cell_spec_inter_array, 'omitnan');
            c_cell_spec_global = mean(cell_spec_global_array, 'omitnan');
            c_sce_purity = mean(sce_purity_array, 'omitnan');
        end
    end
    
    % 4. Lecture du Brainbow (Success Rate)
    chemin_brainbow = fullfile(tous_les_fichiers(i).folder, 'brainbow.mat');
    c_success_rate = NaN;
    if isfile(chemin_brainbow)
        try
            donnees_brainbow = load(chemin_brainbow, 'successrate');
            if isfield(donnees_brainbow, 'successrate')
                c_success_rate = donnees_brainbow.successrate;
            end
        catch
        end
    end
    
    % --- Enregistrement dans les tableaux temporaires ---
    if ~isnan(c_cell_spec_inter)
        temp_sessionName{i} = sessName;
        temp_safeSessionName{i} = matlab.lang.makeValidName(sessName);
        temp_k(i) = c_k;
        temp_silhouette(i) = c_silhouette;
        temp_sce_recall(i) = c_sce_recall;
        temp_cell_spec_inter(i) = c_cell_spec_inter;
        temp_cell_spec_global(i) = c_cell_spec_global;
        temp_sce_purity(i) = c_sce_purity;
        temp_success_rate(i) = c_success_rate;
        temp_paths{i} = tous_les_fichiers(i).folder;
    end
    
    send(D, i); % Met à jour la barre de progression
end
close(h_wait);

% --- Regroupement rapide des résultats dans la structure SessionsData ---
SessionsData = struct();
for i = 1:N_fichiers
    safeName = temp_safeSessionName{i};
    if ~isempty(safeName)
        if ~isfield(SessionsData, safeName)
            SessionsData.(safeName).OriginalName = temp_sessionName{i};
            SessionsData.(safeName).K = [];
            SessionsData.(safeName).Silh = [];
            SessionsData.(safeName).Recall = [];
            SessionsData.(safeName).CellSpecInter = [];
            SessionsData.(safeName).CellSpecGlobal = [];
            SessionsData.(safeName).SCEPurity = [];
            SessionsData.(safeName).SuccessRate = [];
            SessionsData.(safeName).Paths = {};
            SessionsData.(safeName).Scores = [];
        end
        
        SessionsData.(safeName).K(end+1) = temp_k(i);
        SessionsData.(safeName).Silh(end+1) = temp_silhouette(i);
        SessionsData.(safeName).Recall(end+1) = temp_sce_recall(i);
        SessionsData.(safeName).CellSpecInter(end+1) = temp_cell_spec_inter(i);
        SessionsData.(safeName).CellSpecGlobal(end+1) = temp_cell_spec_global(i);
        SessionsData.(safeName).SCEPurity(end+1) = temp_sce_purity(i);
        SessionsData.(safeName).SuccessRate(end+1) = temp_success_rate(i);
        SessionsData.(safeName).Paths{end+1} = temp_paths{i};
    end
end

%% 5. CALCUL DES SCORES ET PRÉPARATION DE LA TABLE EXCEL
sessionFields = fieldnames(SessionsData);
fprintf('--------------------------------------------------\n');

Excel_Session = {};
Excel_Rang = [];
Excel_K = [];
Excel_Score = [];
Excel_Silh = [];
Excel_Recall = [];
Excel_CellSpecInter = [];
Excel_CellSpecGlobal = [];
Excel_SCEPurity = [];
Excel_SuccessRate = []; 
Excel_Path = {};

for s = 1:length(sessionFields)
    sess = SessionsData.(sessionFields{s});
    if length(sess.K) < 2, continue; end
    
    % 1. Silhouette : De [-1, 1] vers [0, 1]
    norm_sil = (sess.Silh + 1) / 2;
    % 2. Recall : Si sess.Recall est en % (0-100), on le divise par 100 pour l'avoir sur [0, 1]
    norm_rec = sess.Recall / 100;
    % 3. Specificity et Purity : Sont déjà des ratios entre 0 et 1
    norm_cell_spec = sess.CellSpecInter;
    norm_purity = sess.SCEPurity;
    
    % SECURITE : On s'assure qu'aucun score ne dépasse 1 ou ne descend sous 0
    norm_sil = max(0, min(1, norm_sil));
    norm_rec = max(0, min(1, norm_rec));
    norm_cell_spec = max(0, min(1, norm_cell_spec));
    norm_purity = max(0, min(1, norm_purity));

    % Calcul du Score Global
    Scores = (weight_silhouette * norm_sil) + ...
             (weight_recall * norm_rec) + ...
             (weight_cell_specificity * norm_cell_spec) + ...
             (weight_sce_purity * norm_purity);
             
    % SAUVEGARDE DES SCORES POUR LES GRAPHIQUES
    SessionsData.(sessionFields{s}).Scores = Scores;
    [sortedScores, sortIdx] = sort(Scores, 'descend');
    top_k = min(nbtop, length(sortIdx)); % Affichage/Export des N meilleurs
    
    fprintf('SESSION : %s\n', sess.OriginalName);
    for i = 1:top_k
        orig_idx = sortIdx(i);
        fprintf(' ⭐ #%d -> K = %d\n', i, sess.K(orig_idx));
        fprintf(' Chemin : %s\n', sess.Paths{orig_idx});
        fprintf(' (Score: %.2f | Silh: %.3f | Recall: %.1f%% | Cell Spec: %.3f | Succes: %.2f)\n', ...
            sortedScores(i), sess.Silh(orig_idx), sess.Recall(orig_idx), sess.CellSpecInter(orig_idx), sess.SuccessRate(orig_idx));
            
        Excel_Session{end+1, 1} = sess.OriginalName;
        Excel_Rang(end+1, 1) = i;
        Excel_K(end+1, 1) = sess.K(orig_idx);
        Excel_Score(end+1, 1) = sortedScores(i);
        Excel_Silh(end+1, 1) = sess.Silh(orig_idx);
        Excel_Recall(end+1, 1) = sess.Recall(orig_idx);
        Excel_CellSpecInter(end+1, 1) = sess.CellSpecInter(orig_idx);
        Excel_CellSpecGlobal(end+1, 1) = sess.CellSpecGlobal(orig_idx);
        Excel_SCEPurity(end+1, 1) = sess.SCEPurity(orig_idx);
        Excel_SuccessRate(end+1, 1) = sess.SuccessRate(orig_idx); 
        Excel_Path{end+1, 1} = sess.Paths{orig_idx};
    end
    fprintf('--------------------------------------------------\n');
end
%%

%% 6. SAUVEGARDE DANS LE FICHIER EXCEL (analysis_brainbow_best.xlsx)
fprintf('\n--------------------------------------------------\n');
fprintf('Génération du fichier analysis_brainbow_best.xlsx...\n');

path_global = fullfile(dossier_mere, 'analysis_brainbow_best.xlsx');

% --- En-tête mis à jour (ajout de Score_Global et Chemin_Dossier) ---
header_global = {'date', 'identifier', 'Score_Global','SCE', 'sce freq', 'SCEok', 'SCEok freq', ...
'Assemblies', ' % colored', 'sil all clusters','SClok','SClokfinal', 'cell/assembly', ...
'median cell/assembly', 'Ncell', 'frames','sampling rate', 'nb frame not running', ...
'% not running', 'average speed', 'max speed', 'nblap', 'sce_n_cells_threshold', ...
'mean distance assemblies', 'mean distance all cells', 'n_kmean_surrogate', 'file_num','NCl_ini',  'Chemin_Dossier'};

% Si le fichier n'existe pas, on l'initialise avec l'en-tête
if ~isfile(path_global)
    writecell(header_global, path_global, 'Range', 'A1');
    numline_global = 2;
else
    data_excel = readcell(path_global);
    numline_global = size(data_excel, 1) + 1;
end

compteur_excel = 0;

for s = 1:length(sessionFields)
    sess = SessionsData.(sessionFields{s});
    if length(sess.K) < 2 || ~isfield(sess, 'Scores'), continue; end
    
    % On trie les K par score décroissant
    [~, sortIdx] = sort(sess.Scores, 'descend');
    top_k = min(nbtop, length(sortIdx)); % On prend les nbtop (ex: 3) meilleurs
    
    % --- BOUCLE SUR LES N MEILLEURS K ---
    for i = 1:top_k
        orig_idx = sortIdx(i); 
        
        best_k_path = sess.Paths{orig_idx};
        matFilePath = fullfile(best_k_path, 'results.mat');
        brainbowmatFilePath = fullfile(best_k_path, 'brainbow.mat');
        
        % On charge les deux fichiers .mat correspondants pour ce K
        if isfile(matFilePath) && isfile(brainbowmatFilePath)
            try
                % Chargement des variables nécessaires depuis results.mat
                load(matFilePath, 'assemblyortho');
                
                % Calcul de writeassemblyortho et idortho
                if isempty(assemblyortho)
                    writeassemblyortho = 0;
                    idortho = [];
                else
                    writeassemblyortho = length(assemblyortho);
                    idortho = cellfun(@length, assemblyortho);
                end
                
                % Chargement des variables biologiques depuis brainbow.mat
                load(brainbowmatFilePath, ...
                    'daytime', 'identifier', 'NRace', 'freqsce', 'NRaceOK', 'freqsceok', ...
                    'successrate', 'sCl', 'NClOK', 'mean_sClOK', 'mean_cell_assembly', ...
                    'NCell', 'Nz', 'sampling_rate', 'WinRest', 'average_speed', 'max_speed', ...
                    'nb_lap', 'sce_n_cells_threshold', 'mean_distance_assembly', ...
                    'meandistance_gcamp', 'kmeans_surrogate', 'file_num', 'NClini');
                
                % Ajout du "K" au nom de l'identifier pour le différencier dans Excel
                identifier_avec_k = sprintf('%s_K%d', identifier, sess.K(orig_idx));
                
                % --- Construction de la ligne de données ---
                data_global = {daytime, identifier_avec_k, sess.Scores(orig_idx), NRace, freqsce, NRaceOK, freqsceok, ...
                    writeassemblyortho, successrate, mean(sCl), mean(sCl(1:NClOK)), mean_sClOK, num2str(idortho), ...
                    median(idortho), NCell, Nz, sampling_rate, length(WinRest), length(WinRest)/Nz*100, ...
                    average_speed, max_speed, nb_lap, sce_n_cells_threshold, mean_distance_assembly, ...
                    meandistance_gcamp, kmeans_surrogate, file_num, NClini,  best_k_path};
                    
                % Écriture dans le fichier Excel
                writecell(data_global, path_global, 'Range', ['A' num2str(numline_global)]);
                
                numline_global = numline_global + 1;
                compteur_excel = compteur_excel + 1;
                
            catch exception
                fprintf('Erreur de chargement sur le K=%d de la session %s\n', sess.K(orig_idx), sess.OriginalName);
            end
        end
    end
end

fprintf('✅ %d lignes (top %d par session) ont été ajoutées dans le fichier Excel :\n-> %s\n', compteur_excel, nbtop, path_global);
%%

%% 7. GRAPHIQUES : Silhouette et Score Global vs K (TRIÉS PAR K CROISSANT)
figure('Name', 'Analyse des métriques par K', 'Color', 'w', 'Position', [100 100 1200 600]);

sessionFields = fieldnames(SessionsData);
colors = lines(length(sessionFields));

% --- GRAPHIQUE 1 : Silhouette moyenne en fonction de K ---
figure
hold on;
for s = 1:length(sessionFields)
    sess = SessionsData.(sessionFields{s});
    if length(sess.K) >= 2
        [K_sorted, idx_sort] = sort(sess.K);
        Silh_sorted = sess.Silh(idx_sort);
        plot(K_sorted, Silh_sorted, '-o', 'Color', colors(s,:), 'LineWidth', 2.5, ...
            'MarkerSize', 7, 'MarkerFaceColor', colors(s,:)*0.8, ...
            'DisplayName', sess.OriginalName);
    end
end
xlabel('Nombre de clusters (K)', 'FontWeight', 'bold');
ylabel('Score Silhouette', 'FontWeight', 'bold');
title('Silhouette moyenne selon K');
grid on; grid minor;
%legend('Location', 'best', 'Interpreter', 'none', 'FontSize', 8);
hold off;
nom_graph = fullfile(dossier_mere, 'silhouette.png');
saveas(gcf, nom_graph);
%%
% --- GRAPHIQUE 2 : Score Global en fonction de K ---
figure
hold on;

% 1. On trace d'abord toutes les sessions individuelles
all_k_for_mean = [];
all_scores_for_mean = [];

for s = 1:length(sessionFields)
    sess = SessionsData.(sessionFields{s});
    if isfield(sess, 'Scores') && length(sess.K) >= 2
        [K_sorted, idx_sort] = sort(sess.K);
        Scores_sorted = sess.Scores(idx_sort);
        
        % On accumule les données pour le calcul de la moyenne plus tard
        all_k_for_mean = [all_k_for_mean, K_sorted];
        all_scores_for_mean = [all_scores_for_mean, Scores_sorted];
        
        % Tracé individuel
        plot(K_sorted, Scores_sorted, '-s', 'Color', [colors(s,:) 0.4], 'LineWidth', 1.5, ...
            'MarkerSize', 6, 'MarkerFaceColor', colors(s,:), ...
            'DisplayName', sess.OriginalName);
    end
end

% 2. Calcul et tracé de la courbe moyenne
if ~isempty(all_k_for_mean)
    unique_k = unique(all_k_for_mean);
    mean_scores = zeros(1, length(unique_k));
    
    for idx_k = 1:length(unique_k)
        current_k_val = unique_k(idx_k);
        % Moyenne de tous les scores obtenus pour ce K précis
        mean_scores(idx_k) = mean(all_scores_for_mean(all_k_for_mean == current_k_val), 'omitnan');
    end
    
    % Tracé de la moyenne en noir et épais pour qu'elle ressorte
    plot(unique_k, mean_scores, '-k*', 'LineWidth', 4, 'MarkerSize', 10, ...
        'DisplayName', 'Moyenne Globale');
end

xlabel('Nombre de clusters (K)', 'FontWeight', 'bold');
ylabel('Score Global', 'FontWeight', 'bold');
title('Score Global selon K');
grid on; grid minor;
% legend('Location', 'best', 'Interpreter', 'none', 'FontSize', 8);
hold off;

nom_graph = fullfile(dossier_mere, 'Score_global_K.png');
saveas(gcf, nom_graph);
%%


% --- GRAPHIQUE 3 : Silhouette vs Score Global (scatter) ---
figure;
subplot(1, 2, 1);
hold on;
all_silh = [];
all_scores = [];
for s = 1:length(sessionFields)
    sess = SessionsData.(sessionFields{s});
    if length(sess.K) >= 2
        all_silh = [all_silh, sess.Silh];
        all_scores = [all_scores, sess.Scores];
    end
end
scatter(all_silh, all_scores, 80, 'filled', 'MarkerEdgeColor', 'k');
xlabel('Score Silhouette', 'FontWeight', 'bold');
ylabel('Score Global', 'FontWeight', 'bold');
title('Silhouette vs Score Global');
grid on;
hold off;

% --- GRAPHIQUE 4 : Top 3 K recommandés ---
subplot(1, 2, 2);
top3_k = [];
top3_scores = [];
for s = 1:length(sessionFields)
    sess = SessionsData.(sessionFields{s});
    if length(sess.K) >= 2 && isfield(sess, 'Scores')
        [~, sortIdx] = sort(sess.Scores, 'descend');
        top3_k = [top3_k, sess.K(sortIdx(1:min(3, end)))];
        top3_scores = [top3_scores, sess.Scores(sortIdx(1:min(3, end)))];
    end
end
if ~isempty(top3_k)
    bar(1:length(top3_k), top3_scores, 'FaceColor', [0.2 0.6 0.8]);
    xticks(1:length(top3_k));
    xticklabels(arrayfun(@num2str, top3_k, 'UniformOutput', false));
    xlabel('Top 3 K recommandés', 'FontWeight', 'bold');
    ylabel('Score Global', 'FontWeight', 'bold');
    title('Meilleurs K (toutes sessions)');
    grid on;
end

% Sauvegarde globale
nom_graph = fullfile(dossier_mere, 'Analyse_Performance_K.png');
saveas(gcf, nom_graph);
fprintf('Graphique sauvegardé : %s\n', nom_graph);


%% 8. FONCTION LOCALE (POUR LA WAITBAR DU PARFOR)
% Doit obligatoirement être tout à la fin du fichier
function updateWaitbar(h, total)
    persistent count;
    if isempty(count), count = 0; end
    count = count + 1;
    if isgraphics(h)
        waitbar(count / total, h, sprintf('Analyse Parallèle : %d/%d', count, total));
    end
    if count == total, count = []; end % Reset à la fin
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
