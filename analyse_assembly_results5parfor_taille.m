%% ========================================================================
% SCRIPT : analyse_assembly_results5parfor_taille.m
% EXTRACT_BEST_K_FROM_RESULTS + GRAPHIQUES + SUCCESS RATE + TAILLE D'ASSEMBLEE
%
% Objectif : trouver le meilleur K par session, en incluant désormais un
% critère de taille d'assemblée (exprimé en % du réseau imagé, plus robuste
% que le nombre absolu de cellules car il s'adapte à la taille du champ
% imagé, qui varie fortement d'une session à l'autre - cf discussion sur
% Malvache et al. : ~20-30 cellules sur ~200 cellules totales, soit ~10-15%
% du réseau).
% =========================================================================

clear; clc; close all;

%% 1. PARAMETRES ET CHEMINS
dossier_mere = 'E:\';
nb_top       = 5;      % nb de meilleurs K à sortir par session
k_start      = 4;

%% 2. POIDS DES CRITERES POUR LE SCORE GLOBAL
weight_silhouette       = 0;     % volontairement à 0 : la silhouette diminue quasi
                                  % monotonement avec K sur ce dataset, elle biaiserait
                                  % systématiquement vers les K faibles
weight_recall           = 0.30;  % proportion de SCE clusterisés (SCEok / SCE)
weight_cell_specificity = 0.20;  % proportion des activations d'une cellule qui tombent
                                  % dans les SCE de son propre cluster, parmi tous les
                                  % SCE valides
weight_sce_purity       = 0.15;  % proportion d'un SCE donné composée des cellules du
                                  % cluster censé l'expliquer, parmi toutes les cellules
                                  % du raster
weight_size             = 0.35;  % NOUVEAU : score de taille d'assemblée (voir sizeScore)

% --- Cible de taille d'assemblée, en % du réseau total (Ncell) ---
% Fenêtre "idéale" 10-15% (cohérent avec Malvache et al.), avec décroissance
% linéaire jusqu'à 0 en dehors de la fenêtre élargie [size_pct_floor ; size_pct_ceiling]
size_pct_target_min = 10;   % %
size_pct_target_max = 15;   % %
size_pct_floor      = 4;    % en dessous, score de taille = 0
size_pct_ceiling    = 30;   % au dessus, score de taille = 0

%% 3. RECHERCHE DE TOUS LES FICHIERS RESULTS.MAT
fprintf('Recherche des fichiers results.mat...\n');
tous_les_fichiers = dir(fullfile(dossier_mere, '**', 'results.mat'));
if isempty(tous_les_fichiers)
    error('Aucun fichier results.mat trouvé.');
end

sessions_autorisees = {'411582_230320_plane0','411582_230320_plane1','443669_220614_plane0','443669_220614_plane1', ...
    '443695_220728_plane0','443695_220809_plane1','444113_220927_plane1','444113_220928_plane1', ...
    '444119_220919_plane0','444119_220919_plane1','444152_221129_plane0','444152_221129_plane1', ...
    '444175_221121_plane0','444175_221121_plane1','444175_221125_plane0','444175_221125_plane1', ...
    '4461781_1062025_plane0','4461781_1062025_plane1','4461981_2062025_plane0','4461781_4062025_plane0', ...
    '4461781_4062025_plane1'};

% Filtrage vectoriel direct et réécriture de la variable tous_les_fichiers
tous_les_fichiers = tous_les_fichiers(contains({tous_les_fichiers.folder}, sessions_autorisees));
fprintf('%d fichiers trouvés. Analyse en cours...\n', length(tous_les_fichiers));

%% 4. EXTRACTION ET CALCULS DES METRIQUES (PARALLELISE)
fprintf('Démarrage du pool parallèle...\n');
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool;
end
Nfichiers = length(tous_les_fichiers);

% --- Pré-allocation des variables pour le parfor ---
temp_sessionName     = cell(Nfichiers, 1);
temp_safeSessionName = cell(Nfichiers, 1);
temp_k               = nan(Nfichiers, 1);
temp_silhouette      = nan(Nfichiers, 1);
temp_sce_recall      = nan(Nfichiers, 1);
temp_cellspec_inter  = nan(Nfichiers, 1);
temp_cellspec_global = nan(Nfichiers, 1);
temp_sce_purity      = nan(Nfichiers, 1);
temp_success_rate    = nan(Nfichiers, 1);
temp_mean_cell_assembly = nan(Nfichiers, 1);   % NOUVEAU
temp_ncell           = nan(Nfichiers, 1);      % NOUVEAU
temp_pct_network     = nan(Nfichiers, 1);      % NOUVEAU
temp_paths           = cell(Nfichiers, 1);

hwait = waitbar(0, 'Extraction des données (Parallel)...', 'Name', 'Progression');
D = parallel.pool.DataQueue;
afterEach(D, @(~) updateWaitbar(hwait, Nfichiers));

fprintf('Début de l''extraction en parallèle...\n');
parfor i = 1:Nfichiers
    % --- Initialisation locale ---
    matFilePath = fullfile(tous_les_fichiers(i).folder, tous_les_fichiers(i).name);
    parts = strsplit(tous_les_fichiers(i).folder, filesep);
    if length(parts) < 2, send(D,i); continue; end
    kFolderName = parts{end};
    sessName    = parts{end-1};
    knumstr = regexp(kFolderName, '\d+', 'match');
    if isempty(knumstr), send(D,i); continue; end
    ck = str2double(knumstr{1});
    if ck < k_start, send(D,i); continue; end

    % --- Chargement des données results.mat ---
    try
        data = load(matFilePath, 'RasterRace', 'NRace', 'NRaceOK', 'meansClOK', 'assemblyortho');
    catch
        send(D,i); continue;
    end

    % 1. Silhouette
    if isfield(data,'meansClOK') && isnumeric(data.meansClOK) && ~isnan(data.meansClOK)
        c_silhouette = data.meansClOK;
    else
        send(D,i); continue;
    end

    % 2. Recall
    if isfield(data,'NRaceOK') && isfield(data,'NRace') && data.NRace > 0
        c_sce_recall = data.NRaceOK / data.NRace * 100;
    else
        c_sce_recall = 0;
    end

    % 3. Specificité / Pureté + Taille moyenne d'assemblée
    c_cellspec_inter  = NaN;
    c_cellspec_global = NaN;
    c_sce_purity      = NaN;
    c_mean_cell_assembly = NaN;

    if isfield(data,'assemblyortho') && ~isempty(data.assemblyortho) && isfield(data,'RasterRace')
        nb_clusters = length(data.assemblyortho);
        nb_SCE = size(data.RasterRace, 2);

        % --- taille moyenne d'assemblée ---
        idortho = cellfun(@length, data.assemblyortho);
        c_mean_cell_assembly = mean(idortho, 'omitnan');

        if nb_clusters > 0 && nb_SCE > 0
            participationSCE = zeros(nb_clusters, nb_SCE);
            for c = 1:nb_clusters
                cellules_c = data.assemblyortho{c};
                if ~isempty(cellules_c)
                    participationSCE(c,:) = sum(data.RasterRace(cellules_c,:,1), 1);
                end
            end
            [valmax, SCElabels] = max(participationSCE, [], 1);
            SCElabels(valmax==0) = 0;

            cellspecinter_array  = nan(1, nb_clusters);
            cellspecglobal_array = nan(1, nb_clusters);
            scepurity_array      = nan(1, nb_clusters);

            for c = 1:nb_clusters
                cellules_du_cluster = data.assemblyortho{c};
                idxSCE_du_cluster   = (SCElabels == c);
                if sum(idxSCE_du_cluster) > 0 && ~isempty(cellules_du_cluster)
                    tirs_in_box = sum(data.RasterRace(cellules_du_cluster, idxSCE_du_cluster, :), 'all');

                    tirs_totaux_colonne = sum(data.RasterRace(:, idxSCE_du_cluster, :), 'all');
                    if tirs_totaux_colonne > 0
                        scepurity_array(c) = tirs_in_box / tirs_totaux_colonne;
                    end

                    idxSCE_valides = (SCElabels ~= 0);
                    tirs_totaux_ligne_inter = sum(data.RasterRace(cellules_du_cluster, idxSCE_valides, :), 'all');
                    if tirs_totaux_ligne_inter > 0
                        cellspecinter_array(c) = tirs_in_box / tirs_totaux_ligne_inter;
                    end

                    tirs_totaux_ligne_global = sum(data.RasterRace(cellules_du_cluster, :, :), 'all');
                    if tirs_totaux_ligne_global > 0
                        cellspecglobal_array(c) = tirs_in_box / tirs_totaux_ligne_global;
                    end
                end
            end
            c_cellspec_inter  = mean(cellspecinter_array, 'omitnan');
            c_cellspec_global = mean(cellspecglobal_array, 'omitnan');
            c_sce_purity      = mean(scepurity_array, 'omitnan');
        end
    end

    % 4. Lecture du Brainbow : Success Rate + Ncell (NOUVEAU)
    chemin_brainbow = fullfile(tous_les_fichiers(i).folder, 'brainbow.mat');
    c_success_rate = NaN;
    c_ncell        = NaN;
    if isfile(chemin_brainbow)
        try
            donnees_brainbow = load(chemin_brainbow, 'successrate', 'NCell');
            if isfield(donnees_brainbow, 'successrate')
                c_success_rate = donnees_brainbow.successrate;
            end
            if isfield(donnees_brainbow, 'NCell')
                c_ncell = donnees_brainbow.NCell;
            end
        catch
        end
    end

    c_pct_network = NaN;
    if ~isnan(c_mean_cell_assembly) && ~isnan(c_ncell) && c_ncell > 0
        c_pct_network = c_mean_cell_assembly / c_ncell * 100;
    end

    % --- Enregistrement dans les tableaux temporaires ---
    if ~isnan(c_cellspec_inter)
        temp_sessionName{i}     = sessName;
        temp_safeSessionName{i} = matlab.lang.makeValidName(sessName);
        temp_k(i)               = ck;
        temp_silhouette(i)      = c_silhouette;
        temp_sce_recall(i)      = c_sce_recall;
        temp_cellspec_inter(i)  = c_cellspec_inter;
        temp_cellspec_global(i) = c_cellspec_global;
        temp_sce_purity(i)      = c_sce_purity;
        temp_success_rate(i)    = c_success_rate;
        temp_mean_cell_assembly(i) = c_mean_cell_assembly;
        temp_ncell(i)           = c_ncell;
        temp_pct_network(i)     = c_pct_network;
        temp_paths{i}           = tous_les_fichiers(i).folder;
    end
    send(D, i);
end
close(hwait);

%% 5. REGROUPEMENT DES RESULTATS DANS LA STRUCTURE SessionsData
SessionsData = struct();
for i = 1:Nfichiers
    safeName = temp_safeSessionName{i};
    if isempty(safeName), continue; end
    if ~isfield(SessionsData, safeName)
        SessionsData.(safeName).OriginalName      = temp_sessionName{i};
        SessionsData.(safeName).K                 = [];
        SessionsData.(safeName).Silh              = [];
        SessionsData.(safeName).Recall            = [];
        SessionsData.(safeName).CellSpecInter     = [];
        SessionsData.(safeName).CellSpecGlobal    = [];
        SessionsData.(safeName).SCEPurity         = [];
        SessionsData.(safeName).SuccessRate       = [];
        SessionsData.(safeName).MeanCellAssembly  = [];
        SessionsData.(safeName).Ncell             = [];
        SessionsData.(safeName).PctNetwork        = [];
        SessionsData.(safeName).Paths             = {};
        SessionsData.(safeName).Scores            = [];
    end
    SessionsData.(safeName).K(end+1)                = temp_k(i);
    SessionsData.(safeName).Silh(end+1)             = temp_silhouette(i);
    SessionsData.(safeName).Recall(end+1)           = temp_sce_recall(i);
    SessionsData.(safeName).CellSpecInter(end+1)    = temp_cellspec_inter(i);
    SessionsData.(safeName).CellSpecGlobal(end+1)   = temp_cellspec_global(i);
    SessionsData.(safeName).SCEPurity(end+1)        = temp_sce_purity(i);
    SessionsData.(safeName).SuccessRate(end+1)      = temp_success_rate(i);
    SessionsData.(safeName).MeanCellAssembly(end+1) = temp_mean_cell_assembly(i);
    SessionsData.(safeName).Ncell(end+1)            = temp_ncell(i);
    SessionsData.(safeName).PctNetwork(end+1)       = temp_pct_network(i);
    SessionsData.(safeName).Paths{end+1}            = temp_paths{i};
end

%% 6. CALCUL DES SCORES ET PREPARATION DE LA TABLE EXCEL
sessionFields = fieldnames(SessionsData);
fprintf('--------------------------------------------------\n');

Excel_Session = {}; Excel_Rang = {}; Excel_K = {}; Excel_Score = {};
Excel_Silh = {}; Excel_Recall = {}; Excel_CellSpecInter = {}; Excel_CellSpecGlobal = {};
Excel_SCEPurity = {}; Excel_SuccessRate = {}; Excel_MeanCellAssembly = {}; Excel_PctNetwork = {};
Excel_SizeScore = {}; Excel_Path = {};

for s = 1:length(sessionFields)
    sess = SessionsData.(sessionFields{s});
    if length(sess.K) < 2, continue; end

    % 1. Silhouette : de [-1,1] vers [0,1]
    norm_sil = (sess.Silh + 1) / 2;
    % 2. Recall : de [0,100] vers [0,1]
    norm_rec = sess.Recall / 100;
    % 3. Specificité et Pureté : déjà des ratios entre 0 et 1
    norm_cellspec = sess.CellSpecInter;
    norm_purity   = sess.SCEPurity;
    % 4. Taille d'assemblée (NOUVEAU) : score en cloche autour de la cible %
    norm_size = arrayfun(@(p) sizeScore(p, size_pct_target_min, size_pct_target_max, ...
                                         size_pct_floor, size_pct_ceiling), sess.PctNetwork);

    % Sécurité : aucun score ne dépasse 1 ou ne descend sous 0
    norm_sil      = max(0, min(1, norm_sil));
    norm_rec      = max(0, min(1, norm_rec));
    norm_cellspec = max(0, min(1, norm_cellspec));
    norm_purity   = max(0, min(1, norm_purity));
    norm_size     = max(0, min(1, norm_size));

    % --- Calcul du Score Global ---
    Scores = weight_silhouette       * norm_sil ...
           + weight_recall           * norm_rec ...
           + weight_cell_specificity * norm_cellspec ...
           + weight_sce_purity       * norm_purity ...
           + weight_size             * norm_size;

    SessionsData.(sessionFields{s}).Scores = Scores;

    [sortedScores, sortIdx] = sort(Scores, 'descend');
    topk = min(nb_top, length(sortIdx));

    fprintf('SESSION %s\n', sess.OriginalName);
    for i = 1:topk
        origidx = sortIdx(i);
        fprintf('  %d - K %d\n', i, sess.K(origidx));
        fprintf('     Chemin %s\n', sess.Paths{origidx});
        fprintf('     Score %.2f | Silh %.3f | Recall %.1f | CellSpec %.3f | Taille %.1f cellules (%.1f%% réseau) | Succes %.2f\n', ...
            sortedScores(i), sess.Silh(origidx), sess.Recall(origidx), sess.CellSpecInter(origidx), ...
            sess.MeanCellAssembly(origidx), sess.PctNetwork(origidx), sess.SuccessRate(origidx));

        Excel_Session{end+1,1}          = sess.OriginalName;
        Excel_Rang{end+1,1}             = i;
        Excel_K{end+1,1}                = sess.K(origidx);
        Excel_Score{end+1,1}            = sortedScores(i);
        Excel_Silh{end+1,1}             = sess.Silh(origidx);
        Excel_Recall{end+1,1}           = sess.Recall(origidx);
        Excel_CellSpecInter{end+1,1}    = sess.CellSpecInter(origidx);
        Excel_CellSpecGlobal{end+1,1}   = sess.CellSpecGlobal(origidx);
        Excel_SCEPurity{end+1,1}        = sess.SCEPurity(origidx);
        Excel_SuccessRate{end+1,1}      = sess.SuccessRate(origidx);
        Excel_MeanCellAssembly{end+1,1} = sess.MeanCellAssembly(origidx);
        Excel_PctNetwork{end+1,1}       = sess.PctNetwork(origidx);
        Excel_SizeScore{end+1,1}        = norm_size(origidx);
        Excel_Path{end+1,1}             = sess.Paths{origidx};
    end
    fprintf('--------------------------------------------------\n');
end

%% 7. SAUVEGARDE DANS LE FICHIER EXCEL analysis_brainbow_best.xlsx
fprintf('Génération du fichier analysis_brainbow_best.xlsx...\n');
path_global = fullfile(dossier_mere, 'analysis_brainbow_best.xlsx');

header_global = {'Session','Rang','K','ScoreGlobal','Silhouette','Recall', ...
    'CellSpecInter','CellSpecGlobal','SCEPurity','SuccessRate', ...
    'MeanCellAssembly','PctReseau','SizeScore','CheminDossier'};

dataOut = [Excel_Session, Excel_Rang, Excel_K, Excel_Score, Excel_Silh, Excel_Recall, ...
    Excel_CellSpecInter, Excel_CellSpecGlobal, Excel_SCEPurity, Excel_SuccessRate, ...
    Excel_MeanCellAssembly, Excel_PctNetwork, Excel_SizeScore, Excel_Path];

writecell([header_global; dataOut], path_global, 'Range', 'A1');
fprintf('%d lignes (top %d par session) écrites dans %s\n', size(dataOut,1), nb_top, path_global);

%% 8. GRAPHIQUES : Taille d'assemblée / Score Global vs K
figure('Name','Analyse des métriques par K','Color','w','Position',[100 100 1300 650]);
colors = lines(length(sessionFields));

% --- Graphique 1 : Taille moyenne d'assemblée (% réseau) vs K ---
subplot(1,2,1); hold on;
for s = 1:length(sessionFields)
    sess = SessionsData.(sessionFields{s});
    if length(sess.K) < 2, continue; end
    [Ksorted, idxsort] = sort(sess.K);
    plot(Ksorted, sess.PctNetwork(idxsort), '-o', 'Color', colors(s,:), 'LineWidth', 1.5, ...
        'MarkerSize', 5, 'DisplayName', sess.OriginalName);
end
yline(size_pct_target_min, '--g', 'Cible min'); yline(size_pct_target_max, '--g', 'Cible max');
xlabel('Nombre de clusters K','FontWeight','bold');
ylabel('Taille assemblée (% du réseau)','FontWeight','bold');
title('Taille relative des assemblées vs K'); grid on;
hold off;

% --- Graphique 2 : Score Global vs K (moyenne toutes sessions) ---
subplot(1,2,2); hold on;
allK = []; allScores = [];
for s = 1:length(sessionFields)
    sess = SessionsData.(sessionFields{s});
    if isfield(sess,'Scores') && length(sess.K) >= 2
        [Ksorted, idxsort] = sort(sess.K);
        Scoressorted = sess.Scores(idxsort);
        allK = [allK, Ksorted]; allScores = [allScores, Scoressorted];
        plot(Ksorted, Scoressorted, '-s', 'Color', [colors(s,:) 0.4], 'LineWidth', 1, 'MarkerSize',4, ...
            'HandleVisibility','off');
    end
end
if ~isempty(allK)
    uniqueK = unique(allK);
    meanScores = arrayfun(@(k) mean(allScores(allK==k), 'omitnan'), uniqueK);
    plot(uniqueK, meanScores, '-k', 'LineWidth', 3, 'MarkerSize', 8, 'DisplayName','Moyenne Globale');
end
xlabel('Nombre de clusters K','FontWeight','bold');
ylabel('Score Global','FontWeight','bold');
title('Score Global (avec critère de taille) vs K'); grid on; legend('Location','best');
hold off;

nomgraph = fullfile(dossier_mere, 'AnalysePerformanceK_taille.png');
saveas(gcf, nomgraph);
fprintf('Graphique sauvegardé : %s\n', nomgraph);


%% =========================================================================
% FONCTIONS LOCALES (doivent être en fin de fichier)
% =========================================================================

function score = sizeScore(pct_network, target_min, target_max, floor_pct, ceiling_pct)
% Score en plateau/cloche : vaut 1 dans [target_min ; target_max],
% décroît linéairement vers 0 en dessous de floor_pct ou au-dessus de ceiling_pct.
    if isnan(pct_network)
        score = NaN;
    elseif pct_network >= target_min && pct_network <= target_max
        score = 1.0;
    elseif pct_network < target_min
        score = max(0, (pct_network - floor_pct) / (target_min - floor_pct));
    else
        score = max(0, (ceiling_pct - pct_network) / (ceiling_pct - target_max));
    end
end

function updateWaitbar(h, total)
    persistent count
    if isempty(count), count = 0; end
    count = count + 1;
    if isgraphics(h)
        waitbar(count/total, h, sprintf('Analyse %d/%d', count, total));
    end
    if count >= total, count = []; end
end
