%% ========================================================================
% SCRIPT : analyse_assembly_results5parfor_taille.m
% CHOIX DU MEILLEUR K PAR SESSION (règle aveugle aux couleurs)
%
% Pour chaque session, chaque K est noté sur des critères de QUALITÉ des
% assemblées, calculés uniquement à partir de l'activité. Le K qui a le
% meilleur score global est retenu. Aucune information de couleur
% (successrate, brainbow) n'entre dans le score : elle est seulement
% exportée pour information.
%
% Les poids donnent un tiers de l'influence à chacune des trois dimensions :
%   1. séparation des clusters de SCE        -> silhouette
%   2. couverture des SCE                    -> recall
%   3. concentration de l'activité           -> spécificité cellulaire + pureté
%      des assemblées                           des SCE (même numérateur,
%                                               donc un sixième chacun)
%
% Sorties (dans dossier_mere) :
%   analysis_brainbow_best.xlsx      : métriques du meilleur K par session
%   AnalysePerformanceK_taille.png   : taille et score en fonction de K
%   best_k_par_session.mat           : table best_k (chemins, identifiants)
%                                      + regle_choix_K (la règle utilisée)
%
% RÈGLE À FIGER : une fois les paramètres des sections 1 et 2 choisis, les
% noter avec la date et ne plus les modifier en fonction des résultats
% clonaux. Les ajuster d'après le diagnostic de la section 6b (aveugle aux
% couleurs) est en revanche légitime.
% =========================================================================

clear; clc; close all;

%% 1. PARAMETRES ET CHEMINS
dossier_mere  = 'E:\Data\Aurelie\analysis\Final\assembly\nocues\';
nb_top        = 1;       % nb de meilleurs K exportés par session dans l'Excel
k_start       = 9;       % plus petit K considéré (à justifier sans les couleurs)
date_fixation = '2026-10-27';   % date à laquelle cette règle a été figée

%% 2. POIDS DES CRITERES POUR LE SCORE GLOBAL
% Un tiers par dimension (voir en-tête). La somme vaut 1.
weight_silhouette       = 1/3;   % séparation des clusters de SCE
weight_recall           = 1/3;   % proportion de SCE expliqués (NRaceOK / NRace)
weight_cell_specificity = 1/6;   % part des activations d'une cellule qui tombent dans
                                 % les SCE de son propre cluster (parmi les SCE valides)
weight_sce_purity       = 1/6;   % part de l'activité d'un SCE portée par les cellules
                                 % du cluster qui l'explique
weight_size             = 0;     % taille d'assemblée : désactivée. C'est un a priori sur
                                 % la réponse (cible 10-15 % du réseau), pas un critère
                                 % de qualité. Calculée et exportée pour information.

% Mise à l'échelle des critères avant pondération :
%   'session' : dans chaque session, 0 = pire K, 1 = meilleur K, pour chaque
%               critère. Chaque critère a alors réellement l'influence de son
%               poids. (Contrepartie : un critère qui varie très peu reçoit
%               quand même toute son influence.)
%   'brute'   : échelle d'origine ((silhouette+1)/2, recall/100, ratios). Un
%               critère qui varie peu d'un K à l'autre pèse alors peu, quel
%               que soit son poids.
normalisation = 'session';

% Cible de taille (utilisée seulement si weight_size > 0), en % du réseau
size_pct_target_min = 10;   % note = 1 entre target_min et target_max
size_pct_target_max = 15;
size_pct_floor      = 4;    % note = 0 en dessous (décroissance linéaire)
size_pct_ceiling    = 30;   % note = 0 au dessus (décroissance linéaire)

%% 3. RECHERCHE DE TOUS LES FICHIERS RESULTS.MAT
fprintf('Recherche des fichiers results.mat...\n');
tous_les_fichiers = dir(fullfile(dossier_mere, '**', 'results.mat'));
if isempty(tous_les_fichiers)
    error('Aucun fichier results.mat trouvé.');
end

% Optionnel : restreindre à une liste de sessions
% sessions_autorisees = {'411582_230320_plane0','411582_230320_plane1','443669_220614_plane0','443669_220614_plane1', ...
%    '444113_220927_plane1','444119_220919_plane0','444119_220919_plane1','444119_220923_plane0','444119_220923_plane1',...
%    '444152_221128_plane1','444152_221129_plane0','444152_221129_plane1', '444175_221121_plane0','444175_221121_plane1','444175_221125_plane0','444175_221125_plane1', ...
%    '446178_1062025_plane0','446178_1062025_plane1','446178_11062025_plane0','446178_11062025_plane1','446198_2062025_plane0'};
% 
% tous_les_fichiers = tous_les_fichiers(contains({tous_les_fichiers.folder}, sessions_autorisees));
fprintf('%d fichiers results.mat trouvés.\n', length(tous_les_fichiers));

%% 4. EXTRACTION ET CALCULS DES METRIQUES (PARALLELISE)
fprintf('Démarrage du pool parallèle...\n');
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool;
end
Nfichiers = length(tous_les_fichiers);

temp_sessionName        = cell(Nfichiers, 1);
temp_safeSessionName    = cell(Nfichiers, 1);
temp_k                  = nan(Nfichiers, 1);
temp_silhouette         = nan(Nfichiers, 1);
temp_sce_recall         = nan(Nfichiers, 1);
temp_cellspec_inter     = nan(Nfichiers, 1);
temp_cellspec_global    = nan(Nfichiers, 1);
temp_sce_purity         = nan(Nfichiers, 1);
temp_success_rate       = nan(Nfichiers, 1);
temp_mean_cell_assembly = nan(Nfichiers, 1);
temp_ncell              = nan(Nfichiers, 1);
temp_pct_network        = nan(Nfichiers, 1);
temp_paths              = cell(Nfichiers, 1);

hwait = waitbar(0, 'Extraction des données (Parallel)...', 'Name', 'Progression');
D = parallel.pool.DataQueue;
afterEach(D, @(~) updateWaitbar(hwait, Nfichiers));

parfor i = 1:Nfichiers
    matFilePath = fullfile(tous_les_fichiers(i).folder, tous_les_fichiers(i).name);
    parts = strsplit(tous_les_fichiers(i).folder, filesep);
    if length(parts) < 2, send(D,i); continue; end
    kFolderName = parts{end};
    sessName    = parts{end-1};
    knumstr = regexp(kFolderName, '\d+', 'match');
    if isempty(knumstr), send(D,i); continue; end
    ck = str2double(knumstr{1});
    if ck < k_start, send(D,i); continue; end

    try
        data = load(matFilePath, 'RasterRace', 'NRace', 'NRaceOK', 'mean_sClOK', 'assemblyortho');
    catch
        send(D,i); continue;
    end

    % 1. Silhouette
    if isfield(data,'mean_sClOK') && isnumeric(data.mean_sClOK) && ~isnan(data.mean_sClOK)
        c_silhouette = data.mean_sClOK;
    else
        send(D,i); continue;
    end

    % 2. Recall
    if isfield(data,'NRaceOK') && isfield(data,'NRace') && data.NRace > 0
        c_sce_recall = data.NRaceOK / data.NRace * 100;
    else
        c_sce_recall = 0;
    end

    % 3. Spécificité / pureté + taille moyenne d'assemblée
    c_cellspec_inter     = NaN;
    c_cellspec_global    = NaN;
    c_sce_purity         = NaN;
    c_mean_cell_assembly = NaN;

    if isfield(data,'assemblyortho') && ~isempty(data.assemblyortho) && isfield(data,'RasterRace')
        nb_clusters = length(data.assemblyortho);
        nb_SCE = size(data.RasterRace, 2);
        c_mean_cell_assembly = mean(cellfun(@length, data.assemblyortho), 'omitnan');

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

    % 4. Brainbow : successrate (INFORMATION SEULEMENT, jamais dans le score) + NCell
    chemin_brainbow = fullfile(tous_les_fichiers(i).folder, 'brainbow.mat');
    c_success_rate = NaN;
    c_ncell        = NaN;
    if isfile(chemin_brainbow)
        try
            donnees_brainbow = load(chemin_brainbow, 'successrate', 'NCell');
            if isfield(donnees_brainbow, 'successrate'), c_success_rate = donnees_brainbow.successrate; end
            if isfield(donnees_brainbow, 'NCell'),       c_ncell        = donnees_brainbow.NCell;       end
        catch
        end
    end

    c_pct_network = NaN;
    if ~isnan(c_mean_cell_assembly) && ~isnan(c_ncell) && c_ncell > 0
        c_pct_network = c_mean_cell_assembly / c_ncell * 100;
    end

    if ~isnan(c_cellspec_inter)
        temp_sessionName{i}        = sessName;
        temp_safeSessionName{i}    = matlab.lang.makeValidName(sessName);
        temp_k(i)                  = ck;
        temp_silhouette(i)         = c_silhouette;
        temp_sce_recall(i)         = c_sce_recall;
        temp_cellspec_inter(i)     = c_cellspec_inter;
        temp_cellspec_global(i)    = c_cellspec_global;
        temp_sce_purity(i)         = c_sce_purity;
        temp_success_rate(i)       = c_success_rate;
        temp_mean_cell_assembly(i) = c_mean_cell_assembly;
        temp_ncell(i)              = c_ncell;
        temp_pct_network(i)        = c_pct_network;
        temp_paths{i}              = tous_les_fichiers(i).folder;
    end
    send(D, i);
end
close(hwait);

%% 5. REGROUPEMENT DES RESULTATS PAR SESSION
SessionsData = struct();
for i = 1:Nfichiers
    safeName = temp_safeSessionName{i};
    if isempty(safeName), continue; end
    if ~isfield(SessionsData, safeName)
        SessionsData.(safeName) = struct('OriginalName', temp_sessionName{i}, ...
            'K', [], 'Silh', [], 'Recall', [], 'CellSpecInter', [], 'CellSpecGlobal', [], ...
            'SCEPurity', [], 'SuccessRate', [], 'MeanCellAssembly', [], 'Ncell', [], ...
            'PctNetwork', [], 'Paths', {{}}, 'Scores', [], 'BestK', NaN);
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

%% 6. CALCUL DES SCORES ET TABLE EXCEL
sessionFields = fieldnames(SessionsData);
w = [weight_silhouette, weight_recall, weight_cell_specificity, weight_sce_purity];
noms_crit = {'silhouette', 'recall', 'specificite', 'purete'};
fprintf('--------------------------------------------------\n');

Excel_Session = {}; Excel_Rang = {}; Excel_K = {}; Excel_Score = {};
Excel_Silh = {}; Excel_Recall = {}; Excel_CellSpecInter = {}; Excel_CellSpecGlobal = {};
Excel_SCEPurity = {}; Excel_SuccessRate = {}; Excel_MeanCellAssembly = {}; Excel_PctNetwork = {};
Excel_SizeScore = {}; Excel_Path = {};

% Pour le diagnostic de la section 6b
influence = zeros(0, 4);   % par session : poids x étendue du critère entre les K
borne     = {};            % par session : 'basse', 'haute' ou 'interieur'
choix_par_critere = zeros(0, 4);   % 1 si le critère seul aurait choisi le même K

for s = 1:length(sessionFields)
    sess = SessionsData.(sessionFields{s});
    if length(sess.K) < 2, continue; end

    % Critères bruts, une ligne par K (plus grand = meilleur pour les quatre)
    crit = [sess.Silh(:), sess.Recall(:), sess.CellSpecInter(:), sess.SCEPurity(:)];

    switch normalisation
        case 'session'
            lo = min(crit, [], 1, 'omitnan');
            hi = max(crit, [], 1, 'omitnan');
            norm_c = (crit - lo) ./ max(hi - lo, eps);
        case 'brute'
            norm_c = [(crit(:,1) + 1) / 2, crit(:,2) / 100, crit(:,3), crit(:,4)];
        otherwise
            error('normalisation doit valoir ''session'' ou ''brute''.');
    end
    norm_c = max(0, min(1, norm_c));

    norm_size = arrayfun(@(p) sizeScore(p, size_pct_target_min, size_pct_target_max, ...
        size_pct_floor, size_pct_ceiling), sess.PctNetwork(:));

    % Score global (NaN si un critère manque pour ce K)
    Scores = norm_c * w';
    if weight_size > 0          % évite 0 x NaN = NaN quand NCell manque
        Scores = Scores + weight_size * norm_size;
    end
    Scores = Scores';
    SessionsData.(sessionFields{s}).Scores = Scores;

    % Tri : les K sans score (NaN) passent en dernier
    [sortedScores, sortIdx] = sort(Scores, 'descend', 'MissingPlacement', 'last');
    if isnan(sortedScores(1))
        warning('Session %s : aucun K avec un score complet, ignorée.', sess.OriginalName);
        continue
    end
    bestK = sess.K(sortIdx(1));
    SessionsData.(sessionFields{s}).BestK = bestK;

    % --- Diagnostic ---
    influence(end+1, :) = w .* (max(norm_c, [], 1, 'omitnan') - min(norm_c, [], 1, 'omitnan')); %#ok<SAGROW>
    if bestK == min(sess.K)
        borne{end+1} = 'basse'; %#ok<SAGROW>
    elseif bestK == max(sess.K)
        borne{end+1} = 'haute'; %#ok<SAGROW>
    else
        borne{end+1} = 'interieur'; %#ok<SAGROW>
    end
    [~, imax] = max(norm_c, [], 1);            % K choisi par chaque critère seul
    choix_par_critere(end+1, :) = sess.K(imax) == bestK; %#ok<SAGROW>

    % --- Affichage et Excel ---
    topk = min(nb_top, sum(~isnan(Scores)));
    fprintf('SESSION %s\n', sess.OriginalName);
    for i = 1:topk
        origidx = sortIdx(i);
        fprintf('  %d - K %d | Score %.2f | Silh %.3f | Recall %.1f | CellSpec %.3f | Pureté %.3f | Taille %.1f cell. (%.1f%% réseau)\n', ...
            i, sess.K(origidx), sortedScores(i), sess.Silh(origidx), sess.Recall(origidx), ...
            sess.CellSpecInter(origidx), sess.SCEPurity(origidx), ...
            sess.MeanCellAssembly(origidx), sess.PctNetwork(origidx));

        Excel_Session{end+1,1}          = sess.OriginalName; %#ok<SAGROW>
        Excel_Rang{end+1,1}             = i; %#ok<SAGROW>
        Excel_K{end+1,1}                = sess.K(origidx); %#ok<SAGROW>
        Excel_Score{end+1,1}            = sortedScores(i); %#ok<SAGROW>
        Excel_Silh{end+1,1}             = sess.Silh(origidx); %#ok<SAGROW>
        Excel_Recall{end+1,1}           = sess.Recall(origidx); %#ok<SAGROW>
        Excel_CellSpecInter{end+1,1}    = sess.CellSpecInter(origidx); %#ok<SAGROW>
        Excel_CellSpecGlobal{end+1,1}   = sess.CellSpecGlobal(origidx); %#ok<SAGROW>
        Excel_SCEPurity{end+1,1}        = sess.SCEPurity(origidx); %#ok<SAGROW>
        Excel_SuccessRate{end+1,1}      = sess.SuccessRate(origidx); %#ok<SAGROW>
        Excel_MeanCellAssembly{end+1,1} = sess.MeanCellAssembly(origidx); %#ok<SAGROW>
        Excel_PctNetwork{end+1,1}       = sess.PctNetwork(origidx); %#ok<SAGROW>
        Excel_SizeScore{end+1,1}        = norm_size(origidx); %#ok<SAGROW>
        Excel_Path{end+1,1}             = sess.Paths{origidx}; %#ok<SAGROW>
    end
end
fprintf('--------------------------------------------------\n');

%% 6b. DIAGNOSTIC DE LA RÈGLE (aveugle aux couleurs)
% Si le meilleur K tombe souvent sur une borne, un critère qui varie
% toujours dans le même sens avec K (ex. silhouette) décide à la place des
% données, et la borne (k_start) devient le vrai choix.
n_sess = numel(borne);
fprintf('\n=== Diagnostic de la règle de choix de K (%d sessions) ===\n', n_sess);
fprintf('Normalisation : %s | poids : silh %.2f, recall %.2f, spécif. %.2f, pureté %.2f, taille %.2f\n', ...
    normalisation, w, weight_size);
fprintf('Meilleur K sur la borne basse (K = %d ou plus petit K dispo) : %d / %d\n', ...
    k_start, sum(strcmp(borne, 'basse')), n_sess);
fprintf('Meilleur K sur la borne haute (plus grand K disponible)      : %d / %d\n', ...
    sum(strcmp(borne, 'haute')), n_sess);
fprintf('Influence effective médiane (poids x étendue entre K), et accord avec le choix final :\n');
for c = 1:4
    fprintf('   %-12s influence %.3f | aurait choisi le même K : %d / %d sessions\n', ...
        noms_crit{c}, median(influence(:, c), 'omitnan'), sum(choix_par_critere(:, c)), n_sess);
end

%% 7. EXCEL
path_global = fullfile(dossier_mere, 'analysis_brainbow_best.xlsx');
header_global = {'Session','Rang','K','ScoreGlobal','Silhouette','Recall', ...
    'CellSpecInter','CellSpecGlobal','SCEPurity','SuccessRate', ...
    'MeanCellAssembly','PctReseau','SizeScore','CheminDossier'};
dataOut = [Excel_Session, Excel_Rang, Excel_K, Excel_Score, Excel_Silh, Excel_Recall, ...
    Excel_CellSpecInter, Excel_CellSpecGlobal, Excel_SCEPurity, Excel_SuccessRate, ...
    Excel_MeanCellAssembly, Excel_PctNetwork, Excel_SizeScore, Excel_Path];
writecell([header_global; dataOut], path_global, 'Range', 'A1');
fprintf('\n%d lignes (top %d par session) écrites dans %s\n', size(dataOut,1), nb_top, path_global);

%% 7b. SAUVEGARDE .MAT : meilleur K par session
if isempty(Excel_Rang)
    warning('Aucune session classée : best_k_par_session.mat non créé.');
else
    estPremier = cell2mat(Excel_Rang) == 1;

    best_k = table();
    best_k.session         = Excel_Session(estPremier);
    best_k.K               = cell2mat(Excel_K(estPremier));
    best_k.score           = cell2mat(Excel_Score(estPremier));
    best_k.pct_reseau      = cell2mat(Excel_PctNetwork(estPremier));
    best_k.dossier         = Excel_Path(estPremier);
    best_k.chemin_results  = fullfile(best_k.dossier, 'results.mat');
    best_k.chemin_brainbow = fullfile(best_k.dossier, 'brainbow.mat');
    best_k.existe          = isfile(best_k.chemin_results);

    % Identifiants : 'P2M_444119_220919_plane0_26_09_...' -> '444119', '220919', 'plane0'
    tok = regexp(best_k.session, '(\d+)_(\d+)_(plane\d+)', 'tokens', 'once');
    if any(cellfun(@isempty, tok))
        error('Nom de session non reconnu : %s', ...
            strjoin(best_k.session(cellfun(@isempty, tok)), ', '));
    end
    best_k.animal = cellfun(@(t) t{1}, tok, 'UniformOutput', false);
    best_k.date   = cellfun(@(t) t{2}, tok, 'UniformOutput', false);
    best_k.plan   = cellfun(@(t) t{3}, tok, 'UniformOutput', false);

    % Même enregistrement analysé plusieurs fois -> avertissement
    cle = strcat(best_k.animal, '_', best_k.date, '_', best_k.plan);
    [cle_u, ~, ic_cle] = unique(cle);
    n_cle = accumarray(ic_cle, 1);
    for d = find(n_cle > 1)'
        warning('Plusieurs dossiers pour %s :\n   %s', cle_u{d}, ...
            strjoin(best_k.session(ic_cle == d), '\n   '));
    end

    liste_chemins = best_k.chemin_results;

    regle_choix_K = struct( ...
        'date_fixation',           date_fixation, ...
        'normalisation',           normalisation, ...
        'weight_silhouette',       weight_silhouette, ...
        'weight_recall',           weight_recall, ...
        'weight_cell_specificity', weight_cell_specificity, ...
        'weight_sce_purity',       weight_sce_purity, ...
        'weight_size',             weight_size, ...
        'size_pct_target_min',     size_pct_target_min, ...
        'size_pct_target_max',     size_pct_target_max, ...
        'size_pct_floor',          size_pct_floor, ...
        'size_pct_ceiling',        size_pct_ceiling, ...
        'k_start',                 k_start);

    path_mat = fullfile(dossier_mere, 'best_k_par_session.mat');
    save(path_mat, 'best_k', 'liste_chemins', 'regle_choix_K');
    fprintf('%d sessions écrites dans %s\n', height(best_k), path_mat);
    if any(~best_k.existe)
        warning('%d chemin(s) results.mat introuvable(s) : voir best_k.existe.', sum(~best_k.existe));
    end
end

%% 8. GRAPHIQUES : taille d'assemblée et score global en fonction de K
figure('Name','Analyse des métriques par K','Color','w','Position',[100 100 1300 650]);
colors = lines(length(sessionFields));

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

subplot(1,2,2); hold on;
allK = []; allScores = [];
for s = 1:length(sessionFields)
    sess = SessionsData.(sessionFields{s});
    if isempty(sess.Scores) || length(sess.K) < 2, continue; end
    [Ksorted, idxsort] = sort(sess.K);
    Scoressorted = sess.Scores(idxsort);
    allK = [allK, Ksorted]; allScores = [allScores, Scoressorted]; %#ok<AGROW>
    plot(Ksorted, Scoressorted, '-', 'Color', [colors(s,:) 0.4], 'LineWidth', 1, 'HandleVisibility','off');
    if ~isnan(sess.BestK)   % K retenu pour cette session
        plot(sess.BestK, sess.Scores(sess.K == sess.BestK), 'o', 'MarkerFaceColor', colors(s,:), ...
            'MarkerEdgeColor', 'none', 'MarkerSize', 6, 'HandleVisibility','off');
    end
end
if ~isempty(allK)
    uniqueK = unique(allK);
    meanScores = arrayfun(@(k) mean(allScores(allK==k), 'omitnan'), uniqueK);
    plot(uniqueK, meanScores, '-k', 'LineWidth', 3, 'DisplayName', 'Moyenne des sessions');
end
xlabel('Nombre de clusters K','FontWeight','bold');
ylabel(sprintf('Score global (normalisation : %s)', normalisation),'FontWeight','bold');
title('Score global vs K (points = K retenu)'); grid on; legend('Location','best');
hold off;

nomgraph = fullfile(dossier_mere, 'AnalysePerformanceK_taille.png');
saveas(gcf, nomgraph);
fprintf('Graphique sauvegardé : %s\n', nomgraph);


%% =========================================================================
% FONCTIONS LOCALES
% =========================================================================
function score = sizeScore(pct_network, target_min, target_max, floor_pct, ceiling_pct)
% Note de taille : 1 entre target_min et target_max, décroît linéairement
% jusqu'à 0 à floor_pct (en dessous) et à ceiling_pct (au dessus).
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
