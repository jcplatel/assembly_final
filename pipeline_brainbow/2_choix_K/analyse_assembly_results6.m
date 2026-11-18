function best_k = analyse_assembly_results6(sessions, config)
%% ========================================================================
% ÉTAPE 2 : CHOIX DU MEILLEUR K PAR SESSION (règle aveugle aux couleurs)
% = analyse_assembly_results5parfor_taille, adapté au nouveau pipeline.
% Les critères, les poids et la normalisation sont INCHANGÉS. Changements :
%  - la liste des sessions vient de la table ; chaque K est lu dans
%    <session>\kXX\assemblees.mat
%  - spécificité et pureté calculées sur Race (cellules x SCE) au lieu de
%    RasterRace (cellules x TRAMES : un SCE y était découpé sur 1 à 3 trames)
%  - plus aucune lecture de brainbow.mat : l'étape 2 ne voit AUCUNE couleur
%    (la colonne SuccessRate a disparu de l'Excel)
%  - le K retenu est écrit dans <session>\choix_K.mat, lu par l'étape 3
%  - une session sans K utilisable reçoit un statut explicite (« exclue : ... »)
%
% Sorties (dans config.PathSave) :
%   analysisALL.xlsx                 une ligne par session et par K (ancien exportdata)
%   analysis_brainbow_best.xlsx      métriques du meilleur K par session
%   AnalysePerformanceK_taille.png   taille et score en fonction de K
%   best_k_par_session.mat           table best_k + regle_choix_K
% =========================================================================

%% 1-2. PARAMETRES (config_brainbow.m, section choixK)
p = config.choixK;
dossier_mere            = config.PathSave;
nb_top                  = p.nb_top;
k_start                 = p.k_start;
weight_silhouette       = p.weight_silhouette;
weight_recall           = p.weight_recall;
weight_cell_specificity = p.weight_cell_specificity;
weight_sce_purity       = p.weight_sce_purity;
weight_size             = p.weight_size;
normalisation           = p.normalisation;
size_pct_target_min     = p.size_pct_target_min;
size_pct_target_max     = p.size_pct_target_max;
size_pct_floor          = p.size_pct_floor;
size_pct_ceiling        = p.size_pct_ceiling;
regle_choix_K           = p;
regle_choix_K.liste_K   = config.liste_K;
% K fixé (facultatif) : le même K pour toutes les sessions, au lieu du score
% (analyse de sensibilité : le résultat clonal dépend-il du K ?)
K_fixe = [];
if isfield(p, 'K_fixe') && ~isempty(p.K_fixe), K_fixe = p.K_fixe; end
if ~isempty(K_fixe)
    fprintf('K FIXÉ à %d pour toutes les sessions (config.choixK.K_fixe)\n', K_fixe);
end
% Règle du choix :
%   'exces' : K où la silhouette des clusters significatifs dépasse le plus
%             le hasard, au total : somme sur les clusters retenus de
%             (silhouette du cluster - seuil de son rang dans les surrogates).
%             Récompense à la fois le nombre de clusters réels et leur netteté.
%             Simulations (444119 simulé, 4 ou 6 assemblées implantées) :
%             assemblées retrouvées 2,2/4 et 4,8/6, contre 1,6/4 et 3,2/6 pour
%             'score', qui tombe presque toujours sur le plus petit K.
%   'score' : ancien score pondéré (silhouette, recall, spécificité, pureté).
%   'consensus' : pas de K unique : assemblées de consensus entre tous les K
%             de liste_K (>= k_start), voir consensus_session. Écrit
%             <session>\consensus\assemblees.mat ; choix_K.mat : best_K = 0.
regle_K = 'exces';
if isfield(p, 'regle') && ~isempty(p.regle), regle_K = p.regle; end
if ~ismember(regle_K, {'exces', 'score', 'consensus'})
    error('config.choixK.regle doit valoir ''exces'', ''score'' ou ''consensus''.');
end
if strcmp(regle_K, 'consensus') && ~isempty(K_fixe)
    error('config.choixK.K_fixe et la règle ''consensus'' ne vont pas ensemble.');
end
K_consensus = config.liste_K(config.liste_K >= k_start);
% réglages du consensus : ceux de config_brainbow, sinon les valeurs choisies sur simulation
defaut_cons = struct('normalisation', 'union', 'seuil', 0.7, 'min_K', 2);
if ~isfield(config.choixK, 'consensus'), config.choixK.consensus = struct(); end
for champ = fieldnames(defaut_cons)'
    if ~isfield(config.choixK.consensus, champ{1})
        config.choixK.consensus.(champ{1}) = defaut_cons.(champ{1});
    end
end
if strcmp(regle_K, 'consensus')
    pc = config.choixK.consensus;
    fprintf(['Consensus entre K = %s : co-appartenance ''%s'' >= %.2f, paires ensemble dans ' ...
        'au moins %d K\n'], mat2str(K_consensus), pc.normalisation, pc.seuil, pc.min_K);
end
n_cons = [];          % assemblées de consensus par session (affichage)
fprintf('Choix du K : règle ''%s''\n', regle_K);

w = [weight_silhouette, weight_recall, weight_cell_specificity, weight_sce_purity];
noms_crit = {'silhouette', 'recall', 'specificite', 'purete'};
nSess = height(sessions);

%% 3-5. MÉTRIQUES DE CHAQUE K, SESSION PAR SESSION
SessionsData = cell(nSess, 1);
lignes_analysisALL = {};          % une ligne par session et par K (ancien exportdata)
for s = 1:nSess
    sess = struct('OriginalName', sessions.nom{s}, 'K', [], 'Silh', [], 'Recall', [], ...
        'CellSpecInter', [], 'CellSpecGlobal', [], 'SCEPurity', [], 'MeanCellAssembly', [], ...
        'Ncell', [], 'PctNetwork', [], 'Exces', [], 'Paths', {{}}, 'Scores', [], 'BestK', NaN, 'statut', '');
    try
        pre = load(fullfile(sessions.dossier_source{s}, 'pretraitement.mat'), 'Race');
        Race = double(pre.Race);
        NCell = size(Race, 1);
        for ck = config.liste_K
            fichier_K = fichier_assemblees(sessions.dossier{s}, sessions.dossier_source{s}, ck);   % v1 ou v2
            if ~isfile(fichier_K), continue; end
            noms_var = intersect({'NRace', 'NRaceOK', 'mean_sClOK', 'assemblyortho', 'resume_K', ...
                'sCl', 'seuils_rang', 'NClOK'}, who('-file', fichier_K));
            data = load(fichier_K, noms_var{:});
            lignes_analysisALL{end+1, 1} = cell2table(struct2cell(data.resume_K)', ...
                'VariableNames', fieldnames(data.resume_K)'); %#ok<AGROW>
            if ck < k_start, continue; end

            % 1. Silhouette (NaN = aucune assemblée à ce K)
            if isnan(data.mean_sClOK), continue; end
            c_silhouette = data.mean_sClOK;
            % 2. Recall
            c_sce_recall = data.NRaceOK / data.NRace * 100;
            % 3. Spécificité / pureté + taille moyenne d'assemblée
            [c_cellspec_inter, c_cellspec_global, c_sce_purity] = specificite_purete(data.assemblyortho, Race);
            if isnan(c_cellspec_inter), continue; end
            c_mean_cell_assembly = mean(cellfun(@length, data.assemblyortho));

            sess.K(end+1)                = ck;
            sess.Silh(end+1)             = c_silhouette;
            sess.Recall(end+1)           = c_sce_recall;
            sess.CellSpecInter(end+1)    = c_cellspec_inter;
            sess.CellSpecGlobal(end+1)   = c_cellspec_global;
            sess.SCEPurity(end+1)        = c_sce_purity;
            sess.MeanCellAssembly(end+1) = c_mean_cell_assembly;
            sess.Ncell(end+1)            = NCell;
            sess.PctNetwork(end+1)       = c_mean_cell_assembly / NCell * 100;
            % silhouette en excès du hasard, sommée sur les clusters significatifs
            c_exces = NaN;
            if all(isfield(data, {'sCl', 'seuils_rang', 'NClOK'}))
                n_ok = min([data.NClOK, numel(data.sCl), numel(data.seuils_rang)]);
                c_exces = sum(data.sCl(1:n_ok) - data.seuils_rang(1:n_ok));
            end
            sess.Exces(end+1)            = c_exces;
            sess.Paths{end+1}            = fileparts(fichier_K);
        end
    catch err
        sess.statut = ['erreur : ' err.message];
    end
    SessionsData{s} = sess;
end

%% 6. CALCUL DES SCORES, choix_K.mat DE CHAQUE SESSION
fprintf('--------------------------------------------------\n');
Excel_Session = {}; Excel_Rang = {}; Excel_K = {}; Excel_Score = {};
Excel_Silh = {}; Excel_Recall = {}; Excel_CellSpecInter = {}; Excel_CellSpecGlobal = {};
Excel_SCEPurity = {}; Excel_MeanCellAssembly = {}; Excel_PctNetwork = {};
Excel_SizeScore = {}; Excel_Path = {};

influence = zeros(0, 4);           % par session : poids x étendue du critère entre les K
borne     = {};                    % 'basse', 'haute' ou 'interieur'
choix_par_critere = zeros(0, 4);   % 1 si le critère seul aurait choisi le même K
statuts = cell(nSess, 1);
best_K_sessions = nan(nSess, 1);
score_sessions  = nan(nSess, 1);
pct_sessions    = nan(nSess, 1);

for s = 1:nSess
    sess = SessionsData{s};
    if ~isfolder(sessions.dossier{s}), mkdir(sessions.dossier{s}); end
    fichier_choix = fullfile(sessions.dossier{s}, 'choix_K.mat');
    if isfile(fichier_choix), delete(fichier_choix); end    % jamais d'ancien choix
    fichier_cons = fullfile(sessions.dossier{s}, 'consensus', 'assemblees.mat');
    if isfile(fichier_cons), delete(fichier_cons); end      % ni d'ancien consensus
    best_K = NaN;
    tableau_K = table();

    if ~isempty(sess.statut)
        statut_choix = sess.statut;                           % erreur de lecture
    elseif isempty(sess.K) || (numel(sess.K) < 2 && isempty(K_fixe) && strcmp(regle_K, 'score'))
        statut_choix = sprintf('exclue : %d K avec des assemblées (il en faut 2)', numel(sess.K));
    else
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
        norm_size = arrayfun(@(pc) sizeScore(pc, size_pct_target_min, size_pct_target_max, ...
            size_pct_floor, size_pct_ceiling), sess.PctNetwork(:));

        % Score global (NaN si un critère manque pour ce K)
        Scores = norm_c * w';
        if weight_size > 0          % évite 0 x NaN = NaN
            Scores = Scores + weight_size * norm_size;
        end
        Scores = Scores';
        sess.Scores = Scores;
        [sortedScores, sortIdx] = sort(Scores, 'descend', 'MissingPlacement', 'last');
        if strcmp(regle_K, 'exces')       % à égalité, le plus petit K (sort est stable)
            [~, sortIdx] = sort(sess.Exces, 'descend', 'MissingPlacement', 'last');
            sortedScores = Scores(sortIdx);
        end

        tableau_K = table(sess.K(:), sess.Silh(:), sess.Recall(:), sess.CellSpecInter(:), ...
            sess.CellSpecGlobal(:), sess.SCEPurity(:), sess.MeanCellAssembly(:), ...
            sess.PctNetwork(:), Scores(:), sess.Exces(:), 'VariableNames', {'K', 'Silhouette', 'Recall', ...
            'CellSpecInter', 'CellSpecGlobal', 'SCEPurity', 'MeanCellAssembly', 'PctReseau', 'Score', 'Exces'});

        if ~isempty(K_fixe) && ~ismember(K_fixe, sess.K)
            statut_choix = sprintf('exclue : aucune assemblée au K fixé (K = %d)', K_fixe);
        elseif isempty(K_fixe) && strcmp(regle_K, 'score') && isnan(sortedScores(1))
            statut_choix = 'exclue : aucun K avec un score complet';
        elseif isempty(K_fixe) && strcmp(regle_K, 'exces') && all(isnan(sess.Exces))
            statut_choix = 'exclue : silhouette en excès indisponible (étape 1 à refaire)';
        elseif strcmp(regle_K, 'consensus')
            [statut_choix, rc] = consensus_session(sessions.nom{s}, sessions.dossier_source{s}, ...
                sessions.dossier{s}, K_consensus, config, sessions.graine(s));
            if strcmp(statut_choix, 'ok')
                best_K = 0;                       % 0 = consensus
                sess.BestK = 0;
                best_K_sessions(s) = 0;
                n_cons(end+1) = rc.n_asm; %#ok<AGROW>
                fprintf('SESSION %s\n  consensus de %d K : %d assemblée(s) [%s] ; vues à : %s\n', ...
                    sess.OriginalName, rc.n_K, rc.n_asm, rc.tailles, rc.K_presence);
                Excel_Session{end+1,1} = sess.OriginalName; %#ok<AGROW>
                Excel_Rang{end+1,1} = 1; %#ok<AGROW>
                Excel_K{end+1,1} = 0; %#ok<AGROW>
                Excel_Score{end+1,1} = NaN; %#ok<AGROW>
                Excel_Silh{end+1,1} = NaN; %#ok<AGROW>
                Excel_Recall{end+1,1} = NaN; %#ok<AGROW>
                Excel_CellSpecInter{end+1,1} = NaN; %#ok<AGROW>
                Excel_CellSpecGlobal{end+1,1} = NaN; %#ok<AGROW>
                Excel_SCEPurity{end+1,1} = NaN; %#ok<AGROW>
                tl = str2num(rc.tailles); %#ok<ST2NM>
                Excel_MeanCellAssembly{end+1,1} = mean(tl); %#ok<AGROW>
                Excel_PctNetwork{end+1,1} = mean(tl) / max(sess.Ncell(1), 1) * 100; %#ok<AGROW>
                Excel_SizeScore{end+1,1} = NaN; %#ok<AGROW>
                Excel_Path{end+1,1} = fileparts(fichier_cons); %#ok<AGROW>
            end
        else
            statut_choix = 'ok';
            if ~isempty(K_fixe)          % K imposé : il passe en tête
                i_fixe = find(sess.K == K_fixe, 1);
                sortIdx = [i_fixe, sortIdx(sortIdx ~= i_fixe)];
                sortedScores = Scores(sortIdx);
            end
            best_K = sess.K(sortIdx(1));
            sess.BestK = best_K;
            best_K_sessions(s) = best_K;
            score_sessions(s)  = sortedScores(1);
            pct_sessions(s)    = sess.PctNetwork(sortIdx(1));

            % --- Diagnostic ---
            influence(end+1, :) = w .* (max(norm_c, [], 1, 'omitnan') - min(norm_c, [], 1, 'omitnan')); %#ok<AGROW>
            if best_K == min(sess.K)
                borne{end+1} = 'basse'; %#ok<AGROW>
            elseif best_K == max(sess.K)
                borne{end+1} = 'haute'; %#ok<AGROW>
            else
                borne{end+1} = 'interieur'; %#ok<AGROW>
            end
            [~, imax] = max(norm_c, [], 1);
            choix_par_critere(end+1, :) = sess.K(imax) == best_K; %#ok<AGROW>

            % --- Affichage et Excel ---
            topk = min(nb_top, sum(~isnan(Scores)));
            fprintf('SESSION %s\n', sess.OriginalName);
            for i = 1:topk
                origidx = sortIdx(i);
                fprintf('  %d - K %d | Excès %.3f | Score %.2f | Silh %.3f | Recall %.1f | CellSpec %.3f | Pureté %.3f | Taille %.1f cell. (%.1f%% réseau)\n', ...
                    i, sess.K(origidx), sess.Exces(origidx), sortedScores(i), sess.Silh(origidx), sess.Recall(origidx), ...
                    sess.CellSpecInter(origidx), sess.SCEPurity(origidx), ...
                    sess.MeanCellAssembly(origidx), sess.PctNetwork(origidx));
                Excel_Session{end+1,1}          = sess.OriginalName; %#ok<AGROW>
                Excel_Rang{end+1,1}             = i; %#ok<AGROW>
                Excel_K{end+1,1}                = sess.K(origidx); %#ok<AGROW>
                Excel_Score{end+1,1}            = sortedScores(i); %#ok<AGROW>
                Excel_Silh{end+1,1}             = sess.Silh(origidx); %#ok<AGROW>
                Excel_Recall{end+1,1}           = sess.Recall(origidx); %#ok<AGROW>
                Excel_CellSpecInter{end+1,1}    = sess.CellSpecInter(origidx); %#ok<AGROW>
                Excel_CellSpecGlobal{end+1,1}   = sess.CellSpecGlobal(origidx); %#ok<AGROW>
                Excel_SCEPurity{end+1,1}        = sess.SCEPurity(origidx); %#ok<AGROW>
                Excel_MeanCellAssembly{end+1,1} = sess.MeanCellAssembly(origidx); %#ok<AGROW>
                Excel_PctNetwork{end+1,1}       = sess.PctNetwork(origidx); %#ok<AGROW>
                Excel_SizeScore{end+1,1}        = norm_size(origidx); %#ok<AGROW>
                Excel_Path{end+1,1}             = sess.Paths{origidx}; %#ok<AGROW>
            end
        end
    end
    statuts{s} = statut_choix;
    SessionsData{s} = sess;
    save(fichier_choix, 'best_K', 'statut_choix', 'tableau_K', 'regle_choix_K');
    if ~strcmp(statut_choix, 'ok')
        fprintf('SESSION %s : %s\n', sess.OriginalName, statut_choix);
    end
end
fprintf('--------------------------------------------------\n');

%% 6b. DIAGNOSTIC DE LA RÈGLE (aveugle aux couleurs)
% Si le meilleur K tombe souvent sur une borne, un critère qui varie toujours
% dans le même sens avec K décide à la place des données.
n_sess = numel(borne);
if strcmp(regle_K, 'consensus')
    fprintf('\n=== Consensus : %d sessions avec au moins une assemblée ; assemblées par session : médiane %g [%g - %g] ===\n', ...
        numel(n_cons), median(n_cons), min([n_cons, NaN]), max([n_cons, NaN]));
end
fprintf('\n=== Diagnostic de la règle de choix de K (%d sessions) ===\n', n_sess);
fprintf('Règle : %s', regle_K);
if strcmp(regle_K, 'exces')
    fprintf(' (silhouette en excès du hasard, sommée sur les clusters significatifs)\n');
else
    fprintf('\n');
end
fprintf('Score (affiché) : normalisation %s | poids : silh %.2f, recall %.2f, spécif. %.2f, pureté %.2f, taille %.2f\n', ...
    normalisation, w, weight_size);
fprintf('Meilleur K sur la borne basse (plus petit K disponible) : %d / %d\n', sum(strcmp(borne, 'basse')), n_sess);
fprintf('Meilleur K sur la borne haute (plus grand K disponible) : %d / %d\n', sum(strcmp(borne, 'haute')), n_sess);
fprintf('Influence effective médiane (poids x étendue entre K), et accord avec le choix final :\n');
for c = 1:4
    fprintf('   %-12s influence %.3f | aurait choisi le même K : %d / %d sessions\n', ...
        noms_crit{c}, median(influence(:, c), 'omitnan'), sum(choix_par_critere(:, c)), n_sess);
end

%% 7. EXCEL
path_global = fullfile(dossier_mere, 'analysis_brainbow_best.xlsx');
header_global = {'Session', 'Rang', 'K', 'ScoreGlobal', 'Silhouette', 'Recall', ...
    'CellSpecInter', 'CellSpecGlobal', 'SCEPurity', 'MeanCellAssembly', 'PctReseau', ...
    'SizeScore', 'CheminDossier'};
dataOut = [Excel_Session, Excel_Rang, Excel_K, Excel_Score, Excel_Silh, Excel_Recall, ...
    Excel_CellSpecInter, Excel_CellSpecGlobal, Excel_SCEPurity, ...
    Excel_MeanCellAssembly, Excel_PctNetwork, Excel_SizeScore, Excel_Path];
if isfile(path_global), delete(path_global); end
writecell([header_global; dataOut], path_global, 'Range', 'A1');

%% 7a. analysisALL.xlsx : toutes les sessions et tous les K (ancien exportdata),
% réécrit entièrement à chaque lancement (plus de doublons)
if ~isempty(lignes_analysisALL)
    path_all = fullfile(dossier_mere, 'analysisALL.xlsx');
    if isfile(path_all), delete(path_all); end
    writetable(vertcat(lignes_analysisALL{:}), path_all);
end

%% 7b. TABLE best_k (pour information ; l'étape 3 lit choix_K.mat)
best_k = table(sessions.nom, best_K_sessions, score_sessions, pct_sessions, statuts, sessions.dossier, ...
    'VariableNames', {'session', 'K', 'score', 'pct_reseau', 'statut', 'dossier'});
save(fullfile(dossier_mere, 'best_k_par_session.mat'), 'best_k', 'regle_choix_K');
fprintf('%d sessions avec un K retenu, %d exclues (voir best_k.statut)\n', ...
    sum(strcmp(statuts, 'ok')), sum(~strcmp(statuts, 'ok')));

%% 8. GRAPHIQUES : taille d'assemblée et score global en fonction de K
fig = figure('Name', 'Analyse des métriques par K', 'Color', 'w', 'Position', [100 100 1300 650], ...
    'Visible', 'off');
colors = lines(nSess);
subplot(1, 2, 1); hold on;
for s = 1:nSess
    sess = SessionsData{s};
    if numel(sess.K) < 2, continue; end
    [Ksorted, idxsort] = sort(sess.K);
    plot(Ksorted, sess.PctNetwork(idxsort), '-o', 'Color', colors(s,:), 'LineWidth', 1.5, ...
        'MarkerSize', 5, 'DisplayName', sess.OriginalName);
end
yline(size_pct_target_min, '--g', 'Cible min'); yline(size_pct_target_max, '--g', 'Cible max');
xlabel('Nombre de clusters K', 'FontWeight', 'bold');
ylabel('Taille assemblée (% du réseau)', 'FontWeight', 'bold');
title('Taille relative des assemblées vs K'); grid on; hold off;

subplot(1, 2, 2); hold on;
allK = []; allScores = [];
for s = 1:nSess
    sess = SessionsData{s};
    if isempty(sess.Scores), continue; end
    [Ksorted, idxsort] = sort(sess.K);
    Scoressorted = sess.Scores(idxsort);
    allK = [allK, Ksorted]; allScores = [allScores, Scoressorted]; %#ok<AGROW>
    plot(Ksorted, Scoressorted, '-', 'Color', [colors(s,:) 0.4], 'LineWidth', 1, 'HandleVisibility', 'off');
    if ~isnan(sess.BestK) && any(sess.K == sess.BestK)
        plot(sess.BestK, sess.Scores(sess.K == sess.BestK), 'o', 'MarkerFaceColor', colors(s,:), ...
            'MarkerEdgeColor', 'none', 'MarkerSize', 6, 'HandleVisibility', 'off');
    end
end
if ~isempty(allK)
    uniqueK = unique(allK);
    meanScores = arrayfun(@(k) mean(allScores(allK == k), 'omitnan'), uniqueK);
    plot(uniqueK, meanScores, '-k', 'LineWidth', 3, 'DisplayName', 'Moyenne des sessions');
end
xlabel('Nombre de clusters K', 'FontWeight', 'bold');
ylabel(sprintf('Score global (normalisation : %s)', normalisation), 'FontWeight', 'bold');
title('Score global vs K (points = K retenu)'); grid on; legend('Location', 'best'); hold off;
exportgraphics(fig, fullfile(dossier_mere, 'AnalysePerformanceK_taille.png'), 'Resolution', 150);
close(fig)
end


%% =========================================================================
% FONCTIONS LOCALES
% =========================================================================
function [c_cellspec_inter, c_cellspec_global, c_sce_purity] = specificite_purete(assemblyortho, Race)
% Même calcul que dans analyse_assembly_results5, mais sur Race (cellules x SCE).
% Chaque SCE est attribué à l'assemblée qui y a le plus de cellules actives.
c_cellspec_inter = NaN; c_cellspec_global = NaN; c_sce_purity = NaN;
nb_clusters = numel(assemblyortho);
nb_SCE = size(Race, 2);
if nb_clusters == 0 || nb_SCE == 0, return; end

participationSCE = zeros(nb_clusters, nb_SCE);
for c = 1:nb_clusters
    participationSCE(c, :) = sum(Race(assemblyortho{c}, :), 1);
end
[valmax, SCElabels] = max(participationSCE, [], 1);
SCElabels(valmax == 0) = 0;
idxSCE_valides = (SCElabels ~= 0);

cellspecinter_array  = nan(1, nb_clusters);
cellspecglobal_array = nan(1, nb_clusters);
scepurity_array      = nan(1, nb_clusters);
for c = 1:nb_clusters
    cellules_du_cluster = assemblyortho{c};
    idxSCE_du_cluster   = (SCElabels == c);
    if sum(idxSCE_du_cluster) > 0 && ~isempty(cellules_du_cluster)
        tirs_in_box = sum(Race(cellules_du_cluster, idxSCE_du_cluster), 'all');
        tirs_totaux_colonne = sum(Race(:, idxSCE_du_cluster), 'all');
        if tirs_totaux_colonne > 0
            scepurity_array(c) = tirs_in_box / tirs_totaux_colonne;
        end
        tirs_totaux_ligne_inter = sum(Race(cellules_du_cluster, idxSCE_valides), 'all');
        if tirs_totaux_ligne_inter > 0
            cellspecinter_array(c) = tirs_in_box / tirs_totaux_ligne_inter;
        end
        tirs_totaux_ligne_global = sum(Race(cellules_du_cluster, :), 'all');
        if tirs_totaux_ligne_global > 0
            cellspecglobal_array(c) = tirs_in_box / tirs_totaux_ligne_global;
        end
    end
end
c_cellspec_inter  = mean(cellspecinter_array, 'omitnan');
c_cellspec_global = mean(cellspecglobal_array, 'omitnan');
c_sce_purity      = mean(scepurity_array, 'omitnan');
end


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
