function main_brainbow(config)
%% =========================================================================
%  MAIN_BRAINBOW : assemblées (Malvache et al. 2016 modifié) + analyse clonale
%  -------------------------------------------------------------------------
%  ÉTAPE 1  fonctionnel  par session : prétraitement, puis assemblées pour
%                        chaque K. AUCUNE couleur n'est chargée.
%  ÉTAPE 2  choix du K   par session, sur l'activité seulement
%  ÉTAPE 3  clonal       par session, sur le K retenu : SEULE étape qui
%                        charge les couleurs
%  ÉTAPE 3b paires       (config.faire_etape3b) synchronie au repos des
%                        voisines de même couleur contre couleurs
%                        différentes, sans les assemblées (= main_paires),
%                        puis regroupement de toutes les sessions
%  ÉTAPE 4  agrégation   toutes les sessions, tous les animaux ; refuse de
%                        tourner tant qu'une session n'est pas finie
%  ÉTAPE 5  séquences    (config.faire_etape5) cellules activées en séquence
%                        pendant la course (méthode de sequence_2026_01_27),
%                        puis distance dans la séquence même couleur contre
%                        couleurs différentes (= main_sequences), puis
%                        regroupement de toutes les sessions
%
%  Lancer :   main_brainbow              (paramètres de config_brainbow.m)
%             main_brainbow(config)      (paramètres modifiés, ex. test)
%
%  Une seule session (ou un seul animal) :
%             config = config_brainbow();
%             config.seulement = {'444119_220919_plane0_sans_indice'};   % ou {'444119'}
%             main_brainbow(config)
%  Étapes 1 et 3 pour ces sessions seulement ; l'étape 2 est refaite pour
%  toutes les sessions déjà analysées (les fichiers de résumé restent
%  complets) ; l'étape 4 est sautée.
%
%  Dossiers (voir dossiers_analyse) :
%    config.PathSave            étape 1, un sous-dossier par session ; l'analyse
%                               principale (étapes 2 à 4) au même endroit
%    config.analyse = 'nom'     étapes 2 à 4 dans <PathSave>\analyses\nom, sur
%                               l'étape 1 de PathSave (rien n'est recopié)
%
%  Tout le travail se fait dans des fonctions : aucune variable ne peut
%  rester d'une session, d'un K ou d'un animal à l'autre.
%  =========================================================================

if nargin < 1, config = config_brainbow(); end
close all

% Les tirages du pipeline changent le générateur aléatoire de MATLAB :
% son état d'origine est remis à la fin, même en cas d'erreur.
etat_rng = rng;
remettre_rng = onCleanup(@() rng(etat_rng)); %#ok<NASGU>

% Les fonctions de ce dossier passent avant d'anciennes versions du même nom
dossier_pipeline = fileparts(mfilename('fullpath'));
addpath(genpath(dossier_pipeline));

%% Dossiers : étape 1 dans config.PathSave ; étapes 2 à 4 au même endroit
% (analyse principale) ou dans <PathSave>\analyses\<config.analyse>
% (autre analyse sur la même étape 1, ex. une autre règle de choix du K).
% Voir dossiers_analyse.
config = dossiers_analyse(config);

%% Sessions à traiter : une ligne par session dans la table
sessions = lire_table_sessions(config);
a_traiter = 1:height(sessions);
partiel = isfield(config, 'seulement') && ~isempty(config.seulement);
if partiel
    seulement = cellstr(config.seulement);
    a_traiter = find(ismember(sessions.nom, seulement) | ismember(sessions.id_animal, seulement))';
    if isempty(a_traiter)
        error('config.seulement: no session of the table (inclure = 1) matches %s.', ...
            strjoin(seulement, ', '));
    end
end
nSess = numel(a_traiter);
if ~isfolder(config.PathSave), mkdir(config.PathSave); end

log_file = fullfile(config.PathSave, sprintf('log_analyse_%s.txt', datestr(now, 'yy_mm_dd_HH_MM')));
historique = {};
historique = ecrire_log(log_file, historique, sprintf('=== START OF ANALYSIS: %s (%d sessions) ===', ...
    datestr(now), nSess));
historique = ecrire_log(log_file, historique, sprintf('Step 1 folder: %s', config.PathEtape1));
if ~strcmp(config.PathEtape1, config.PathSave)
    historique = ecrire_log(log_file, historique, sprintf('Analysis "%s" (steps 2-4): %s', ...
        config.analyse, config.PathSave));
end
temps_global = tic;

if isempty(gcp('nocreate')), parpool; end
echecs = cell(0, 3);   % {étape, session, erreur}

%% ========== ÉTAPE 1 : prétraitement + assemblées pour chaque K ==========
if config.faire_etape1
    for i_sess = 1:nSess
        session = table2struct(sessions(a_traiter(i_sess), :));
        session.dossier = session.dossier_source;   % l'étape 1 écrit dans le dossier de l'étape 1
        temps_fichier = tic;
        signaler(config, 0.88 * (i_sess - 1) / nSess, ...
            sprintf('Step 1: session %d/%d (%s), preprocessing', i_sess, nSess, session.nom));
        config_s = config;      % avancement K par K (interface)
        if isfield(config, 'rappel_progression') && ~isempty(config.rappel_progression)
            config_s.progression_K = @(iK, nK, K) signaler(config, 0.88 * ((i_sess - 1) + iK / nK) / nSess, ...
                sprintf('Step 1: session %d/%d (%s), K = %d (%d/%d)', i_sess, nSess, session.nom, K, iK + 1, nK));
        end
        try
            etape1_fonctionnel(session, config_s);
            msg = sprintf('Step 1 - %d/%d [%s] done (%s)', i_sess, nSess, session.nom, ...
                duree(toc(temps_fichier)));
            ecrire_analysisALL(sessions, config);    % Excel mis à jour après chaque session
        catch err
            echecs(end+1, :) = {'step 1', session.nom, err.message}; %#ok<AGROW>
            msg = sprintf('Step 1 - %d/%d [%s] FAILED: %s', i_sess, nSess, session.nom, err.message);
        end
        historique = ecrire_log(log_file, historique, msg);
        appeler_rappel(config, session.nom, 1, msg);   % interface : affichage, arrêt demandé
    end
end

%% ========== ÉTAPE 2 : choix du K de chaque session ======================
if config.faire_etape2
    % Avec config.seulement : toutes les sessions déjà analysées (dossier
    % présent), pour que les fichiers de résumé restent complets
    sessions_etape2 = sessions;
    if partiel
        sessions_etape2 = sessions(cellfun(@isfolder, sessions.dossier_source), :);
    end
    signaler(config, 0.88, 'Step 2: K choice');
    try
        analyse_assembly_results6(sessions_etape2, config);
        historique = ecrire_log(log_file, historique, 'Step 2 - K choice done');
    catch err
        echecs(end+1, :) = {'step 2', 'all', err.message};
        historique = ecrire_log(log_file, historique, ['Step 2 - FAILED: ' err.message]);
    end
end

%% ========== ÉTAPE 3 : statistiques clonales sur le K retenu =============
id_execution = '';    % identifie les résultats de CETTE exécution de l'étape 3
faire_3b = isfield(config, 'faire_etape3b') && config.faire_etape3b;
fin_3 = 0.98; if faire_3b, fin_3 = 0.94; end      % barre d'avancement
if config.faire_etape3
    id_execution = datestr(now, 'yy_mm_dd_HH_MM_SS');
    for i_sess = 1:nSess
        session = table2struct(sessions(a_traiter(i_sess), :));
        signaler(config, 0.90 + (fin_3 - 0.90) * (i_sess - 1) / nSess, ...
            sprintf('Step 3: session %d/%d (%s), clonal statistics', i_sess, nSess, session.nom));
        try
            statut = etape3_clonal(session, config, id_execution);
            msg = sprintf('Step 3 - %d/%d [%s]: %s', i_sess, nSess, session.nom, statut);
        catch err
            echecs(end+1, :) = {'step 3', session.nom, err.message}; %#ok<AGROW>
            msg = sprintf('Step 3 - %d/%d [%s] FAILED: %s', i_sess, nSess, session.nom, err.message);
        end
        historique = ecrire_log(log_file, historique, msg);
        appeler_rappel(config, session.nom, 3, msg);
    end
end

%% ========== ÉTAPE 3b : synchronie des paires selon la couleur ===========
if faire_3b
    for i_sess = 1:nSess
        session = table2struct(sessions(a_traiter(i_sess), :));
        signaler(config, fin_3 + (0.98 - fin_3) * (i_sess - 1) / nSess, ...
            sprintf('Step 3b: session %d/%d (%s), pair synchrony', i_sess, nSess, session.nom));
        t0 = tic;
        if ~isfile(fullfile(session.dossier_source, 'pretraitement.mat'))
            msg = sprintf('Step 3b - %d/%d [%s] skipped: no preprocessing (step 1)', i_sess, nSess, session.nom);
        else
            try
                etape3b_paires(session, config);
                msg = sprintf('Step 3b - %d/%d [%s] done (%s)', i_sess, nSess, session.nom, duree(toc(t0)));
            catch err
                echecs(end+1, :) = {'step 3b', session.nom, err.message}; %#ok<AGROW>
                msg = sprintf('Step 3b - %d/%d [%s] FAILED: %s', i_sess, nSess, session.nom, err.message);
            end
        end
        historique = ecrire_log(log_file, historique, msg);
        appeler_rappel(config, session.nom, 3, msg);
    end
    % Regroupement : toutes les sessions qui ont des résultats de paires
    signaler(config, 0.98, 'Step 3b: pairs, all sessions');
    try
        agregation_paires(sessions, config);
        historique = ecrire_log(log_file, historique, sprintf('Step 3b - pairs aggregated (%s)', ...
            fullfile(config.PathSave, 'paires_agregation')));
    catch err
        echecs(end+1, :) = {'step 3b', 'all', err.message};
        historique = ecrire_log(log_file, historique, ['Step 3b - aggregation FAILED: ' err.message]);
    end
end

%% ========== ÉTAPE 4 : agrégation ========================================
fichier_resume = '';     % phrases récapitulatives, affichées tout à la fin
lignes_resume  = {};
if config.faire_etape4 && partiel
    historique = ecrire_log(log_file, historique, ['Step 4 skipped (config.seulement): ' ...
        'run main_brainbow without config.seulement once all sessions are done.']);
elseif config.faire_etape4
    if ~isempty(echecs) && ~config.agreg.autoriser_partiel
        historique = ecrire_log(log_file, historique, ...
            'Step 4 NOT run: some sessions failed (list below). Fix them, then run again.');
    else
        signaler(config, 0.98, 'Step 4: aggregation');
        try
            [fichier_resume, lignes_resume] = agregation_stats_clonales3(sessions, config, id_execution);
            historique = ecrire_log(log_file, historique, 'Step 4 - aggregation done');
        catch err
            echecs(end+1, :) = {'step 4', 'all', err.message};
            historique = ecrire_log(log_file, historique, ['Step 4 - FAILED: ' err.message]);
        end
    end
end

%% ========== ÉTAPE 5 : séquences pendant la course, puis couleurs ========
if isfield(config, 'faire_etape5') && config.faire_etape5
    for i_sess = 1:nSess
        session = table2struct(sessions(a_traiter(i_sess), :));
        signaler(config, 0.99, sprintf('Step 5: session %d/%d (%s), sequences during running', ...
            i_sess, nSess, session.nom));
        t0 = tic;
        try
            statut = etape5_sequences(session, config);
            msg = sprintf('Step 5 - %d/%d [%s] (%s): %s', i_sess, nSess, session.nom, duree(toc(t0)), statut);
        catch err
            echecs(end+1, :) = {'step 5', session.nom, err.message}; %#ok<AGROW>
            msg = sprintf('Step 5 - %d/%d [%s] FAILED: %s', i_sess, nSess, session.nom, err.message);
        end
        historique = ecrire_log(log_file, historique, msg);
        appeler_rappel(config, session.nom, 5, msg);
    end
    signaler(config, 0.995, 'Step 5: sequences, all sessions');
    try
        agregation_sequences(sessions, config);
        historique = ecrire_log(log_file, historique, sprintf('Step 5 - sequences aggregated (%s)', ...
            fullfile(config.PathSave, 'sequences_agregation')));
    catch err
        echecs(end+1, :) = {'step 5', 'all', err.message};
        historique = ecrire_log(log_file, historique, ['Step 5 - aggregation FAILED: ' err.message]);
    end
end

%% Bilan
signaler(config, 1, sprintf('Finished (%d failure(s))', size(echecs, 1)));
historique = ecrire_log(log_file, historique, sprintf('=== END (total duration %s) ===', ...
    duree(toc(temps_global))));
if ~isempty(echecs)
    fprintf(2, '\n%d failure(s):\n', size(echecs, 1));
    disp(cell2table(echecs, 'VariableNames', {'step', 'session', 'error'}))
end
% Tous les paramètres de l'analyse, lus dans les fichiers (étapes 1 à 4)
try
    parametres_analyse(sessions, config, fullfile(config.PathSave, 'parametres_analyse.txt'));
catch err
    warning('parametres_analyse.txt non écrit : %s', err.message);
end
% Phrases récapitulatives de l'étape 4 : en dernier, dans la console, le
% journal et l'interface (fenêtre). Aussi dans agregation\resume_resultats.txt
if ~isempty(lignes_resume)
    fprintf('\n');
    for k = 1:numel(lignes_resume)
        historique = ecrire_log(log_file, historique, lignes_resume{k});
    end
    if isfield(config, 'rappel_resume') && ~isempty(config.rappel_resume)
        try
            config.rappel_resume(lignes_resume, fichier_resume);
        catch
        end
    end
end
end


function signaler(config, fraction, message)
% Avancement pour l'interface (config.rappel_progression(fraction, message)).
% Une erreur d'affichage n'arrête jamais l'analyse.
if ~isfield(config, 'rappel_progression') || isempty(config.rappel_progression), return; end
try
    config.rappel_progression(fraction, message);
catch
end
end


function appeler_rappel(config, nom, etape, msg)
% config.rappel_session(nom, etape, msg) : utilisé par interface_brainbow pour
% afficher la session qui vient d'être faite. Une erreur 'brainbow:arret'
% (arrêt demandé dans l'interface) arrête l'analyse ; les autres sont ignorées.
if ~isfield(config, 'rappel_session') || isempty(config.rappel_session), return; end
try
    config.rappel_session(nom, etape, msg);
catch err
    if strcmp(err.identifier, 'brainbow:arret'), rethrow(err); end
end
end


function t = duree(secondes)
t = char(string(duration(0, 0, secondes, 'Format', 'hh:mm:ss')));
end
