function main_soeurs(config)
%% =========================================================================
%  MAIN_SOEURS : la similarité d'activité d'une paire augmente-t-elle avec
%  sa probabilité d'être sœurs ? (étape 3d, 7 nov.)
%  -------------------------------------------------------------------------
%  Pourquoi : la couleur n'est pas le clone. Parmi les paires de même
%  couleur, seule une petite part sont de vraies sœurs (beaucoup de clones
%  indépendants partagent une couleur), et cette part baisse avec la
%  distance. Comparer « même couleur contre couleur différente » dilue donc
%  fortement tout effet clonal.
%  1. Pour chaque session (etape_soeurs_session) : comptes de paires de même
%     couleur par distance (couleurs et positions seulement) ; similarité
%     d'activité de chaque paire, distance et niveau de chaque cellule
%     retirés (activité seulement) : STTC au repos, même assemblée.
%  2. Toutes les sessions (carte_soeurs) : probabilité d'être sœurs pi selon
%     la couleur et la distance = part de l'excès de paires de même couleur ;
%     pi = 0 au-delà de 200 µm (taille maximale d'un clone). Direction des
%     bandes de couleur de chaque session, à vérifier à l'œil
%     (figure_angles_bandes).
%  3. Pour chaque session (test_soeurs_session) : pente de la similarité
%     sur pi, à couleur égale (seules les paires de même couleur proches et
%     lointaines la déterminent) ; couleurs mélangées pour le p de session.
%  4. Toutes les sessions (agregation_soeurs) : un test par mesure,
%     Wilcoxon sur les animaux ; effet chez les vraies sœurs en % ; figure.
%  5. Chaque assemblée (coeurs_assemblees, agregation_coeurs) : score de
%     cœur clonal (paires de sœurs attendues parmi ses membres), p contre
%     les couleurs mélangées et contre des voisines, FDR, verdict.
%  La co-appartenance aux assemblées demande les étapes 2 et 3 (K retenu).
%
%  Lancer :   main_soeurs
%             config = config_brainbow(); config.seulement = {'444119'}; main_soeurs(config)
%  ATTENTION : avec config.seulement, pi est estimé sur ces sessions-là
%  seulement (moins précis) ; pour le résultat final, toutes les sessions.
%  Réglages : config.soeurs dans config_brainbow.m (sinon options_soeurs).
%  Sorties : <session>\soeurs\soeurs.mat, coeurs.mat ; <PathSave>\soeurs_agregation\
%  (carte_soeurs.mat, angles_bandes.csv, angles_bandes_<animal>.png,
%  soeurs_sessions.csv, resume_soeurs.txt, soeurs.png,
%  coeurs_assemblees.csv, resume_coeurs.txt, coeurs_assemblees.png).
%  =========================================================================

if nargin < 1, config = config_brainbow(); end
dossier_pipeline = fileparts(mfilename('fullpath'));
addpath(genpath(dossier_pipeline));
etat_rng = rng;
remettre_rng = onCleanup(@() rng(etat_rng)); %#ok<NASGU>

config = dossiers_analyse(config);
sessions = lire_table_sessions(config);
a_traiter = 1:height(sessions);
if isfield(config, 'seulement') && ~isempty(config.seulement)
    seulement = cellstr(config.seulement);
    a_traiter = find(ismember(sessions.nom, seulement) | ismember(sessions.id_animal, seulement))';
    if isempty(a_traiter)
        error('config.seulement : aucune session ne correspond à %s.', strjoin(seulement, ', '));
    end
end
opts = options_soeurs(config);
fprintf('=== SŒURS : %d session(s) ; mesures %s ===\n', numel(a_traiter), strjoin(opts.mesures, ', '));

% 1. comptes (couleurs) et résidus (activité), session par session
echecs = cell(0, 2);
ok = false(1, height(sessions));
for k = 1:numel(a_traiter)
    session = table2struct(sessions(a_traiter(k), :));
    t0 = tic;
    if ~isfile(fullfile(session.dossier_source, 'pretraitement.mat'))
        fprintf('%d/%d [%s] : pas de prétraitement (étape 1 à faire), session sautée\n', ...
            k, numel(a_traiter), session.nom);
        continue
    end
    try
        statut = etape_soeurs_session(session, config);
        ok(a_traiter(k)) = true;
        fprintf('%d/%d [%s] (%.0f s) : %s\n', k, numel(a_traiter), session.nom, toc(t0), statut);
    catch err
        echecs(end+1, :) = {session.nom, err.message}; %#ok<AGROW>
        fprintf(2, '%d/%d [%s] : ÉCHEC : %s\n', k, numel(a_traiter), session.nom, err.message);
    end
end

% 2. probabilité d'être sœurs, sur les sessions traitées
K = carte_soeurs(sessions(ok, :), config);
fprintf('Probabilité d''être sœurs estimée sur %d session(s).\n', K.n_sessions);
try
    figure_angles_bandes(sessions(ok, :), fullfile(config.PathSave, 'soeurs_agregation'));
catch err
    warning('figure_angles_bandes non exécutée : %s', err.message);
end

% 3. test, session par session
for k = find(ok)
    session = table2struct(sessions(k, :));
    try
        test_soeurs_session(session, config, K);
    catch err
        echecs(end+1, :) = {session.nom, ['test : ' err.message]}; %#ok<AGROW>
        fprintf(2, '[%s] : ÉCHEC du test : %s\n', session.nom, err.message);
    end
end

% 4. toutes les sessions
agregation_soeurs(sessions(ok, :), config, K);

% 5. cœur clonal de chaque assemblée (étapes 2-3 nécessaires)
for k = find(ok)
    session = table2struct(sessions(k, :));
    try
        coeurs_assemblees(session, config, K);
    catch err
        echecs(end+1, :) = {session.nom, ['cœurs : ' err.message]}; %#ok<AGROW>
        fprintf(2, '[%s] : ÉCHEC des cœurs clonaux : %s\n', session.nom, err.message);
    end
end
agregation_coeurs(sessions(ok, :), config);
if ~isempty(echecs)
    fprintf(2, '\n%d échec(s) :\n', size(echecs, 1));
    disp(cell2table(echecs, 'VariableNames', {'session', 'erreur'}))
end
end
