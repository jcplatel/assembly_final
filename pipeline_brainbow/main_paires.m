function main_paires(config)
%% =========================================================================
%  MAIN_PAIRES : les cellules voisines de même couleur sont-elles plus
%  synchrones au repos que les voisines de couleurs différentes ?
%  -------------------------------------------------------------------------
%  Indépendant des assemblées : n'utilise que le prétraitement de l'étape 1
%  (pretraitement.mat), les contours suite2p et les couleurs.
%
%  Lancer :   main_paires                      (toutes les sessions de la table)
%             config = config_brainbow();
%             config.seulement = {'444119_220919_plane0_sans_indice'};   % ou {'444119'}
%             main_paires(config)
%  Réglages : config.paires dans config_brainbow.m (à figer avant de regarder).
%
%  Par session : <session>\paires\paires.mat et paires.png
%  Toutes sessions : <PathSave>\paires_agregation\ (tableaux et figure)
%  Avec config.analyse = 'nom' : dans <PathSave>\analyses\nom\ (voir dossiers_analyse)
%  =========================================================================

if nargin < 1, config = config_brainbow(); end
dossier_pipeline = fileparts(mfilename('fullpath'));
addpath(genpath(dossier_pipeline));
etat_rng = rng;
remettre_rng = onCleanup(@() rng(etat_rng)); %#ok<NASGU>

config = dossiers_analyse(config);      % étape 1 dans PathSave ; résultats dans l'analyse choisie
sessions = lire_table_sessions(config);
a_traiter = 1:height(sessions);
if isfield(config, 'seulement') && ~isempty(config.seulement)
    seulement = cellstr(config.seulement);
    a_traiter = find(ismember(sessions.nom, seulement) | ismember(sessions.id_animal, seulement))';
    if isempty(a_traiter)
        error('config.seulement : aucune session ne correspond à %s.', strjoin(seulement, ', '));
    end
end

fprintf('=== PAIRES : %d session(s) ; mesures %s ; permutations %s ===\n', numel(a_traiter), ...
    strjoin(config.paires.mesures, ', '), strjoin(config.paires.permutations, ', '));
echecs = cell(0, 2);
for k = 1:numel(a_traiter)
    session = table2struct(sessions(a_traiter(k), :));
    t0 = tic;
    if ~isfile(fullfile(session.dossier_source, 'pretraitement.mat'))
        fprintf('%d/%d [%s] : pas de prétraitement (étape 1 à faire), session sautée\n', ...
            k, numel(a_traiter), session.nom);
        continue
    end
    try
        etape3b_paires(session, config);
        fprintf('%d/%d [%s] : ok (%.0f s)\n', k, numel(a_traiter), session.nom, toc(t0));
    catch err
        echecs(end+1, :) = {session.nom, err.message}; %#ok<AGROW>
        fprintf(2, '%d/%d [%s] : ÉCHEC : %s\n', k, numel(a_traiter), session.nom, err.message);
    end
end

agregation_paires(sessions, config);
if ~isempty(echecs)
    fprintf(2, '\n%d échec(s) :\n', size(echecs, 1));
    disp(cell2table(echecs, 'VariableNames', {'session', 'erreur'}))
end
end
