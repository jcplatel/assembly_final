function main_hierarchique(config)
%% =========================================================================
%  MAIN_HIERARCHIQUE : clustering hiérarchique des cellules, puis couleur
%  de chaque groupe de l'arbre (étape 3c, exploratoire, 7 nov.)
%  -------------------------------------------------------------------------
%  Des petits groupes de cellules qui s'activent ensemble aux gros groupes :
%  jusqu'à quelle taille trouve-t-on des groupes avec une couleur en excès ?
%  Trois arbres, sur les mêmes cellules (au moins config.hier.min_sce SCE) :
%    'sce'  : co-participation aux SCE (corrélation des lignes de Race) ;
%    'sttc' : STTC au repos, distance et niveau de chaque cellule retirés
%             (comme l'étape 3b) ;
%    'fca'  : Functional Clustering Algorithm (Feldt et al. 2009, utilisé
%             par Feldt Muldoon et al. 2013), plus ses clusters seuls.
%  Deux hasards, sur le MÊME arbre : couleurs mélangées dans le champ
%  (principal) ; par blocs de config.hier.taille_bloc_um (contrôle spatial).
%  N'utilise ni les assemblées ni le choix du K.
%
%  Lancer :   main_hierarchique                (toutes les sessions de la table)
%             config = config_brainbow();
%             config.seulement = {'444119'};   % une session ou un animal
%             main_hierarchique(config)
%  Réglages : config.hier dans config_brainbow.m (sinon options_hierarchique).
%
%  Par session : <session>\hierarchique\hierarchique.mat et arbres.png
%  Toutes sessions : <PathSave>\hierarchique_agregation\ (tableaux, résumé,
%  hierarchique.png)
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
opts = options_hierarchique(config);
fprintf('=== HIÉRARCHIQUE : %d session(s) ; arbres %s ; hasards %s ; %d mélanges ===\n', numel(a_traiter), ...
    strjoin(opts.arbres, ', '), strjoin(opts.permutations, ', '), opts.nPerm);
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
        statut = etape3c_hierarchique(session, config);
        fprintf('%d/%d [%s] (%.0f s) : %s\n', k, numel(a_traiter), session.nom, toc(t0), statut);
    catch err
        echecs(end+1, :) = {session.nom, err.message}; %#ok<AGROW>
        fprintf(2, '%d/%d [%s] : ÉCHEC : %s\n', k, numel(a_traiter), session.nom, err.message);
    end
end

agregation_hierarchique(sessions, config);
if ~isempty(echecs)
    fprintf(2, '\n%d échec(s) :\n', size(echecs, 1));
    disp(cell2table(echecs, 'VariableNames', {'session', 'erreur'}))
end
end
