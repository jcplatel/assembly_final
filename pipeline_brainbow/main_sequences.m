function main_sequences(config)
%% =========================================================================
%  MAIN_SEQUENCES : y a-t-il des cellules activées en séquence pendant la
%  course ? Puis : les cellules de même couleur sont-elles plus proches dans
%  la séquence ? (étape 5)
%  -------------------------------------------------------------------------
%  Méthode de l'utilisateur (sequence_2026_01_27) : corrélation croisée de
%  chaque cellule avec les débuts de course, seuil d'Otsu, décalage < 0 ;
%  ordre appris sur la moitié des courses et vérifié sur l'autre moitié.
%  Couleurs : distance dans la séquence (ordre de toutes les courses), paires
%  de même couleur contre couleurs différentes, couleurs mélangées entre les
%  cellules colorées de la séquence ; Wilcoxon sur les animaux.
%  N'utilise ni les assemblées ni le choix du K.
%
%  Lancer :   main_sequences                   (toutes les sessions de la table)
%             config = config_brainbow();
%             config.seulement = {'444119'};   % une session ou un animal
%             main_sequences(config)
%  Réglages : config.seq dans config_brainbow.m (sinon options_sequences).
%
%  Par session : <session>\sequences\sequences.mat et sequence_session.png
%  Toutes sessions : <PathSave>\sequences_agregation\ (tableau, résumé,
%  sequences.png, sequence_colour.png)
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

fprintf('=== SÉQUENCES : %d session(s) ===\n', numel(a_traiter));
echecs = cell(0, 2);
for k = 1:numel(a_traiter)
    session = table2struct(sessions(a_traiter(k), :));
    t0 = tic;
    try
        statut = etape5_sequences(session, config);
        fprintf('%d/%d [%s] (%.0f s) : %s\n', k, numel(a_traiter), session.nom, toc(t0), statut);
    catch err
        echecs(end+1, :) = {session.nom, err.message}; %#ok<AGROW>
        fprintf(2, '%d/%d [%s] : ÉCHEC : %s\n', k, numel(a_traiter), session.nom, err.message);
    end
end

agregation_sequences(sessions, config);
if ~isempty(echecs)
    fprintf(2, '\n%d échec(s) :\n', size(echecs, 1));
    disp(cell2table(echecs, 'VariableNames', {'session', 'erreur'}))
end
end
