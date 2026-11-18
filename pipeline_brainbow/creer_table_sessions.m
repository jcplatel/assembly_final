function creer_table_sessions(config, id_condition)
%% Crée la table des sessions à partir de ton ancienne liste « matfile »
% (ex. nocues_data_recording_all.mat). À lancer UNE fois :
%     config = config_brainbow();
%     creer_table_sessions(config, 'sans_indice')
% Puis ouvrir le fichier Excel et vérifier, AVANT de lancer main_brainbow :
%   - id_animal, id_session, id_champ : les dates n'ont pas toutes le même
%     format (220919 / 11062025) ;
%   - id_champ : le plan (plane0, plane1...). C'est une étiquette : depuis
%     le 5 nov., chaque session est traitée comme son propre champ de vue
%     (un même plan revu un autre jour n'est en général pas à la même
%     profondeur) ;
%   - inclure = 0 pour une session à écarter ;
%   - fichier_colorcell : seulement pour forcer un fichier couleur ;
%   - um_par_px : si tu veux des blocs spatiaux en µm.
% Pour ajouter des sessions plus tard (ou une autre condition) :
% ajouter_sessions (ne touche pas aux lignes existantes).

if nargin < 2, id_condition = 'sans_indice'; end
if isfile(config.fichier_sessions)
    error('%s existe déjà : je ne l''écrase pas.', config.fichier_sessions);
end

L = load(config.fichier_liste_ancienne, 'matfile');   % liste dans une struct : pas de variable « matfile »
liste = L.matfile;
n = numel(liste);
T = lignes_sessions(liste, id_condition);
dossier = fileparts(config.fichier_sessions);
if ~isfolder(dossier), mkdir(dossier); end
writetable(T, config.fichier_sessions);
fprintf('Table écrite : %s (%d sessions). À vérifier avant de lancer main_brainbow.\n', ...
    config.fichier_sessions, n);
end
