function ajouter_sessions(config, fichiers, id_condition)
%AJOUTER_SESSIONS  Ajoute de nouvelles sessions à la table des sessions,
%   sans toucher aux lignes existantes.
%
%   config = config_brainbow();
%   ajouter_sessions(config)
%       ajoute les fichiers de ta liste « matfile » (config.fichier_liste_ancienne)
%       qui ne sont pas encore dans la table
%   ajouter_sessions(config, {'E:\...\data_1.mat', 'E:\...\data_2.mat'})
%       ajoute ces fichiers (fichiers .mat qui contiennent data_recording)
%   ajouter_sessions(config, fichiers, 'avec_indice')
%       même chose pour une autre condition
%
%   Ensuite : vérifier les nouvelles lignes dans l'Excel (id_animal,
%   id_session, id_champ), puis relancer main_brainbow. Seules les nouvelles
%   sessions passent par l'étape 1 : les autres sont gardées en cache tant
%   que les paramètres n'ont pas changé.

if nargin < 3 || isempty(id_condition)
    if isempty(config.conditions), id_condition = 'sans_indice'; else, id_condition = config.conditions{1}; end
end
if ~isfile(config.fichier_sessions)
    error('Pas encore de table %s : la créer d''abord avec creer_table_sessions(config).', ...
        config.fichier_sessions);
end
if nargin < 2 || isempty(fichiers)
    L = load(config.fichier_liste_ancienne, 'matfile');
    fichiers = L.matfile;
end
if ischar(fichiers) || isstring(fichiers), fichiers = cellstr(fichiers); end
fichiers = cellfun(@(f) char(string(f)), fichiers(:)', 'UniformOutput', false);

% Table actuelle, toutes les lignes (y compris inclure = 0)
opts_lecture = detectImportOptions(config.fichier_sessions);
texte = intersect({'id_animal', 'id_champ', 'id_session', 'id_condition', 'fichier_data', ...
    'dossier_suite2p', 'fichier_colorcell'}, opts_lecture.VariableNames, 'stable');
opts_lecture = setvartype(opts_lecture, texte, 'char');
T = readtable(config.fichier_sessions, opts_lecture);

deja = ismember(lower(fichiers), lower(T.fichier_data));
nouveaux = fichiers(~deja);
if isempty(nouveaux)
    fprintf('Aucune nouvelle session : ces %d fichiers sont déjà dans la table.\n', numel(fichiers));
    return
end
fprintf('%d nouvelle(s) session(s) :\n', numel(nouveaux));
N = lignes_sessions(nouveaux, id_condition);

% Mêmes colonnes que la table (y compris des colonnes que tu aurais ajoutées)
for v = setdiff(T.Properties.VariableNames, N.Properties.VariableNames)
    if isnumeric(T.(v{1}))
        N.(v{1}) = nan(height(N), 1);
    else
        N.(v{1}) = repmat({''}, height(N), 1);
    end
end
T = [T; N(:, T.Properties.VariableNames)];

% Écrit d'abord une copie, puis remplace : la table n'est jamais perdue
[dossier, nom, ext] = fileparts(config.fichier_sessions);
copie = fullfile(dossier, [nom '_tmp' ext]);
writetable(T, copie);
[ok, message] = movefile(copie, config.fichier_sessions, 'f');
if ~ok
    error('Impossible de remplacer %s (ouvert dans Excel ?) : %s\nLa nouvelle table est dans %s.', ...
        config.fichier_sessions, message, copie);
end
fprintf('Table mise à jour : %s (%d sessions dont %d nouvelles). Vérifier les nouvelles lignes.\n', ...
    config.fichier_sessions, height(T), height(N));
end
