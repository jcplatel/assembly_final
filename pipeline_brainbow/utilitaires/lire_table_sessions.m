function sessions = lire_table_sessions(config)
%LIRE_TABLE_SESSIONS  Lit la table des sessions (une ligne par session).
%
%   Colonnes obligatoires :
%     inclure          1 = traitée, 0 = écartée
%     id_animal        ex. '444119'
%     id_champ         ex. 'plane0' ; une étiquette (chaque session est son propre champ de vue)
%     id_session       ex. '220919' (date de l'enregistrement)
%     id_condition     ex. 'sans_indice' / 'avec_indice'
%     fichier_data     fichier .mat qui contient data_recording
%     dossier_suite2p  dossier qui contient Fall.mat, calcium.tif, colorcell_*.mat
%   Colonnes facultatives :
%     fichier_colorcell  pour forcer un fichier couleur (sinon GT > Maxwellnew > Maxwell)
%     um_par_px          taille du pixel (pour des blocs spatiaux en µm)
%
%   Ajoute : nom (= nom des dossiers de la session), dossier (étapes 2 à 4),
%   dossier_source (étape 1), graine.

% Pas encore de table : on la crée depuis l'ancienne liste « matfile »,
% puis on s'arrête pour qu'elle soit vérifiée avant toute analyse.
if ~isfile(config.fichier_sessions)
    if isempty(config.fichier_liste_ancienne) || ~isfile(config.fichier_liste_ancienne)
        error(['Table des sessions introuvable : %s\nEt ancienne liste introuvable : %s\n' ...
            'Renseigner config.fichier_liste_ancienne dans config_brainbow.m, ou créer la table à la main.'], ...
            config.fichier_sessions, config.fichier_liste_ancienne);
    end
    if isempty(config.conditions), condition = 'sans_indice'; else, condition = config.conditions{1}; end
    fprintf('Table des sessions absente : création depuis %s ...\n', config.fichier_liste_ancienne);
    creer_table_sessions(config, condition);
    error(['Table des sessions créée : %s\nOuvre-la et vérifie id_animal, id_session, id_champ ' ...
        '(le plan, ex. plane0) et inclure, puis relance main_brainbow.'], ...
        config.fichier_sessions);
end
opts_lecture = detectImportOptions(config.fichier_sessions);
texte = intersect({'id_animal', 'id_champ', 'id_session', 'id_condition', 'fichier_data', ...
    'dossier_suite2p', 'fichier_colorcell'}, opts_lecture.VariableNames, 'stable');
if ~isempty(texte)
    opts_lecture = setvartype(opts_lecture, texte, 'char');      % '444119' reste du texte
end
nombres = intersect({'inclure', 'um_par_px'}, opts_lecture.VariableNames, 'stable');
if ~isempty(nombres)
    opts_lecture = setvartype(opts_lecture, nombres, 'double');
end
sessions = readtable(config.fichier_sessions, opts_lecture);

obligatoires = {'inclure', 'id_animal', 'id_champ', 'id_session', 'id_condition', ...
    'fichier_data', 'dossier_suite2p'};
manque = setdiff(obligatoires, sessions.Properties.VariableNames);
if ~isempty(manque)
    error('Colonnes absentes de %s : %s', config.fichier_sessions, strjoin(manque, ', '));
end
if ~ismember('fichier_colorcell', sessions.Properties.VariableNames)
    sessions.fichier_colorcell = repmat({''}, height(sessions), 1);
end
if ~ismember('um_par_px', sessions.Properties.VariableNames)
    sessions.um_par_px = nan(height(sessions), 1);
end

% dossier_suite2p : le DOSSIER qui contient Fall.mat (ex. ...\suite2p\plane0).
% Si la table donne le chemin de Fall.mat lui-même, on garde son dossier.
for s = 1:height(sessions)
    d = strtrim(char(sessions.dossier_suite2p{s}));
    [parent, ~, ext] = fileparts(d);
    if strcmpi(ext, '.mat'), d = parent; end
    sessions.dossier_suite2p{s} = d;
end

% Sessions retenues
garde = sessions.inclure == 1;
if ~isempty(config.conditions)
    garde = garde & ismember(sessions.id_condition, config.conditions);
end
sessions = sessions(garde, :);
if isempty(sessions)
    error('Aucune session à traiter dans %s (inclure = 1, conditions : %s).', ...
        config.fichier_sessions, strjoin(config.conditions, ', '));
end

% Nom et dossier de résultats de chaque session (stables : pas de date d'analyse)
sessions.nom = strcat(sessions.id_animal, '_', sessions.id_session, '_', ...
    sessions.id_champ, '_', sessions.id_condition);
% dossier : étapes 2 à 4 (choix_K.mat, clonal, paires) ;
% dossier_source : étape 1 (prétraitement, k<K>). Les deux sont fixés par
% dossiers_analyse (le même dossier pour l'analyse principale).
sessions.dossier = fullfile(config.PathSave, sessions.nom);
source = config.PathSave;
if isfield(config, 'PathEtape1') && ~isempty(config.PathEtape1)
    source = config.PathEtape1;
elseif isfield(config, 'PathSource') && ~isempty(config.PathSource)
    source = config.PathSource;
end
sessions.dossier_source = fullfile(source, sessions.nom);

% Une session ne doit apparaître qu'une fois
[noms_u, ~, ic] = unique(sessions.nom);
n_fois = accumarray(ic, 1);
if any(n_fois > 1)
    error('Sessions en double dans la table : %s', strjoin(noms_u(n_fois > 1), ', '));
end

% Graine propre à chaque session, indépendante de l'ordre des lignes
sessions.graine = cellfun(@(s) config.graine + mod(sum(double(s) .* (1:numel(s))), 1e6), sessions.nom);
end
