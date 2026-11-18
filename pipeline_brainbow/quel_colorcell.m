function T = quel_colorcell(config, toutes)
%% =========================================================================
%  QUEL_COLORCELL : quel fichier de couleurs le pipeline utilise pour chaque
%  session (13 nov., demandé par l'utilisateur)
%  -------------------------------------------------------------------------
%  Même règle que charger_couleurs (la seule fonction du pipeline qui lit
%  les couleurs) : le fichier de la colonne fichier_colorcell de la table
%  s'il est rempli, sinon colorcell_GT, puis colorcell_Maxwellnew, puis
%  colorcell_Maxwell dans le dossier suite2p.
%  Ne modifie rien.
%
%  Lancer :   quel_colorcell                 sessions de la table avec inclure = 1
%             quel_colorcell(config)          idem, avec ce config
%             quel_colorcell(config, true)    TOUTES les lignes de la table
%                                             (inclure 0 ou 1, toutes conditions)
%
%  Colonnes du tableau :
%    source          table / GT / Maxwellnew / Maxwell / AUCUN / INTROUVABLE
%    fichier_utilise chemin complet du fichier lu
%    date_fichier    date de modification de ce fichier
%    n_valeurs       nombre de valeurs dans colorcell
%    n_iscell        nombre de cellules iscell = 1 dans Fall.mat
%    taille_ok       1 si n_valeurs = n_iscell
%    autres_fichiers les autres fichiers colorcell*.mat du dossier suite2p
%                    (pas utilisés)
%  Écrit aussi <dossier de l'étape 1>\verification_couleurs\colorcell_utilise.csv
%  =========================================================================

if nargin < 1 || isempty(config), config = config_brainbow(); end
if nargin < 2 || isempty(toutes), toutes = false; end
addpath(genpath(fileparts(mfilename('fullpath'))));
config = dossiers_analyse(config);

if toutes
    sessions = lire_table_complete(config.fichier_sessions);
else
    sessions = lire_table_sessions(config);
end
if isfield(config, 'seulement') && ~isempty(config.seulement)
    seulement = cellstr(config.seulement);
    sessions = sessions(ismember(sessions.nom, seulement) | ismember(sessions.id_animal, seulement), :);
end

n = height(sessions);
source = repmat({''}, n, 1); fichier_utilise = repmat({''}, n, 1); date_fichier = repmat({''}, n, 1);
n_valeurs = nan(n, 1); n_iscell = nan(n, 1); taille_ok = nan(n, 1);
autres_fichiers = repmat({''}, n, 1);
for s = 1:n
    path = avec_separateur(char(sessions.dossier_suite2p{s}));
    force = strtrim(char(sessions.fichier_colorcell{s}));
    % même ordre que charger_couleurs
    if ~isempty(force)
        f = force; source{s} = 'table';
    elseif isfile([path 'colorcell_GT.mat'])
        f = [path 'colorcell_GT.mat']; source{s} = 'GT';
    elseif isfile([path 'colorcell_Maxwellnew.mat'])
        f = [path 'colorcell_Maxwellnew.mat']; source{s} = 'Maxwellnew';
    elseif isfile([path 'colorcell_Maxwell.mat'])
        f = [path 'colorcell_Maxwell.mat']; source{s} = 'Maxwell';
    else
        f = ''; source{s} = 'AUCUN';
    end
    fichier_utilise{s} = f;
    if ~isempty(f)
        if isfile(f)
            d = dir(f);
            date_fichier{s} = datestr(d.datenum, 'yyyy-mm-dd HH:MM');
            try
                S = load(f, 'colorcell');
                n_valeurs(s) = numel(S.colorcell);
            catch
                source{s} = [source{s} ' (pas de variable colorcell)'];
            end
        else
            source{s} = [source{s} ' INTROUVABLE'];
        end
    end
    if isfile([path 'Fall.mat'])
        try
            F = load([path 'Fall.mat'], 'iscell');
            n_iscell(s) = sum(F.iscell(:, 1) == 1);
        catch
        end
    end
    if ~isnan(n_valeurs(s)) && ~isnan(n_iscell(s))
        taille_ok(s) = double(n_valeurs(s) == n_iscell(s));
    end
    % autres fichiers de couleurs présents dans le dossier (non utilisés)
    L = dir([path 'colorcell*.mat']);
    autres = {L.name};
    [~, nom_utilise, ext] = fileparts(f);
    autres = autres(~strcmpi(autres, [nom_utilise ext]));
    autres_fichiers{s} = strjoin(autres, ', ');
end

T = table(sessions.nom, sessions.inclure, sessions.id_condition, source, fichier_utilise, date_fichier, ...
    n_valeurs, n_iscell, taille_ok, autres_fichiers, ...
    'VariableNames', {'session', 'inclure', 'condition', 'source', 'fichier_utilise', 'date_fichier', ...
    'n_valeurs', 'n_iscell', 'taille_ok', 'autres_fichiers'});

% affichage
fprintf('\n=== FICHIER COULEUR UTILISÉ : %d session(s) ===\n', n);
for s = 1:n
    if taille_ok(s) == 1, t = 'taille ok';
    elseif taille_ok(s) == 0, t = sprintf('TAILLE DIFFÉRENTE (%d valeurs, %d iscell)', n_valeurs(s), n_iscell(s));
    else, t = 'taille non vérifiée';
    end
    fprintf('%-45s inclure %d  %-11s %s  [%s]\n', T.session{s}, T.inclure(s), source{s}, fichier_utilise{s}, t);
end
sources = unique(source);
fprintf('\nRésumé :');
for i = 1:numel(sources)
    fprintf('  %s : %d', sources{i}, sum(strcmp(source, sources{i})));
end
fprintf('\n');
nb = sum(taille_ok == 0);
if nb > 0, fprintf('%d session(s) avec une taille différente de iscell.\n', nb); end

dossier = fullfile(config.PathEtape1, 'verification_couleurs');
try
    if ~isfolder(dossier), mkdir(dossier); end
    writetable(T, fullfile(dossier, 'colorcell_utilise.csv'));
    fprintf('Tableau : %s\n', fullfile(dossier, 'colorcell_utilise.csv'));
catch err
    fprintf(2, 'Tableau non écrit (%s)\n', err.message);
end
end

function sessions = lire_table_complete(fichier)
% Toute la table, sans filtre inclure ni condition (même lecture que
% lire_table_sessions).
opts_lecture = detectImportOptions(fichier);
texte = intersect({'id_animal', 'id_champ', 'id_session', 'id_condition', 'fichier_data', ...
    'dossier_suite2p', 'fichier_colorcell'}, opts_lecture.VariableNames, 'stable');
if ~isempty(texte), opts_lecture = setvartype(opts_lecture, texte, 'char'); end
if ismember('inclure', opts_lecture.VariableNames)
    opts_lecture = setvartype(opts_lecture, 'inclure', 'double');
end
sessions = readtable(fichier, opts_lecture);
if ~ismember('fichier_colorcell', sessions.Properties.VariableNames)
    sessions.fichier_colorcell = repmat({''}, height(sessions), 1);
end
sessions.fichier_colorcell(cellfun(@(x) ~ischar(x), sessions.fichier_colorcell)) = {''};
for s = 1:height(sessions)
    d = strtrim(char(sessions.dossier_suite2p{s}));
    [parent, ~, ext] = fileparts(d);
    if strcmpi(ext, '.mat'), d = parent; end
    sessions.dossier_suite2p{s} = d;
end
sessions = sessions(~cellfun(@isempty, sessions.dossier_suite2p), :);
sessions.nom = strcat(sessions.id_animal, '_', sessions.id_session, '_', ...
    sessions.id_champ, '_', sessions.id_condition);
end
