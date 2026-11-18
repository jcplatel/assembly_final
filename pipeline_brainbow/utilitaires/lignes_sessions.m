function T = lignes_sessions(liste, id_condition)
%LIGNES_SESSIONS  Lignes de la table des sessions pour une liste de fichiers
%   qui contiennent data_recording (utilisé par creer_table_sessions et
%   ajouter_sessions). Les identifiants sont lus dans data_recording.identifier,
%   sinon dans le nom du fichier (forme 444119_220919_plane0).
if ischar(liste) || isstring(liste), liste = cellstr(liste); end
n = numel(liste);
inclure           = ones(n, 1);
id_animal         = cell(n, 1);
id_champ          = cell(n, 1);
id_session        = cell(n, 1);
id_cond           = repmat({id_condition}, n, 1);
fichier_data      = cell(n, 1);
dossier_suite2p   = cell(n, 1);
fichier_colorcell = repmat({''}, n, 1);
um_par_px         = nan(n, 1);

for i = 1:n
    filename = char(string(liste{i}));
    fichier_data{i} = filename;
    S = load(filename, 'data_recording');
    identifier = char(string(S.data_recording.identifier));
    tok = regexp(identifier, '(\d+)_(\d+)_(plane\d+)', 'tokens', 'once');
    if isempty(tok)
        [~, nom] = fileparts(filename);
        tok = regexp(nom, '(\d+)_(\d+)_(plane\d+)', 'tokens', 'once');
    end
    if isempty(tok)
        warning('Identifiants non reconnus pour %s : à remplir à la main (inclure = 0 en attendant).', filename);
        tok = {'?', '?', '?'};
        inclure(i) = 0;
    end
    id_animal{i}  = tok{1};
    id_session{i} = tok{2};
    id_champ{i}   = tok{3};
    dossier_suite2p{i} = char(string(S.data_recording.path));
    fprintf('%d/%d  %s\n', i, n, identifier);
end

T = table(inclure, id_animal, id_champ, id_session, id_cond, fichier_data, dossier_suite2p, ...
    fichier_colorcell, um_par_px, 'VariableNames', {'inclure', 'id_animal', 'id_champ', ...
    'id_session', 'id_condition', 'fichier_data', 'dossier_suite2p', 'fichier_colorcell', 'um_par_px'});
end
