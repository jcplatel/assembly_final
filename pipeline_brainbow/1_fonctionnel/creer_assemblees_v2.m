function creer_assemblees_v2(f1, f2)
%CREER_ASSEMBLEES_V2  Écrit la version 2 des assemblées d'un K (11 nov.)
%   f1 : <étape 1>\<session>\k<K>\assemblees.mat (version 1, RACEORTHO_FDR2)
%   f2 : <étape 1>\<session>\k<K>\assemblees_v2.mat
%   Reprend tout le fichier de la version 1 (clusters de SCE, silhouettes,
%   p des cellules...) et remplace seulement les assemblées par celles de
%   RACEORTHO_V2 (cellules partagées entre assemblées). Le k-means n'est pas
%   refait. Champs ajoutés : assemblee_principale, n_appartenances,
%   version_assemblees = 'v2', source_v1, n_assemblees_v1.
%   Restent ceux de la version 1 (non recalculés) : resume_K,
%   Cellular_Fidelity, Event_Purity, AmpTable et leurs moyennes.
A = load(f1);
dossier_session = fileparts(fileparts(f1));
[~, nom_K] = fileparts(fileparts(f1));
[~, nom_session] = fileparts(dossier_session);
K = sscanf(nom_K, 'k%d');
if isempty(K), K = 0; end
graine = mod(sum(double(nom_session) .* (1:numel(nom_session))), 1e6) + 1000 * K;
A.n_assemblees_v1 = 0;
if isfield(A, 'assemblyortho'), A.n_assemblees_v1 = numel(A.assemblyortho); end
A.source_v1 = f1;
A.version_assemblees = 'v2';
if ~all(isfield(A, {'IDX2', 'NClOK', 'pValues_cells', 'config_asm'})) || isempty(A.pValues_cells)
    A.statut = 'v2 impossible : clusters ou p des cellules absents (v1 gardée)';
    save(f2, '-struct', 'A');
    return
end
P = load(fullfile(dossier_session, 'pretraitement.mat'), 'Race');
R = RACEORTHO_V2(double(P.Race), A.IDX2, A.NClOK, A.pValues_cells, A.config_asm, graine);
for ch = {'assemblyortho', 'cluster_origine', 'PCl_final', 'CellCl_final', 'pValues_race', 'x1', 'RList', ...
        'assemblee_principale', 'n_appartenances'}
    A.(ch{1}) = R.(ch{1});
end
A.statut = R.statut;
if isfield(A, 'NCl'), A.NCl = numel(R.assemblyortho); end
save(f2, '-struct', 'A');
end
