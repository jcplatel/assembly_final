function statut = etape3_clonal(session, config, id_execution)
%% =========================================================================
%  ÉTAPE 3 : statistiques clonales sur le K retenu
%  SEULE étape du pipeline qui charge les couleurs.
%  -------------------------------------------------------------------------
%  Remplace new_test_brainbow (pour une session), la partie couleur de
%  brainbowassemblies2026_05_07, distance_calculation et
%  export_data_brainbow. Écrit dans <session>\clonal\ :
%     stats_clonales.mat   res (variante principale), res_variantes, statut,
%                          dist (distances en µm), resume_session, ...
%     figures              cartes et camemberts
%  L'ancien stats_clonales.mat est supprimé AVANT le calcul : une session en
%  échec ne laisse jamais un ancien résultat derrière elle.
%  Une session exclue pour une raison prévue (pas d'assemblée, ...) écrit
%  un stats_clonales.mat avec statut = 'exclue : raison'.
%  =========================================================================

dossier_clonal = fullfile(session.dossier, 'clonal');
if ~isfolder(dossier_clonal), mkdir(dossier_clonal); end
fichier_clonal = fullfile(dossier_clonal, 'stats_clonales.mat');
if isfile(fichier_clonal), delete(fichier_clonal); end
namefull = avec_separateur(dossier_clonal);

% Identifiants : viennent de la table, plus du nom des dossiers
id_animal    = session.id_animal;
id_champ     = session.id_champ;
id_session   = session.id_session;
id_condition = session.id_condition;

%% 1. K retenu à l'étape 2 (0 = consensus entre les K)
choix = load(fullfile(session.dossier, 'choix_K.mat'));
if startsWith(choix.statut_choix, 'erreur')
    error('Étape 2 en erreur pour cette session : %s', choix.statut_choix);
elseif ~strcmp(choix.statut_choix, 'ok')
    statut = choix.statut_choix;
    enregistrer_exclusion(fichier_clonal, statut, id_execution);
    return
end
id_K = choix.best_K;

%% 2. Assemblées du K retenu (aucune couleur dans ces fichiers)
% K = 0 : assemblées de consensus entre les K (étape 2, règle 'consensus'),
% dans le dossier de l'analyse ; sinon k<K> de l'étape 1
fichier_asm = fichier_assemblees(session.dossier, session.dossier_source, id_K);   % v1 ou v2 selon l'analyse
A = load(fichier_asm, 'assemblyortho', 'PCl_final', 'resume_K');
assembly     = A.assemblyortho;
nbassemblies = numel(assembly);
NCl          = nbassemblies;
PCl_final    = A.PCl_final;
fichier_pre = fullfile(session.dossier_source, 'pretraitement.mat');
pre = load(fichier_pre, 'Race', 'Raster', 'TRace', 'th_detection', 'Fzero');
Race = double(pre.Race);
th_detection = pre.th_detection;
F0_mediane = median(double(pre.Fzero), 2, 'omitnan');
pre = rmfield(pre, 'Fzero');
% diagnostic du prétraitement (cellules hyperactives) : absent des très vieux fichiers
diag_pre = struct();
if ismember('diagnostic', who('-file', fichier_pre))
    Sd = load(fichier_pre, 'diagnostic');
    diag_pre = Sd.diagnostic;
end
if nbassemblies == 0
    statut = 'exclue : aucune assemblée au K retenu';
    enregistrer_exclusion(fichier_clonal, statut, id_execution);
    return
end
% Score de séquentialité : pas calculé par le pipeline. S'il existe un jour
% dans assemblees.mat (variable score_sequence, une valeur par assemblée,
% définie SANS les couleurs), il sera utilisé.
score_sequence = nan(nbassemblies, 1);
contenu = whos('-file', fichier_asm);
if ismember('score_sequence', {contenu.name})
    S = load(fichier_asm, 'score_sequence');
    score_sequence = S.score_sequence(:);
end

%% 3. Géométrie et couleurs
path = avec_separateur(session.dossier_suite2p);
[outline_gcampx, outline_gcampy, gcamp_centroid, calcium] = charger_geometrie(path);
[colorcell, path_colorcell, source_couleur] = charger_couleurs(path, session.fichier_colorcell);
colorcell = couleurs_sans_doublons(colorcell, session);    % doublons de l'étape 1 : non colorés
verifier_alignement(size(Race, 1), colorcell, outline_gcampx, session.nom);

%% 4. Biais de détection : bruit et F0 selon la couleur
[p_detection, p_F0, p_composition_membres, compte_membres, compte_non_membres] = ...
    test_biais_detection(colorcell, assembly, th_detection, F0_mediane, config.clonal.nPerm, session.graine);

%% 5. Statistiques clonales, pour chaque variante (la 1re est le test principal)
taille_bloc_px = config.clonal.taille_bloc_px;
if ~isnan(config.clonal.taille_bloc_um) && ~isnan(session.um_par_px)
    taille_bloc_px = config.clonal.taille_bloc_um / session.um_par_px;
end
% 'territoire' : marge autour de l'enveloppe des cellules de l'assemblée (µm -> pixels)
um_px = session.um_par_px;
if isnan(um_px), um_px = config.um_par_px_defaut; end
marge_um = 10;
if isfield(config.clonal, 'marge_territoire_um'), marge_um = config.clonal.marge_territoire_um; end
D = struct();
D.marge_territoire_px = marge_um / um_px;
D.colorcell      = colorcell;
D.assembly       = assembly;
D.gcamp_centroid = gcamp_centroid;
D.Race           = Race;
D.PCl_final      = PCl_final;
D.score_sequence = score_sequence;
D.id_animal      = id_animal;
D.id_champ       = id_champ;
D.id_session     = id_session;
D.id_condition   = id_condition;
D.id_K           = id_K;
D.source_couleur = source_couleur;
D.taille_bloc_px = taille_bloc_px;
D.activite       = sum(Race > 0, 2);    % nombre de SCE de chaque cellule (permutation 'activite')
D.p_detection    = p_detection;

variantes = config.clonal.variantes;
nV = size(variantes, 1);
res_variantes  = cell(nV, 1);
noms_variantes = cell(nV, 1);
Kcol_variantes = cell(nV, 1);
for v = 1:nV
    [res_variantes{v}, Kcol_variantes{v}] = stats_clonales_champ3(D, config.clonal, ...
        variantes{v, 1}, variantes{v, 2}, session.graine + 10 * v);
    noms_variantes{v} = sprintf('%s_%s', variantes{v, 1}, variantes{v, 2});
end
res = res_variantes{1};

%% 6. Distances (ancien distance_calculation) et ligne de résumé (ancien export_data_brainbow)
size_pixel = session.um_par_px;                 % µm par pixel : table, sinon valeur par défaut
if isnan(size_pixel), size_pixel = config.um_par_px_defaut; end
dist = distance_calculation2(gcamp_centroid, assembly, colorcell, size_pixel, ...
    config.clonal.nshuffle_distance, session.graine + 100);
resume_session = export_data_brainbow2(A.resume_K, res, dist, p_detection, p_F0, ...
    source_couleur, config.clonal.seuilstat);

%% 6b. Une ligne par cellule, pour les contrôles du champ de l'étape 4
% (couleurs regroupées dans l'espace ; détection selon la couleur)
NCell = size(Race, 1);
membre = false(NCell, 1);
for i_asm = 1:nbassemblies
    cel_a = assembly{i_asm}(:);
    membre(cel_a(cel_a >= 1 & cel_a <= NCell)) = true;
end
hyperactive = false(NCell, 1);
if isfield(diag_pre, 'cellules_hyperactives_idx')
    idx_h = diag_pre.cellules_hyperactives_idx(:);
    hyperactive(idx_h(idx_h >= 1 & idx_h <= NCell)) = true;
end
cellules = struct();
cellules.couleur      = colorcell(:);
cellules.xy_um        = double(gcamp_centroid) * um_px;
cellules.bruit        = double(th_detection(:));     % seuil de détection = 3,09 x bruit (dF/F0)
cellules.F0           = double(F0_mediane(:));
cellules.n_transients = full(sum(pre.Raster > 0, 2)); % transients au repos
cellules.n_SCE        = sum(Race > 0, 2);             % participation aux SCE
cellules.NRace        = size(Race, 2);
cellules.membre       = membre;                       % dans une assemblée du K retenu
cellules.hyperactive  = hyperactive;

statut = 'ok';
opts_clonal = config.clonal;                 % réglages utilisés (bouton « Analysis parameters »)
opts_clonal.um_par_px_defaut = config.um_par_px_defaut;
save(fichier_clonal, 'res', 'res_variantes', 'noms_variantes', 'Kcol_variantes', 'statut', ...
    'id_execution', 'p_detection', 'p_F0', 'p_composition_membres', 'compte_membres', 'compte_non_membres', ...
    'path_colorcell', 'source_couleur', 'id_K', 'dist', 'resume_session', 'opts_clonal', 'cellules');

%% 7. Figures
if config.faire_figures
    figures_clonales(colorcell, assembly, calcium, outline_gcampx, outline_gcampy, res, namefull);
end
% Anciens scripts couleur non relus (s'il y en a dans config.scripts_anciens_clonal)
if config.executer_scripts_anciens && ~isempty(config.scripts_anciens_clonal)
    V = struct('colorcell', colorcell, 'assembly', {assembly}, 'assemblyortho', {assembly}, ...
        'nbassemblies', nbassemblies, 'NCl', NCl, 'PCl_final', PCl_final, 'Race', Race, ...
        'Raster', pre.Raster, 'TRace', pre.TRace, 'calcium', calcium, ...
        'outline_gcampx', {outline_gcampx}, 'outline_gcampy', {outline_gcampy}, ...
        'gcamp_centroid', gcamp_centroid, 'namefull', namefull, 'path', path, ...
        'path_colorcell', path_colorcell, 'NCell', size(Race, 1), 'res', res, ...
        'id_animal', id_animal, 'id_champ', id_champ, 'id_session', id_session, 'id_K', id_K);
    executer_scripts_anciens(config.scripts_anciens_clonal, V);
end
end


function enregistrer_exclusion(fichier_clonal, statut, id_execution)
% Session terminée, mais exclue pour une raison prévue à l'avance.
res = table(); %#ok<NASGU>
save(fichier_clonal, 'statut', 'id_execution', 'res');
end
