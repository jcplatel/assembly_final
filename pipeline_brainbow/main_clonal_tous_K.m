function T = main_clonal_tous_K(config)
%% =========================================================================
%  MAIN_CLONAL_TOUS_K : biais clonal de chaque assemblée, pour CHAQUE K
%  (8 nov., demandé par l'utilisateur : « donner le p initial puis aussi le
%  p_minp, et le calculer pour chaque K pour essayer de comprendre »)
%  -------------------------------------------------------------------------
%  EXPLORATOIRE. Le test de l'étape 3 reste celui du K retenu par l'étape 2
%  (choisi SANS les couleurs). Choisir après coup le K où le biais est le
%  plus fort sous-estimerait les p (on garderait le meilleur de ~20 essais).
%  Les mêmes assemblées reviennent d'un K à l'autre : les K ne sont pas des
%  répétitions indépendantes.
%
%  Pour chaque session et chaque K de config.liste_K dont l'étape 1 existe
%  (k<K>\assemblees.mat), plus le consensus (K = 0) s'il existe : même test
%  que l'étape 3, variante principale (config.clonal.variantes, 1re ligne).
%  Pour chaque assemblée :
%    p_initial : p hypergéométrique de la couleur la plus enrichie, NON
%                corrigé pour les 6 couleurs essayées ; sous H0, p_initial
%                < 0,05 arrive jusqu'à ~6 x 5 % du temps (taux_nul_initial) ;
%    p_minp    : le même, corrigé par permutation des couleurs (min-p de
%                Westfall & Young) : le test de l'étape 3 ;
%    taux_nul, taux_nul_initial : probabilité sous H0 que p_minp (resp.
%                p_initial) < 0,05 pour CETTE assemblée (attendu par hasard).
%  Par K : assemblées testables, nombre à p < 0,05 contre attendu, FDR 10 %
%  sur les assemblées testables de toutes les sessions à ce K, Wilcoxon
%  unilatéral sur les animaux (observé - attendu, p_minp ; sans correction
%  pour le nombre de K).
%
%  Lancer :   main_clonal_tous_K
%             config = config_brainbow(); config.seulement = {'444175'}; main_clonal_tous_K(config)
%  Écrit dans <PathSave>\clonal_tous_K\ : clonal_tous_K.csv (une ligne par
%  assemblée et par K), clonal_par_K.csv (une ligne par K), clonal_tous_K.mat,
%  resume_clonal_tous_K.txt, clonal_tous_K.png. Ne modifie rien d'autre.
%  =========================================================================

if nargin < 1, config = config_brainbow(); end
addpath(genpath(fileparts(mfilename('fullpath'))));
etat_rng = rng;
remettre_rng = onCleanup(@() rng(etat_rng)); %#ok<NASGU>
config = dossiers_analyse(config);
sessions = lire_table_sessions(config);
a_traiter = 1:height(sessions);
if isfield(config, 'seulement') && ~isempty(config.seulement)
    seulement = cellstr(config.seulement);
    a_traiter = find(ismember(sessions.nom, seulement) | ismember(sessions.id_animal, seulement))';
end
variante = config.clonal.variantes(1, :);
dossier = fullfile(config.PathSave, 'clonal_tous_K');
if ~isfolder(dossier), mkdir(dossier); end
fprintf('=== BIAIS CLONAL POUR CHAQUE K : %d session(s), variante %s_%s ===\n', numel(a_traiter), ...
    variante{1}, variante{2});

morceaux = {};
for k = a_traiter
    session = table2struct(sessions(k, :));
    fichier_pre = fullfile(session.dossier_source, 'pretraitement.mat');
    if ~isfile(fichier_pre), continue; end
    try
        pre = load(fichier_pre, 'Race');
        Race = double(pre.Race);
        path = avec_separateur(session.dossier_suite2p);
        [outline_gcampx, ~, gcamp_centroid] = charger_geometrie(path);
        [colorcell, ~, source_couleur] = charger_couleurs(path, session.fichier_colorcell);
        colorcell = couleurs_sans_doublons(colorcell, session);
        verifier_alignement(size(Race, 1), colorcell, outline_gcampx, session.nom);
        D = donnees_session(session, config, colorcell, gcamp_centroid, Race, source_couleur);
        n_K = 0;
        for K = [0, config.liste_K(:)']
            fa = fichier_assemblees(session.dossier, session.dossier_source, K);   % v1 ou v2 selon l'analyse
            if ~isfile(fa), continue; end
            A = load(fa, 'assemblyortho', 'PCl_final');
            if isempty(A.assemblyortho), continue; end
            D.assembly = A.assemblyortho;
            D.PCl_final = A.PCl_final;
            D.score_sequence = nan(numel(A.assemblyortho), 1);
            D.id_K = K;
            res = stats_clonales_champ3(D, config.clonal, variante{1}, variante{2}, session.graine + 10);
            morceaux{end+1, 1} = res(:, {'animal', 'session', 'K', 'assemblee', 'n_total', 'n_colorees', ...
                'testable', 'couleur_dominante', 'part_couleur_dominante', 'p_initial', 'taux_nul_initial', ...
                'p_minp', 'taux_nul', 'p_min_atteignable'}); %#ok<AGROW>
            n_K = n_K + 1;
        end
        fprintf('[%s] %d K\n', session.nom, n_K);
    catch err
        fprintf(2, '[%s] : ÉCHEC : %s\n', session.nom, err.message);
    end
end
T = table();
if isempty(morceaux)
    fprintf('Aucune assemblée trouvée (étape 1 faite ?).\n');
    return
end
T = vertcat(morceaux{:});
T.nom_session = T.session;
seuil = config.clonal.seuilstat;
q = config.agreg.fdr_seuil;

% FDR à chaque K, sur les assemblées testables de toutes les sessions
T.fdr = false(height(T), 1);
Ks = unique(T.K)';
for K = Ks
    m = T.K == K & T.testable;
    if any(m), T.fdr(m) = bh_fdr_mask(T.p_minp(m)', q)'; end
end
writetable(T, fullfile(dossier, 'clonal_tous_K.csv'));

% Une ligne par K
P = table();
for K = Ks
    t = T(T.K == K & T.testable, :);
    anim = {}; exces = []; p_w = NaN;
    if ~isempty(t)
        [G, anim] = findgroups(t.animal);
        exces = splitapply(@sum, double(t.p_minp < seuil), G) - splitapply(@sum, t.taux_nul, G);
        if numel(anim) >= 2, p_w = signrank(exces, 0, 'tail', 'right'); end
    end
    P = [P; table(K, numel(unique(T.nom_session(T.K == K))), sum(T.K == K), height(t), ...
        sum(t.p_initial < seuil), sum(t.taux_nul_initial, 'omitnan'), sum(t.p_minp < seuil), ...
        sum(t.taux_nul), sum(t.fdr), numel(anim), sum(exces > 0), p_w, median(t.n_total), ...
        'VariableNames', {'K', 'n_sessions', 'n_assemblees', 'n_testables', 'n_p_initial', ...
        'attendu_p_initial', 'n_p_minp', 'attendu_p_minp', 'n_fdr', 'n_animaux', 'n_animaux_exces', ...
        'p_wilcoxon', 'taille_mediane'})]; %#ok<AGROW>
end
writetable(P, fullfile(dossier, 'clonal_par_K.csv'));

% Résumé
L = {};
L{end+1} = sprintf(['BIAIS CLONAL POUR CHAQUE K (exploratoire ; variante %s_%s ; testables ; p < %g ; ' ...
    'K = 0 : consensus) :'], variante{1}, variante{2}, seuil);
for r = 1:height(P)
    L{end+1} = sprintf(['   K = %2d : %3d assemblées, %3d testables (taille médiane %g) ; p initial < %g : %3d ' ...
        '(attendu %5.1f) ; p_minp < %g : %3d (attendu %5.1f) ; FDR %g %% : %d ; animaux en excès %d/%d, ' ...
        'Wilcoxon p = %s'], P.K(r), P.n_assemblees(r), P.n_testables(r), P.taille_mediane(r), seuil, ...
        P.n_p_initial(r), P.attendu_p_initial(r), seuil, P.n_p_minp(r), P.attendu_p_minp(r), 100 * q, ...
        P.n_fdr(r), P.n_animaux_exces(r), P.n_animaux(r), texte_p(P.p_wilcoxon(r))); %#ok<AGROW>
end
L{end+1} = ['   Attention : ne pas choisir le K d''après ces résultats (meilleur de ~20 essais) ; les mêmes ' ...
    'assemblées reviennent d''un K à l''autre.'];
fid = fopen(fullfile(dossier, 'resume_clonal_tous_K.txt'), 'w', 'n', 'UTF-8');
if fid >= 0
    fprintf(fid, '%s\n', L{:});
    fclose(fid);
end
fprintf('\n%s\n', L{:});
save(fullfile(dossier, 'clonal_tous_K.mat'), 'T', 'P', 'variante', 'seuil', 'q');
try
    figure_clonal_tous_K(dossier);
catch err
    warning('figure_clonal_tous_K non exécutée : %s', err.message);
end
end


function D = donnees_session(session, config, colorcell, gcamp_centroid, Race, source_couleur)
% mêmes entrées que l'étape 3 (etape3_clonal), sans les assemblées
taille_bloc_px = config.clonal.taille_bloc_px;
if ~isnan(config.clonal.taille_bloc_um) && ~isnan(session.um_par_px)
    taille_bloc_px = config.clonal.taille_bloc_um / session.um_par_px;
end
um_px = session.um_par_px;
if isnan(um_px), um_px = config.um_par_px_defaut; end
marge_um = 10;
if isfield(config.clonal, 'marge_territoire_um'), marge_um = config.clonal.marge_territoire_um; end
D = struct();
D.marge_territoire_px = marge_um / um_px;
D.colorcell      = colorcell;
D.gcamp_centroid = gcamp_centroid;
D.Race           = Race;
D.id_animal      = session.id_animal;
D.id_champ       = session.id_champ;
D.id_session     = session.nom;
D.id_condition   = session.id_condition;
D.source_couleur = source_couleur;
D.taille_bloc_px = taille_bloc_px;
D.activite       = sum(Race > 0, 2);
D.p_detection    = NaN;
end


function t = texte_p(p)
if isnan(p), t = 'NA'; elseif p < 0.001, t = '< 0.001'; else, t = sprintf('%.3f', p); end
end
