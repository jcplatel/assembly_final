function R = comparer_reglages(nom_session, liste_K, seuils_fixes, seuil_anciens)
%COMPARER_REGLAGES  Refait l'étape 1 et le choix du K d'UNE session avec
%   plusieurs réglages, pour voir lequel fait gagner ou perdre des assemblées.
%   Aucune couleur n'est lue. Chaque réglage a son propre dossier :
%   <PathSave>\_comparaison_reglages\<reglage>\ ; l'analyse principale
%   n'est pas touchée.
%
%   comparer_reglages('444119_220919_plane0_sans_indice')
%   comparer_reglages('444119_220919_plane0_sans_indice', 9:16)          % moins de K : plus rapide
%   comparer_reglages('444119_220919_plane0_sans_indice', [], [25 14 10], 14)
%
%   Réglages comparés :
%     actuel            config_brainbow tel quel (seuil des SCE calculé)
%     seuilNN           seuil des SCE fixé à NN cellules, pour chaque NN de
%                       seuils_fixes (défaut : 14 et 10)
%     seuilNN_anciens   seuil fixé à seuil_anciens (défaut : 14) + TOUS les
%                       anciens réglages des assemblées : seuil des clusters
%                       'tous', surrogates par permutation de colonnes,
%                       RACEORTHO uniforme avec FDR 5 % par cellule, 1
%                       participation et 1 activation significative suffisent
%   On peut aussi demander 'seuilNN_max' (seuil des clusters 'max') en
%   appelant appliquer_reglage à la main.
%   Les rasters de chaque réglage sont dans le dossier k<K retenu> indiqué.
%   Plus le seuil est bas, plus il y a de SCE et plus le k-means est long.

config0 = config_brainbow();
if nargin >= 2 && ~isempty(liste_K), config0.liste_K = liste_K; end
config0.executer_scripts_anciens = false;
config0.forcer_pretraitement = false;
config0.forcer_assemblees    = false;
dossier_pipeline = fileparts(mfilename('fullpath'));
addpath(genpath(dossier_pipeline));
if isempty(gcp('nocreate')), parpool; end

sessions = lire_table_sessions(config0);
i = find(strcmp(sessions.nom, nom_session), 1);
if isempty(i)
    error('Session %s absente de la table (ou inclure = 0).', nom_session);
end

if nargin < 3 || isempty(seuils_fixes), seuils_fixes = [14 10]; end
if nargin < 4 || isempty(seuil_anciens), seuil_anciens = 14; end
reglages = [{'actuel'}, arrayfun(@(s) sprintf('seuil%d', s), seuils_fixes, 'UniformOutput', false), ...
    {sprintf('seuil%d_anciens', seuil_anciens)}];
base = fullfile(config0.PathSave, '_comparaison_reglages');
lignes = cell(numel(reglages), 9);
for r = 1:numel(reglages)
    cfg = appliquer_reglage(config0, reglages{r});
    cfg.PathSave = fullfile(base, reglages{r});
    S = sessions(i, :);
    S.dossier = {fullfile(cfg.PathSave, S.nom{1})};
    S.dossier_source = S.dossier;
    if ~isfolder(S.dossier{1}), mkdir(S.dossier{1}); end
    fprintf('\n=========== Réglage %s ===========\n', reglages{r});
    t0 = tic;
    try
        etape1_fonctionnel(table2struct(S), cfg);
        analyse_assembly_results6(S, cfg);
        lignes(r, :) = resumer(S.dossier{1}, reglages{r}, toc(t0));
    catch err
        fprintf('ÉCHEC : %s\n', err.message);
        lignes(r, :) = {reglages{r}, NaN, NaN, NaN, ['échec : ' err.message], 0, '', NaN, ''};
    end
end

R = cell2table(lignes, 'VariableNames', {'reglage', 'seuil_SCE', 'n_SCE', 'K_retenu', 'statut', ...
    'n_assemblees', 'tailles', 'n_SCE_dans_assemblees', 'dossier_K'});
fprintf('\n=========== COMPARAISON : %s ===========\n', nom_session);
disp(R(:, 1:8))
fprintf('Rasters : dossier_K de chaque ligne (R.dossier_K).\n');
writetable(R, fullfile(base, ['comparaison_' nom_session '.csv']));
end


function cfg = appliquer_reglage(cfg, nom)
% 'actuel', 'seuilNN', 'seuilNN_anciens' ou 'seuilNN_max'
if strcmp(nom, 'actuel'), return; end
jetons = regexp(nom, '^seuil(\d+)(_anciens|_max)?$', 'tokens', 'once');
if isempty(jetons)
    error('Réglage inconnu : %s', nom);
end
cfg.pre.sce_regle = 'fixe';
cfg.pre.sce_n_cells_threshold = str2double(jetons{1});
suffixe = '';
if numel(jetons) >= 2, suffixe = jetons{2}; end
if strcmp(suffixe, '_max')
    cfg.asm.seuil_clusters = 'max';
elseif strcmp(suffixe, '_anciens')
    cfg.asm.seuil_clusters      = 'tous';
    cfg.asm.null_kmeans         = 'permutation_colonnes';
    cfg.asm.null_cellules       = 'uniforme';
    cfg.asm.fdr_famille         = 'par_cellule';
    cfg.asm.FDR_Q               = 0.05;
    cfg.asm.min_participations  = 1;
    cfg.asm.min_activations_asm = 1;
end
end


function ligne = resumer(dossier, nom, duree_s)
P = load(fullfile(dossier, 'pretraitement.mat'), 'opts', 'Race');
seuil = P.opts.sce_n_cells_threshold_calcule;
n_SCE = size(P.Race, 2);
C = load(fullfile(dossier, 'choix_K.mat'), 'best_K', 'statut_choix');
n_asm = 0; tailles = ''; n_ok = NaN; dossier_K = '';
if strcmp(C.statut_choix, 'ok')
    dossier_K = fullfile(dossier, sprintf('k%d', C.best_K));
    A = load(fullfile(dossier_K, 'assemblees.mat'), 'assemblyortho', 'NRaceOK');
    n_asm = numel(A.assemblyortho);
    tailles = strjoin(arrayfun(@(x) sprintf('%d', x), cellfun(@numel, A.assemblyortho), ...
        'UniformOutput', false), ' ');
    n_ok = A.NRaceOK;
end
ligne = {nom, seuil, n_SCE, C.best_K, C.statut_choix, n_asm, tailles, n_ok, dossier_K};
fprintf('Réglage %s : seuil %d, %d SCE, K = %g, %d assemblées [%s] (%.0f min)\n', ...
    nom, seuil, n_SCE, C.best_K, n_asm, tailles, duree_s / 60);
end
