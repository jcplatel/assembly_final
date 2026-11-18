function R = retrouver_ancien(nom_session, K, seuil, fichier_ancien, variantes)
%RETROUVER_ANCIEN  Quel changement du pipeline fait perdre des assemblées ?
%   Refait UNE session, à UN K, avec les réglages de l'ancien pipeline
%   (brainbow_analysis_2026_09_03 : SCE_clustering, kmeansoptnew,
%   kmeansoptrndnew, RACEORTHO_FDR, preprocessing_6), puis remet les
%   nouveaux réglages UN PAR UN. Aucune couleur n'est lue.
%
%   retrouver_ancien('444119_220919_plane0_sans_indice')            % K = 9, seuil 14
%   retrouver_ancien('444119_220919_plane0_sans_indice', 9, 13)     % autre seuil des SCE
%   retrouver_ancien(nom, 9, 14, 'E:\...\444119_220919\k9\results.mat')
%       compare aussi chaque variante aux assemblées de ce fichier de
%       l'ancien pipeline (variable assemblyortho ; mêmes numéros de
%       cellules, c.-à-d. même Fall.mat)
%   retrouver_ancien(nom, 9, 14, '', {'ancien', 'replicats', 'curveball'})
%       seulement quelques variantes
%
%   Variantes (chacune = 'ancien' + UN changement, sauf 'nouveau') :
%     ancien         tous les anciens réglages, seuil des SCE fixé à « seuil »
%     pretraitement  nouveau prétraitement : bad frames à -3 vrais sigma (avant
%                    ~ -1,8 sigma, ~10 % des trames), élargies et retirées du
%                    repos ; bruit de chaque cellule mesuré dans le repos
%     seuil_calcule  seuil des SCE calculé (percentile 99 des shuffles)
%     replicats      mêmes réplicats et même phase en ligne pour données et
%                    surrogates (avant : 200 avec phase en ligne pour les
%                    données, 20 sans pour les surrogates, 50 surrogates)
%     curveball      surrogates qui gardent l'activité de chaque cellule
%                    (avant : cellules tirées au hasard dans chaque SCE)
%     rang           seuil des clusters 'rang' (avant : 'tous')
%     MinSCE         clusters de moins de 3 SCE exclus (avant : un SCE seul
%                    avait la silhouette 1 et passait toujours)
%     raceortho      RACEORTHO : modèle pondéré, FDR global à FDR_Q
%                    (avant : uniforme, FDR 5 % cellule par cellule)
%     pondere        RACEORTHO : modèle pondéré seulement
%     fdr_global     RACEORTHO : FDR global à FDR_Q seulement
%     minimums       min_participations et min_activations_asm (avant : 1 et 1)
%     nouveau        tous les réglages de config_brainbow, seuil des SCE fixé
%                    à « seuil »
%   Dans l'autre sens, 'nouveau' + UN ancien réglage remis :
%     N_pretraitement  ancien prétraitement     N_permutation  anciens surrogates
%     N_tous           seuil des clusters 'tous' N_uniforme     RACEORTHO uniforme
%     N_fdr_cellule    FDR 5 % par cellule      N_raceortho    les deux derniers
%   Les valeurs « nouvelles » sont celles de config_brainbow.m.
%   Non reproduit : l'ordre de traitement des cellules actives dans plusieurs
%   clusters (RACEORTHO), et NClOK = 1 qui ne donnait aucune assemblée.
%
%   Résultats : <PathSave>\_retrouver_ancien\<variante>\<session>\k<K>\
%   (raster des assemblées comme d'habitude) ; tableau dans
%   <PathSave>\_retrouver_ancien\retrouver_ancien_<session>.csv
%   Durée : quelques minutes par variante (un seul K).

if nargin < 2 || isempty(K), K = 9; end
if numel(K) ~= 1, error('K : un seul nombre de clusters à la fois (ex. 9).'); end
if nargin < 3 || isempty(seuil), seuil = 14; end
if nargin < 4, fichier_ancien = ''; end
toutes = {'ancien', 'pretraitement', 'seuil_calcule', 'replicats', 'curveball', 'rang', ...
    'MinSCE', 'raceortho', 'pondere', 'fdr_global', 'minimums', 'nouveau', ...
    'N_pretraitement', 'N_permutation', 'N_tous', 'N_uniforme', 'N_fdr_cellule', 'N_raceortho'};
if nargin < 5 || isempty(variantes), variantes = toutes; end
variantes = cellstr(variantes);
inconnues = setdiff(variantes, toutes);
if ~isempty(inconnues)
    error('Variante(s) inconnue(s) : %s. Possibles : %s', strjoin(inconnues, ', '), strjoin(toutes, ', '));
end

config0 = config_brainbow('valide');
dossier_pipeline = fileparts(mfilename('fullpath'));
addpath(genpath(dossier_pipeline));
if isempty(gcp('nocreate')), parpool; end
etat_rng = rng;
remettre_rng = onCleanup(@() rng(etat_rng)); %#ok<NASGU>

sessions = lire_table_sessions(config0);
i = find(strcmp(sessions.nom, nom_session), 1);
if isempty(i)
    error('Session %s absente de la table (ou inclure = 0).', nom_session);
end

% Assemblées de l'ancien pipeline (facultatif)
asm_fichier = {};
if ~isempty(fichier_ancien)
    F = load(fichier_ancien, 'assemblyortho');
    asm_fichier = F.assemblyortho;
    fprintf('Ancien fichier : %d assemblées [%s]\n', numel(asm_fichier), ...
        strjoin(arrayfun(@(x) sprintf('%d', x), cellfun(@numel, asm_fichier), 'UniformOutput', false), ' '));
end

base = fullfile(config0.PathSave, '_retrouver_ancien');
lignes = cell(numel(variantes), 13);
asm_ancien = {};
dossiers_faits = {};
for r = 1:numel(variantes)
    v = variantes{r};
    cfg = appliquer_variante(config0, v, seuil);
    cfg.liste_K = K;
    cfg.PathSave = fullfile(base, v);
    cfg.executer_scripts_anciens = false;
    cfg.forcer_pretraitement = false;
    cfg.forcer_assemblees = false;
    cfg.faire_figures = true;
    S = sessions(i, :);
    S.dossier = {fullfile(cfg.PathSave, S.nom{1})};
    S.dossier_source = S.dossier;
    if ~isfolder(S.dossier{1}), mkdir(S.dossier{1}); end
    reprendre_pretraitement(S.dossier{1}, dossiers_faits, cfg.pre);

    fprintf('\n=========== Variante %s (K = %d) ===========\n', v, K);
    t0 = tic;
    try
        etape1_fonctionnel(table2struct(S), cfg);
        [ligne, asm_v] = resumer(S.dossier{1}, v, K, toc(t0));
        lignes(r, :) = ligne;
        dossiers_faits{end+1} = S.dossier{1}; %#ok<AGROW>
        if strcmp(v, 'ancien'), asm_ancien = asm_v; end
        lignes{r, 11} = recouvrement(asm_ancien, asm_v);
        lignes{r, 12} = recouvrement(asm_fichier, asm_v);
    catch err
        fprintf(2, 'ÉCHEC : %s\n', err.message);
        lignes(r, :) = {v, NaN, NaN, NaN, NaN, 0, '', '', NaN, '', '', '', ['échec : ' err.message]};
    end
end

R = cell2table(lignes, 'VariableNames', {'variante', 'seuil_SCE', 'n_SCE', 'SCE_hasard', ...
    'clusters_signif', 'n_assemblees', 'tailles', 'activations', 'duree_min', 'dossier_K', ...
    'comme_ancien', 'comme_fichier', 'statut'});
fprintf('\n=========== %s, K = %d : ancien + UN changement ; N_ = nouveau + UN ancien réglage ===========\n', ...
    nom_session, K);
colonnes = {'variante', 'seuil_SCE', 'n_SCE', 'SCE_hasard', 'clusters_signif', 'n_assemblees', ...
    'tailles', 'activations', 'comme_ancien'};
if ~isempty(asm_fichier), colonnes{end+1} = 'comme_fichier'; end
disp(R(:, colonnes))
fprintf(['clusters_signif : clusters de SCE significatifs sur %d. SCE_hasard : SCE attendus par hasard ' ...
    'au seuil utilisé.\ncomme_ancien : assemblées de la variante ''ancien'' retrouvées (recouvrement ' ...
    'de Jaccard >= 0.5) / leur nombre, et recouvrement moyen.\n'], K);
fprintf('Une variante qui perd des assemblées par rapport à ''ancien'' désigne le changement en cause.\n');
fprintf('Rasters : R.dossier_K de chaque ligne.\n');
writetable(R, fullfile(base, ['retrouver_ancien_' nom_session '.csv']));
end


function cfg = appliquer_variante(config0, v, seuil)
% 'ancien' = tous les anciens réglages ; les autres = 'ancien' + un changement
cfg = config0;
if strcmp(v, 'nouveau') || strncmp(v, 'N_', 2)
    cfg.pre.sce_regle = 'fixe';
    cfg.pre.sce_n_cells_threshold = seuil;
    switch v          % 'nouveau' + UN ancien réglage
        case 'N_pretraitement'
            cfg.pre.pretraitement_ancien = true;
        case 'N_permutation'
            cfg.asm.null_kmeans = 'permutation_colonnes';
        case 'N_tous'
            cfg.asm.seuil_clusters = 'tous';
        case 'N_uniforme'
            cfg.asm.null_cellules = 'uniforme';
        case 'N_fdr_cellule'
            cfg.asm.fdr_famille = 'par_cellule';
            cfg.asm.FDR_Q = 0.05;
            cfg.asm.FDR_Q_activations = 0.05;
        case 'N_raceortho'
            cfg.asm.null_cellules = 'uniforme';
            cfg.asm.fdr_famille = 'par_cellule';
            cfg.asm.FDR_Q = 0.05;
            cfg.asm.FDR_Q_activations = 0.05;
    end
    return
end

% --- Ancien prétraitement (opts de brainbow_analysis_2026_09_03) ---
cfg.pre.MinPeakDistancesce    = 4;
cfg.pre.MinPeakDistance       = 2;
cfg.pre.threshold_peak        = 3.09;
cfg.pre.synchronous_frames    = 2;
cfg.pre.SG_window             = 7;
cfg.pre.motion_correction     = true;
cfg.pre.colorsubstraction     = false;
cfg.pre.retirer_hyperactives  = false;
cfg.pre.sce_regle             = 'fixe';
cfg.pre.sce_n_cells_threshold = seuil;
cfg.pre.pretraitement_ancien  = true;
% --- Anciennes assemblées (SCE_clustering, RACEORTHO_FDR) ---
a = cfg.asm;
a.kmean_iter           = 200;
a.kmean_iter_surrogate = 20;
a.kmeans_surrogate     = 50;
a.kmeans_online        = true;     % données : phase en ligne
a.online_surrogate     = false;    % surrogates : sans
a.null_kmeans          = 'permutation_colonnes';
a.seuil_clusters       = 'tous';
a.percentile_seuil     = 95;
a.MinSCE               = 1;
a.NShuf                = 5000;
a.null_cellules        = 'uniforme';
a.fdr_famille          = 'par_cellule';
a.FDR_Q                = 0.05;
a.FDR_Q_activations    = 0.05;     % avant : le même FDR pour les activations
a.MIN_CELL_SIZE        = 5;
a.min_participations   = 1;
a.min_activations_asm  = 1;

% --- UN changement ---
n = config0.asm;
switch v
    case 'ancien'
    case 'pretraitement'
        cfg.pre.pretraitement_ancien = false;
    case 'seuil_calcule'
        cfg.pre.sce_regle = 'percentile';
        cfg.pre.sce_percentile = 99;
        cfg.pre.sce_n_cells_threshold = 5;          % minimum absolu
    case 'replicats'
        a.kmean_iter           = n.kmean_iter;
        a.kmean_iter_surrogate = n.kmean_iter_surrogate;
        a.kmeans_surrogate     = n.kmeans_surrogate;
        a.kmeans_online        = true;
        if isfield(n, 'kmeans_online'), a.kmeans_online = n.kmeans_online; end
        a.online_surrogate     = a.kmeans_online;     % la même pour les deux
    case 'curveball'
        a.null_kmeans = 'curveball';
    case 'rang'
        a.seuil_clusters = 'rang';
    case 'MinSCE'
        a.MinSCE = n.MinSCE;
    case 'raceortho'
        a.null_cellules = n.null_cellules;
        a.fdr_famille   = n.fdr_famille;
        a.FDR_Q         = n.FDR_Q;
        a.FDR_Q_activations = q_activations(n);
    case 'pondere'
        a.null_cellules = n.null_cellules;
    case 'fdr_global'
        a.fdr_famille   = n.fdr_famille;
        a.FDR_Q         = n.FDR_Q;
        a.FDR_Q_activations = q_activations(n);
    case 'minimums'
        a.min_participations  = n.min_participations;
        a.min_activations_asm = n.min_activations_asm;
end
cfg.asm = a;
end


function q = q_activations(a)
% FDR des activations d'un config.asm (FDR_Q s'il n'y en a pas de séparé)
q = a.FDR_Q;
if isfield(a, 'FDR_Q_activations') && ~isempty(a.FDR_Q_activations), q = a.FDR_Q_activations; end
end


function reprendre_pretraitement(dossier, dossiers_faits, opts)
% Recopie le prétraitement d'une variante déjà faite s'il a les mêmes
% paramètres (évite de le refaire)
f = fullfile(dossier, 'pretraitement.mat');
if memes_parametres(f, 'opts_utilisees', opts), return; end
for d = 1:numel(dossiers_faits)
    g = fullfile(dossiers_faits{d}, 'pretraitement.mat');
    if memes_parametres(g, 'opts_utilisees', opts)
        copyfile(g, f);
        return
    end
end
end


function [ligne, asm] = resumer(dossier, v, K, duree_s)
P = load(fullfile(dossier, 'pretraitement.mat'), 'opts', 'Race', 'diagnostic');
seuil = P.opts.sce_n_cells_threshold_calcule;
n_SCE = size(P.Race, 2);
n_hasard = NaN;
if isfield(P, 'diagnostic') && isfield(P.diagnostic, 'n_SCE_hasard'), n_hasard = P.diagnostic.n_SCE_hasard; end
A = load(fullfile(dossier, sprintf('k%d', K), 'assemblees.mat'), 'assemblyortho', 'NClOK', 'PCl_final');
asm = A.assemblyortho;
tailles = strjoin(arrayfun(@(x) sprintf('%d', x), cellfun(@numel, asm), 'UniformOutput', false), ' ');
act = '';
if ~isempty(A.PCl_final)
    act = strjoin(arrayfun(@(x) sprintf('%d', x), sum(A.PCl_final > 0, 2)', 'UniformOutput', false), ' ');
end
dossier_K = fullfile(dossier, sprintf('k%d', K));
ligne = {v, seuil, n_SCE, round(n_hasard), A.NClOK, numel(asm), tailles, act, ...
    round(duree_s / 60, 1), dossier_K, '', '', 'ok'};
fprintf('Variante %s : seuil %d, %d SCE (~%.0f par hasard), %d clusters significatifs, %d assemblées [%s] (%.1f min)\n', ...
    v, seuil, n_SCE, n_hasard, A.NClOK, numel(asm), tailles, duree_s / 60);
end


function txt = recouvrement(reference, asm)
% Assemblées de la référence retrouvées (Jaccard >= 0.5) et recouvrement moyen
txt = '';
if isempty(reference), return; end
J = zeros(1, numel(reference));
for a = 1:numel(reference)
    for b = 1:numel(asm)
        inter = numel(intersect(reference{a}, asm{b}));
        J(a) = max(J(a), inter / numel(union(reference{a}, asm{b})));
    end
end
txt = sprintf('%d/%d (%.2f)', sum(J >= 0.5), numel(reference), mean(J));
end
