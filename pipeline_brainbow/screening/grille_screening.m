function G = grille_screening(config)
%GRILLE_SCREENING  Les réglages essayés par screening_reglages (13 nov.).
%
%   C'EST LE FICHIER À MODIFIER pour ajouter ou retirer des réglages.
%   Le réglage de départ (« base ») est celui de config (config_brainbow.m).
%
%   Deux sortes de réglages :
%
%   1. PRÉTRAITEMENT (config.pre) : chaque ligne de la liste V ci-dessous
%      refait toute l'étape 1 (détection, SCE, k-means pour chaque K) sur
%      les sessions du screening : c'est LONG. Un seul réglage change à la
%      fois par rapport à la base (« un facteur à la fois »).
%
%   2. ASSEMBLÉES (rapide, calculé à partir du k-means déjà fait) : toutes
%      les combinaisons de
%        G.clusters   règle et percentile du seuil des clusters de SCE
%        G.membres    règle des membres (FDR_Q, participations...)
%        G.versions   'v1' (une assemblée par cellule) / 'v2' (cellules partagées)
%        G.choix_K    un K fixé, la règle 'exces', ou un consensus entre K
%        G.min_cellules, G.min_SCE   critères de l'étape 4 (config.agreg)
%
%   G = grille_screening(config), puis modifier G et le donner à
%   screening_reglages(config, struct('grille', G)) marche aussi.

base = config.pre;
G = struct();
G.liste_K = config.liste_K(:)';

%% ---------- 1. Prétraitement : {facteur, niveau, champs de config.pre changés}
% Une ligne identique à la base est sautée toute seule.
V = cell(0, 3);
% seuil des SCE (nombre minimal de cellules actives ensemble)
V(end+1, :) = {'sce_seuil', 'fdr0.5_inf',  struct('sce_regle', 'fdr', 'sce_fdr', 0.5, 'sce_decalage_fdr_s', Inf)};
V(end+1, :) = {'sce_seuil', 'fdr0.2_inf',  struct('sce_regle', 'fdr', 'sce_fdr', 0.2, 'sce_decalage_fdr_s', Inf)};
V(end+1, :) = {'sce_seuil', 'fdr0.05_1s',  struct('sce_regle', 'fdr', 'sce_fdr', 0.05, 'sce_decalage_fdr_s', 1)};
V(end+1, :) = {'sce_seuil', 'pct95',       struct('sce_regle', 'percentile', 'sce_percentile', 95)};
V(end+1, :) = {'sce_seuil', 'pct99',       struct('sce_regle', 'percentile', 'sce_percentile', 99)};
V(end+1, :) = {'sce_seuil', 'fixe10',      struct('sce_regle', 'fixe', 'sce_n_cells_threshold', 10)};
V(end+1, :) = {'sce_seuil', 'fixe14',      struct('sce_regle', 'fixe', 'sce_n_cells_threshold', 14)};
V(end+1, :) = {'sce_seuil', 'fixe20',      struct('sce_regle', 'fixe', 'sce_n_cells_threshold', 20)};
% seuil de détection des transients (en unités de bruit de la cellule)
V(end+1, :) = {'threshold_peak', '2.5', struct('threshold_peak', 2.5)};
V(end+1, :) = {'threshold_peak', '3.5', struct('threshold_peak', 3.5)};
V(end+1, :) = {'threshold_peak', '4',   struct('threshold_peak', 4)};
% largeur d'un SCE : cellules actives dans synchronous_frames + 1 trames
V(end+1, :) = {'synchronous_frames', '1', struct('synchronous_frames', 1)};
V(end+1, :) = {'synchronous_frames', '3', struct('synchronous_frames', 3)};
V(end+1, :) = {'synchronous_frames', '5', struct('synchronous_frames', 5)};
% distance minimale entre deux SCE (trames)
V(end+1, :) = {'MinPeakDistancesce', '2',  struct('MinPeakDistancesce', 2)};
V(end+1, :) = {'MinPeakDistancesce', '6',  struct('MinPeakDistancesce', 6)};
V(end+1, :) = {'MinPeakDistancesce', '10', struct('MinPeakDistancesce', 10)};
% distance minimale entre deux transients d'une cellule (trames)
V(end+1, :) = {'MinPeakDistance', '3', struct('MinPeakDistance', 3)};
V(end+1, :) = {'MinPeakDistance', '5', struct('MinPeakDistance', 5)};
% lissage Savitzky-Golay (trames, impair)
V(end+1, :) = {'SG_window', '5',  struct('SG_window', 5)};
V(end+1, :) = {'SG_window', '9',  struct('SG_window', 9)};
V(end+1, :) = {'SG_window', '11', struct('SG_window', 11)};
% garde-fou des règles 'percentile' et 'fdr' : jamais moins de ... cellules
V(end+1, :) = {'sce_minimum', '8', struct('sce_minimum', 8)};
% cellules hyperactives retirées
V(end+1, :) = {'retirer_hyperactives', 'oui', struct('retirer_hyperactives', true)};

G.pre = struct('nom', 'base', 'facteur', 'base', 'niveau', 'base', 'modif', struct());
deja = {base};
for i = 1:size(V, 1)
    fusionne = fusion_struct(base, V{i, 3});
    if any(cellfun(@(x) isequaln(x, fusionne), deja)), continue; end    % déjà dans la liste (ou = base)
    deja{end+1} = fusionne; %#ok<AGROW>
    G.pre(end+1) = struct('nom', nom_dossier([V{i, 1} '_' V{i, 2}]), 'facteur', V{i, 1}, ...
        'niveau', V{i, 2}, 'modif', V{i, 3});
end
% niveau de la base pour chaque facteur (affichage)
G.niveau_base = struct();
G.niveau_base.sce_seuil = texte_seuil_sce(base);
for ch = {'threshold_peak', 'synchronous_frames', 'MinPeakDistancesce', 'MinPeakDistance', 'SG_window', ...
        'sce_minimum'}
    if isfield(base, ch{1}), G.niveau_base.(ch{1}) = num2str(base.(ch{1})); end
end
G.niveau_base.retirer_hyperactives = 'non';
if isfield(base, 'retirer_hyperactives') && base.retirer_hyperactives
    G.niveau_base.retirer_hyperactives = 'oui';
end

%% ---------- 2. Clusters de SCE retenus : {règle, percentile des surrogates}
% 'rang' : i-e cluster contre i-e cluster des surrogates ; 'max' : tous
% contre le meilleur cluster des surrogates (strict) ; 'tous' : ancien seuil
G.clusters = { ...
    'rang', 90; ...
    'rang', 95; ...
    'rang', 99; ...
    'max',  95; ...
    'max',  99; ...
    'tous', 95};

%% ---------- 3. Membres des assemblées : champs de config.asm changés ------
% (seulement les champs utilisés APRÈS le k-means : FDR_Q, fdr_famille,
% min_participations, MIN_CELL_SIZE, FDR_Q_activations, min_activations_asm)
G.membres = struct('nom', {'FDR_Q=0.05', 'FDR_Q=0.10', 'FDR_Q=0.20', 'participations>=3'}, ...
    'modif', {struct('FDR_Q', 0.05), struct('FDR_Q', 0.10), struct('FDR_Q', 0.20), ...
    struct('min_participations', 3)});

%% ---------- 4. Version des assemblées -------------------------------------
G.versions = {'v1', 'v2'};

%% ---------- 5. Choix du K --------------------------------------------------
% un K fixé pour toutes les sessions ; la règle 'exces' ; des consensus
% (seuil de co-appartenance x plage de K)
vide = struct('nom', '', 'type', '', 'K', NaN, 'liste_K', [], 'seuil', NaN, 'min_K', NaN, 'normalisation', '');
G.choix_K = vide([]);
for K = G.liste_K
    c = vide; c.nom = sprintf('K=%d', K); c.type = 'fixe'; c.K = K;
    G.choix_K(end+1) = c;
end
c = vide; c.nom = 'exces'; c.type = 'exces'; c.liste_K = G.liste_K;
G.choix_K(end+1) = c;
min_K = 2; normalisation = 'union';
if isfield(config, 'choixK') && isfield(config.choixK, 'consensus')
    if isfield(config.choixK.consensus, 'min_K'), min_K = config.choixK.consensus.min_K; end
    if isfield(config.choixK.consensus, 'normalisation'), normalisation = config.choixK.consensus.normalisation; end
end
plages = {G.liste_K};
nK = numel(G.liste_K);
if nK >= 6                        % aussi la moitié basse et la moitié haute des K
    plages{end+1} = G.liste_K(1:ceil(nK / 2));
    plages{end+1} = G.liste_K(ceil(nK / 2):end);
end
for seuil = [0.5 0.6 0.7]
    for ip = 1:numel(plages)
        c = vide; c.type = 'consensus'; c.liste_K = plages{ip}; c.seuil = seuil;
        c.min_K = min_K; c.normalisation = normalisation;
        c.nom = sprintf('consensus %.1f K%d-%d', seuil, min(plages{ip}), max(plages{ip}));
        G.choix_K(end+1) = c;
    end
end

%% ---------- 6. Critères de l'étape 4 (config.agreg) -----------------------
G.min_cellules = [6 10 15];       % min_cellules_asm (6 = pas de critère : RACEORTHO garde >= 6)
G.min_SCE      = [5 10 20];       % min_SCE_asm
end


function s = fusion_struct(s, modif)
noms = fieldnames(modif);
for k = 1:numel(noms)
    s.(noms{k}) = modif.(noms{k});
end
end


function nom = nom_dossier(nom)
% nom utilisable comme nom de dossier
nom = regexprep(nom, '[^A-Za-z0-9_\-]', 'p');
end


function t = texte_seuil_sce(pre)
t = '?';
if ~isfield(pre, 'sce_regle'), return; end
switch pre.sce_regle
    case 'fixe'
        t = sprintf('fixe%g', pre.sce_n_cells_threshold);
    case 'percentile'
        t = sprintf('pct%g', pre.sce_percentile);
    case 'fdr'
        if isinf(pre.sce_decalage_fdr_s)
            t = sprintf('fdr%g_inf', pre.sce_fdr);
        else
            t = sprintf('fdr%g_%gs', pre.sce_fdr, pre.sce_decalage_fdr_s);
        end
end
end
