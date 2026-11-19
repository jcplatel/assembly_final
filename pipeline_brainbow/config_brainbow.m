function config = config_brainbow(reglage)
%CONFIG_BRAINBOW  Tous les paramètres du pipeline brainbow, au même endroit.
%
%   config = config_brainbow             réglage d'ensemble choisi plus bas
%   config = config_brainbow('valide')   réglages validés
%   config = config_brainbow('ancien_rang')   anciens modèles nuls + règle 'rang'
%
%   À FIGER AVANT DE REGARDER LES RÉSULTATS CLONAUX. Chaque étape enregistre
%   avec ses résultats les paramètres qu'elle a utilisés.
%
%   Les noms sont ceux de tes anciens scripts (opts.*, kmean_iter, NShuf,
%   FDR_Q, MIN_CELL_SIZE, weight_*, nPerm, seuilstat, fdr_seuil, ...).

%% ---------- Réglage d'ensemble ------------------------------------------
% 'valide'      : réglages validés sur 444119 (curveball, 'rang', RACEORTHO
%                 pondéré, FDR global) : sur cette session mélangée, aucune
%                 fausse assemblée. 444119 à K = 9 : 2 assemblées (20, 13).
% 'ancien_rang' : anciens modèles nuls et ancien prétraitement + règle 'rang'
%                 (voir la fin du fichier). 444119 à K = 9 : 4 assemblées
%                 (35, 18, 11, 11), mais sur cette session mélangée, 1,3
%                 fausse assemblée (~15 cellules) en moyenne.
% Les résultats vont dans des dossiers séparés (resultats, resultats_ancien_rang).
if nargin < 1 || isempty(reglage)
    reglage = 'valide';
end
config.reglage = reglage;

%% ---------- Chemins ------------------------------------------------------
if ispc
    racine = 'E:\Data\Aurelie\analysis\Final\assembly\nocues\';
    config.fichier_liste_ancienne = 'E:\Data\Aurelie\data\nocues_data_recording_all.mat';
else
    racine = '/Users/platel/Desktop/exp/brainbow_hippocampus/analysis/nocues/';
    config.fichier_liste_ancienne = '';
end
config.PathSave         = fullfile(racine, 'resultats_reglage1');   % étape 1 : un dossier par session dedans
                                % (15 nov. : nouveau dossier pour le réglage 1 ; avant : 'resultats')
config.analyse          = '';   % '' = analyse principale (étapes 2 à 4 dans PathSave aussi) ;
                                % 'nom' = étapes 2 à 4 dans <PathSave>\analyses\nom, sur la même
                                % étape 1 (ex. tester une autre règle de K). Voir dossiers_analyse.
config.fichier_sessions = fullfile(racine, 'sessions_brainbow.xlsx');
config.conditions       = {'sans_indice'};   % conditions traitées ; {} = toutes celles de la table
config.seulement        = {};   % {} = toutes ; ex. {'444119_220919_plane0_sans_indice'} ou {'444119'}

%% ---------- Étapes à lancer ---------------------------------------------
config.faire_etape1 = true;    % prétraitement + assemblées pour chaque K (long : gardé en cache)
config.faire_etape2 = true;    % choix du K (rapide, toujours refait)
config.faire_etape3 = true;    % statistiques clonales (rapide, toujours refait)
config.faire_etape3b = false;  % synchronie des paires selon la couleur (= main_paires), sans les assemblées
config.faire_etape4 = false;   % agrégation (refusée tant qu'une session n'est pas finie)
config.faire_etape5 = false;   % séquences pendant la course, puis même couleur plus proche dans la séquence (= main_sequences)
% Version des assemblées (11 nov.) : 'v1' = une assemblée par cellule
% (RACEORTHO_FDR2 : une cellule membre de plusieurs clusters est mise dans
% celui où elle participe le plus) ; 'v2' = une cellule reste membre de tous
% les clusters où elle est significative (RACEORTHO_V2). La v2 se calcule à
% partir de l'étape 1 existante (pas de k-means refait), dans une analyse à
% part dont le nom finit par '_v2' (ajouté tout seul). Étapes 2 à 4 à refaire.
config.version_assemblees = 'v2';   % réglage 1 (15 nov.) : cellules partagées, analyse « v2 » ; avant : 'v1'
config.forcer_pretraitement = false;   % true : refait le prétraitement même s'il est à jour
config.forcer_assemblees    = false;   % true : refait les assemblées même si elles sont à jour
config.faire_figures        = true;
config.graine               = 1;       % tous les tirages aléatoires en dépendent

config.um_par_px_defaut     = 705 / 512;   % µm par pixel si la table ne le donne pas (size_pixel)

% Anciens scripts que je n'ai pas relus. Ils sont lancés dans un espace de
% travail qui contient tes variables habituelles ; une erreur ne bloque pas
% le pipeline (simple avertissement). Voir LISEZMOI.md.
% (rastercolor, graphSCE, exportdata, distance_calculation et
% export_data_brainbow sont devenus des fonctions du pipeline.)
config.executer_scripts_anciens = true;
config.scripts_anciens_session = {'raster_rastermap'};   % étape 1, une fois par session
config.scripts_anciens_K       = {};                     % étape 1, pour chaque K
config.scripts_anciens_clonal  = {};                     % étape 3

%% ---------- Étape 1a : prétraitement (= ton ancienne struct opts) -------
config.pre = struct( ...
    'MinPeakDistancesce',    4, ...
    'MinPeakDistance',       2, ...
    'threshold_peak',        3.09, ...  % en unités de bruit de chaque cellule (PAS un z-score)
    'synchronous_frames',    1, ...     % réglage 1 (15 nov.) ; avant : 2
    'sce_n_cells_threshold', 14, ...    % règle 'fixe' seulement : LE seuil
    'sce_minimum',           5, ...     % règles 'percentile' et 'fdr' : jamais moins de 5 cellules (garde-fou, 8 nov.)
    'sce_regle',             'fdr', ...  % règle du seuil des SCE : 'percentile', 'fixe' ou 'fdr' (réglage 1 ; avant : 'fixe')
    'sce_percentile',        95, ...    % 'percentile' : percentile des shuffles (avant : 90, plafonné à 15)
    'sce_fdr',               0.2, ...   % 'fdr' : part maximale de SCE attendus par hasard (réglage 1 ; avant : 0.05)
    'sce_decalage_fdr_s',    Inf, ...   % 'fdr' : décalage max du hasard (s) ; Inf = circulaire (réglage 1 ; avant : 1)
    'SG_window',             7, ...
    'motion_correction',     true, ...
    'colorsubstraction',     false, ...
    'retirer_hyperactives',  true, ...  % true : raster des cellules hyperactives mis à zéro (réglage 1 ; avant : false)
    'retirer_doublons',      true, ...  % même cellule comptée deux fois (traces au repos r > doublons_r_min et
    'doublons_r_min',        0.7, ...   % masques qui se touchent) : la plus petite ROI retirée (7 nov.)
    'retirer_tailles',       true, ...  % ROI trop petites (< roi_min_px pixels) ou trop grandes
    'roi_min_px',            4, ...     % (> roi_max_facteur x surface médiane de la session) retirées (7 nov.)
    'roi_max_facteur',       3);

%% ---------- Étape 1b : assemblées ---------------------------------------
config.liste_K = 4:14;                 % K calculés (réglage 1, 15 nov. ; avant : 4:24)
config.asm.kmean_iter       = 200;     % réplicats k-means des données (50 avant le 5 nov.)
config.asm.kmean_iter_surrogate = 100; % réplicats k-means de chaque surrogate (50 avant le 5 nov. ; temps de calcul :
                                       % c'est ce nombre x kmeans_surrogate qui compte). Peut être
                                       % plus petit que kmean_iter si diagnostic_kmeans montre que
                                       % la silhouette des surrogates ne monte plus au-delà.
config.asm.kmeans_surrogate = 100;     % surrogates (temps de calcul : c'est le poste principal)
config.asm.kmeans_online    = false;   % phase en ligne du k-means (données ET surrogates) : très lente.
                                       % 444119 : 200 réplicats + phase en ligne ou 50 : mêmes assemblées
config.asm.null_kmeans      = 'curveball';   % 'curveball' : garde l'activité de chaque cellule
                                             % 'permutation_colonnes' : ancien modèle
config.asm.MinSCE           = 3;       % un cluster de moins de MinSCE SCE n'est jamais retenu
config.asm.seuil_clusters   = 'tous';  % (réglage 1 ; avant : 'rang')
                                       % 'rang' : i-e cluster contre i-e cluster des surrogates
                                       % 'max'  : tous contre le meilleur (très strict)
                                       % 'tous' : ancien seuil (tous les clusters mélangés)
config.asm.percentile_seuil = 95;      % percentile des surrogates utilisé comme seuil (95 = 5 % d'erreur).
                                       % 99 demande au moins 500 surrogates pour être stable
config.asm.NShuf            = 5000;    % permutations de RACEORTHO (modèle 'pondere')
config.asm.null_cellules    = 'pondere';     % 'pondere' ; 'uniforme' = ancien modèle, calcul exact
config.asm.FDR_Q            = 0.10;    % (réglage 1 ; avant : 0.20) membres des assemblées (FDR global). Sur 444119, cellules
                                       % mélangées : aucun faux membre ni fausse assemblée à 0.10 ou 0.20
config.asm.FDR_Q_activations = 0.10;   % activations des assemblées dans les SCE
config.asm.fdr_famille      = 'globale';     % 'globale' (toute la matrice) ; 'par_cellule' = ancien
config.asm.MIN_CELL_SIZE    = 5;       % taille > MIN_CELL_SIZE : au moins 6 cellules, comme avant
config.asm.min_participations  = 2;    % SCE minimum d'une cellule dans le cluster de son assemblée
config.asm.min_activations_asm = 2;    % activations significatives minimum (avant : 1)

%% ---------- Étape 2 : choix du K (= analyse_assembly_results5, inchangé)
config.choixK.date_fixation           = '2026-10-30';   % règle 'exces' (voir plus bas)
% Règle du choix du K :
% 'exces' : K où les clusters significatifs dépassent le plus le hasard
%           (somme des écarts silhouette - seuil des surrogates du même rang).
%           Sur simulation : 2,2/4 et 4,8/6 assemblées retrouvées, contre 1,6/4 et
%           3,2/6 pour 'score', qui choisissait presque toujours le plus petit K
%           (21 sessions sur 38 au plus petit K sur les vraies données).
% 'score' : ancien score pondéré (poids ci-dessous)
% 'consensus' : pas de K unique (voir config.choixK.consensus plus bas)
config.choixK.regle                   = 'exces';
config.choixK.k_start                 = 4;
config.choixK.K_fixe                  = [];   % [] = K choisi par le score ; ex. 9 = même K partout (sensibilité)
% Règle 'consensus' (config.choixK.regle = 'consensus') : assemblées de
% consensus entre tous les K de liste_K (>= k_start), sans choisir de K.
% Co-appartenance de deux cellules = nb de K où elles sont dans la même
% assemblée / nb de K où l'une ou l'autre est dans une assemblée.
% Simulation (22 sessions, K 4-12) : autant d'assemblées retrouvées qu'avec
% 'exces' (2,8/4, 4,3/6, 6,0/8 contre 2,3/4, 4,8/6, 5,8/8) et aucune fausse
% assemblée (contre 1 à 1,7 par session avec 'exces'). Voir consensus_assemblees.
config.choixK.consensus.normalisation = 'union';   % 'tous' = / nb de K (assembly_core_analysis)
config.choixK.consensus.seuil         = 0.7;       % co-appartenance minimale (fraction)
config.choixK.consensus.min_K         = 2;         % paire ensemble dans au moins 2 K
config.choixK.nb_top                  = 1;
config.choixK.weight_silhouette       = 1/3;
config.choixK.weight_recall           = 1/3;
config.choixK.weight_cell_specificity = 1/6;
config.choixK.weight_sce_purity       = 1/6;
config.choixK.weight_size             = 0;
% 'brute' plutôt que 'session' : sur 10 animaux simulés (4 assemblées
% implantées, aucune couleur utilisée), 'session' a retenu K = 3 (deux
% assemblées fusionnées) 8 fois sur 10, 'brute' le bon K 10 fois sur 10.
% 'session' étire à pleine échelle des écarts minuscules (pureté 0,974
% contre 0,971). Remettre 'session' pour retrouver l'ancienne règle.
config.choixK.normalisation           = 'brute';
config.choixK.size_pct_target_min     = 10;
config.choixK.size_pct_target_max     = 15;
config.choixK.size_pct_floor          = 4;
config.choixK.size_pct_ceiling        = 30;

%% ---------- Étape 3 : statistiques clonales (= stats_clonales_champ2) --
config.clonal.nPerm          = 10000;
config.clonal.seuilstat      = 0.05;   % sert aussi de seuil (alpha) à l'agrégation
config.clonal.nmin_couleur   = 4;
config.clonal.nmin_assemblee = 3;
% Variantes {référence, permutation}. La PREMIÈRE ligne est le test principal.
% ('assemblees','globale' retiré le 30 oct. : il ne répond pas à la question posée)
% (3 nov.) contrôle spatial = même référence que le principal, seule la
% permutation change (avant : 'assemblees','blocs', qui changeait les deux)
config.clonal.variantes = { ...
    'champ', 'globale'; ...     % principal : excès d'une couleur par rapport au champ
    'champ', 'blocs'; ...       % contrôle spatial : couleurs mélangées entre cellules proches
    'champ', 'territoire'};     % contrôle « sélection locale ? » : par rapport à la surface de l'assemblée
config.clonal.marge_territoire_um = 10;   % 'territoire' : marge autour de l'enveloppe des cellules (µm)
% Variante possible, pas lancée pour l'instant : {'assemblees', 'activite'} (couleurs
% mélangées entre membres d'activité proche : contrôle du biais de détection)
config.clonal.taille_bloc_px = 100;    % côté des carrés pour 'blocs' (pixels)
config.clonal.taille_bloc_um = NaN;    % si renseigné ET um_par_px dans la table : remplace taille_bloc_px
config.clonal.nshuffle_distance = 10000;   % mélanges du test « chaque couleur est-elle regroupée ? »

%% ---------- Étape 3b : synchronie des paires selon la couleur (main_paires)
% Les voisines de même couleur sont-elles plus synchrones au repos que les
% voisines de couleurs différentes, à distance égale ? (pas d'assemblées)
config.paires.mesures         = {'sttc', 'dff', 'sce'};   % la PREMIÈRE est la mesure principale :
                                                  % 'sttc' : événements au repos (dépend peu de l'activité de chaque cellule)
                                                  % 'dff'  : corrélation du ΔF/F au repos (sensible à l'activité)
                                                  % 'sce'  : co-participation aux mêmes SCE (phi, 8 nov.)
config.paires.dt_sttc_s       = 0.2;      % fenêtre de coïncidence du STTC (s) : +/- 2 trames à 10 Hz,
                                          % comme un SCE (pics dans 3 trames consécutives)
config.paires.dt_sensibilite_s = [0.1 0.5];   % STTC refait à ces fenêtres (sensibilité)
config.paires.distances       = 'disques';  % (5 nov.) 'disques' : toutes les paires à moins de R µm
                                          % (cumulatif, plus de paires, moyenne plus stable) ;
                                          % 'tranches' : ancienne analyse par tranches (bornes_um)
config.paires.rayons_um       = [50 100 200 300 500 Inf];   % 'disques' : rayons R (µm)
config.paires.rayon_principal_um = 100;   % 'disques' : LE test principal (choix de l'utilisateur, 5 nov. ;
                                          % avant : 200, taille des clones Emx ~200 x 500 µm) ;
                                          % les autres rayons = profil, sans conclusion propre
config.paires.bornes_um       = [0 20 40 80 160 Inf];   % 'tranches' : tranches de distance entre centres (µm)
config.paires.pas_distance_um = 10;       % finesse de la correction de distance
config.paires.permutations    = {'couplage', 'globale', 'blocs'};   % la PREMIÈRE est le test principal
% 'couplage' : couleurs permutées entre cellules de couplage au réseau proche (20 niveaux).
%              Principal : sur simulation, si une couleur est plus couplée au réseau sans
%              lien entre sœurs, 'globale' donne 56 % de faux positifs, 'couplage' 10 %,
%              pour une puissance proche face à de vraies sœurs.
% 'globale'  : couleurs permutées entre toutes les cellules colorées
% 'blocs'    : couleurs permutées dans des carrés de taille_bloc_um (densité locale gardée)
% 'activite' : couleurs permutées entre cellules de taux d'événements proche (10 niveaux)
config.paires.taille_bloc_um  = 100;      % blocs spatiaux de la permutation 'blocs'
config.paires.nPerm           = 2000;
config.paires.min_paires      = 10;       % paires minimum (même couleur ET différentes) pour tester une case

%% ---------- Étape 3c : clustering hiérarchique (main_hierarchique) ------
% Exploratoire (7 nov.) : arbre des cellules (sans couleur), puis couleur de
% chaque groupe de l'arbre selon sa taille ; valeurs par défaut dans
% options_hierarchique si ces lignes manquent
config.hier.arbres         = {'sce', 'sttc', 'sttc_100ms', 'dff', 'fca'};
                               % SCE ; STTC résiduel +/- 200 ms ; +/- 100 ms ; Pearson dF/F résiduel ; FCA
config.hier.permutations   = {'globale', 'blocs'};  % couleurs mélangées dans le champ ; par blocs (contrôle)
config.hier.taille_bloc_um = 100;                   % blocs du contrôle spatial
config.hier.liaison        = 'average';             % lien moyen
config.hier.min_sce        = 3;                     % cellules dans au moins 3 SCE (mêmes cellules pour les 2 arbres)
config.hier.alpha          = 0.05;                  % groupe avec une couleur en excès : min des p hypergéométriques
config.hier.niveaux        = [5 10 20 40];         % coupes de l'arbre : groupes de 5, 10, 20, 40 cellules en moyenne
config.hier.min_colorees   = 3;                     % groupe testé : au moins 3 cellules colorées
config.hier.fdr_q          = 0.10;                  % FDR des groupes d'un niveau (comme les assemblées biaisées)
config.hier.nPerm          = 1000;                  % mélanges des couleurs sur le même arbre
config.hier.fca_sigma_s    = 10;                    % FCA : décalage des surrogates (s), comme l'article
config.hier.fca_nsurr      = 500;                   % FCA : surrogates par comparaison (5000-10000 dans l'article)

%% ---------- Étape 3d : probabilité d'être sœurs (main_soeurs) -----------
% 7 nov. : la similarité d'activité d'une paire augmente-t-elle avec sa
% probabilité d'être sœurs (estimée à partir des couleurs et des positions) ?
% Valeurs par défaut dans options_soeurs si ces lignes manquent
config.soeurs.bornes_um     = [0 25 50 75 100 150 200 300 400 600 Inf];  % tranches de distance pour pi
config.soeurs.pi_max_um     = 200;     % pi = 0 au-delà (taille maximale d'un clone)
config.soeurs.angle_um      = [25 200];   % paires utilisées pour trouver la direction des bandes
config.soeurs.carte_pas_um  = 50;      % carte 2D alignée sur les bandes (le long, en travers) : pas
config.soeurs.carte_max_um  = 500;     % carte 2D : jusqu'à
config.soeurs.min_observees = 5;       % pi d'une couleur : si moins de 5 paires de cette couleur, celui de toutes les couleurs
config.soeurs.mesures       = {'sttc', 'sce', 'coassemblee'};   % STTC au repos ; co-participation aux SCE ; même assemblée (K retenu)
config.soeurs.nPerm         = 500;     % mélanges des couleurs par session
config.soeurs.categories_pi = [0 0.05 0.1 0.2 Inf];      % figure dose-réponse
config.soeurs.nPerm_coeur   = 1000;    % cœur clonal de chaque assemblée : tirages par p
config.soeurs.k_voisins     = 5;       % p « voisines » : 5 plus proches cellules colorées de chaque membre
config.soeurs.min_colorees  = 3;       % assemblée testable : au moins 3 membres colorés
config.soeurs.fdr_q         = 0.10;    % FDR du verdict (comme les assemblées biaisées)

%% ---------- Étape 5 : séquences pendant la course (main_sequences) -------
% Méthode de sequence_2026_01_27 (code de l'utilisateur), puis distance dans
% la séquence même couleur / couleurs différentes ; valeurs par défaut dans
% options_sequences si ces lignes manquent
config.seq.lissage_vitesse        = 20;    % trames : vitesse lissée (gaussienne)
config.seq.seuil_vitesse          = 2;     % cm/s : course au-dessus
config.seq.max_lag                = 100;   % trames : décalage maximal de la corrélation croisée
config.seq.fraction_apprentissage = 0.5;   % courses tirées au hasard pour ranger ; les autres vérifient l'ordre
config.seq.psth_avant             = 100;   % trames avant le début de course (figure)
config.seq.psth_apres             = 300;   % trames après (figure)
config.seq.nPerm_couleur          = 2000;  % mélanges des couleurs (même couleur / couleurs différentes dans la séquence)
config.seq.min_paires_meme        = 5;     % paires minimum de même couleur ET de couleurs différentes
config.seq.min_courses            = 5;     % au moins 5 courses (choix de l'utilisateur, 6 nov.) : 3 pour ranger, 2 pour vérifier
config.seq.validation             = 'kfold';   % (19 nov., choix de l'utilisateur ; avant : 'moitie' = une coupure 50/50)
config.seq.kfold_k                = 5;     % 5-fold : 80 % des courses pour ranger, 20 % pour vérifier
config.seq.kfold_repetitions      = 20;    % 5 x 20 = 100 coupures
config.seq.kfold_min_courses      = 10;    % en dessous : 50/50 x 50 (100 coupures)
config.seq.kfold_nPerm            = 1000;  % permutations pour la p du rho moyen de chaque session
config.seq.kfold_seuil_valide     = 0.5;   % test de couleur : cellules retenues ET retrouvées dans >= 50 % des coupures

%% ---------- Étape 4 : agrégation (= agregation_stats_clonales2) -------
config.agreg.fdr_seuil            = 0.10;
config.agreg.seuil_non_biaise     = 0.5;
config.agreg.min_cellules_asm     = 15;      % réglage 1 (15 nov.) ; avant : 10
config.agreg.min_SCE_asm          = 5;       % réglage 1 (15 nov.) ; avant : 20
config.agreg.fdr_seuil_figure     = 0.10;
config.agreg.fdr_filtre_pmin      = 0.001;   % sensibilité (10 nov.) : FDR sur les assemblées à p_min_atteignable sous ce seuil
config.agreg.figure_non_testables = false;
config.agreg.autoriser_partiel    = false;   % true : agrège quand même, résultats marqués PROVISOIRE
config.agreg.test_champ           = false;   % B : test global du champ (non utilisé, 3 nov.)
config.agreg.analyses_exploratoires = false; % D, E : taille, compacité... (effet de puissance)

%% ---------- Réglage 'ancien_rang' -----------------------------------------
% Reproduit retrouver_ancien, variante 'rang' : ancien pipeline
% (brainbow_analysis_2026_09_03) + seuil des clusters 'rang'.
switch config.reglage
    case 'valide'
    case 'ancien_rang'
        config.PathSave = fullfile(racine, 'resultats_ancien_rang');
        config.pre.sce_regle             = 'fixe';
        config.pre.sce_n_cells_threshold = 14;
        config.pre.pretraitement_ancien  = true;    % bad frames, repos et bruit comme preprocessing_6
        config.pre.retirer_doublons      = false;   % comme l'ancien pipeline
        config.pre.retirer_tailles       = false;
        config.asm.kmean_iter            = 200;     % données, avec phase en ligne
        config.asm.kmean_iter_surrogate  = 20;      % surrogates, sans phase en ligne
        config.asm.kmeans_surrogate      = 50;
        config.asm.kmeans_online         = true;
        config.asm.online_surrogate      = false;
        config.asm.null_kmeans           = 'permutation_colonnes';
        config.asm.seuil_clusters        = 'rang';
        config.asm.percentile_seuil      = 95;
        config.asm.MinSCE                = 1;
        config.asm.null_cellules         = 'uniforme';
        config.asm.fdr_famille           = 'par_cellule';
        config.asm.FDR_Q                 = 0.05;
        config.asm.FDR_Q_activations     = 0.05;
        config.asm.MIN_CELL_SIZE         = 5;
        config.asm.min_participations    = 1;
        config.asm.min_activations_asm   = 1;
    otherwise
        error('config_brainbow : réglage ''%s'' inconnu (''valide'' ou ''ancien_rang'').', config.reglage);
end
end
