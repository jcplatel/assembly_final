function config = config_reglage1(population)
%% =========================================================================
%  CONFIG_REGLAGE1 : ton config_brainbow.m + le réglage 1 du screening fdr
%  (15 nov. : fdr0p2_sync1_hyper | clusters tous 95 | FDR_Q=0.10 | v2 |
%  exces | >= 15 cellules, >= 5 SCE)
%  -------------------------------------------------------------------------
%  Part de TON config_brainbow.m (chemins, table des sessions, tout le
%  reste) et ne remplace que les lignes du réglage 1. Ton config_brainbow.m
%  n'est pas modifié.
%
%  Utilisation :
%     config = config_reglage1();                  % toutes les cellules
%     config = config_reglage1('sequence');        % cellules de la séquence de course
%     config = config_reglage1('hors_sequence');   % les autres cellules
%  puis, comme d'habitude : main_brainbow(config), screening_reglages(config, opts), ...
%
%  À REMPLIR AVANT main_brainbow : config.PathSave (un NOUVEAU dossier, voir
%  plus bas), parce que l'étape 1 change avec ces réglages.
%  =========================================================================

config = config_brainbow();

%% Dossier de l'étape 1 : à remplir pour main_brainbow (un nouveau dossier
% par population). Pas nécessaire pour screening_reglages.
% config.PathSave = 'E:\Data\Aurelie\analysis\Final\assembly\nocues\resultats_reglage1';

%% Prétraitement
config.pre.sce_regle            = 'fdr';   % seuil des SCE selon la part de hasard
config.pre.sce_fdr              = 0.2;
config.pre.sce_decalage_fdr_s   = Inf;
config.pre.synchronous_frames   = 1;
config.pre.retirer_hyperactives = true;

%% Assemblées
config.asm.seuil_clusters   = 'tous';
config.asm.percentile_seuil = 95;
config.asm.FDR_Q            = 0.1;
config.version_assemblees   = 'v2';        % sans orthogonalisation (analyse « ..._v2 »)

%% Choix du K
config.liste_K          = 4:14;
config.choixK.regle     = 'exces';
config.choixK.K_fixe    = [];
config.choixK.k_start   = 4;

%% Étape 4 : assemblées gardées
config.agreg.min_cellules_asm = 15;
config.agreg.min_SCE_asm      = 5;

%% Population de cellules (facultatif, pipeline 2026.11.15a ou plus récent)
if nargin >= 1 && ~isempty(population) && ~strcmp(population, 'toutes')
    if ~ismember(population, {'sequence', 'hors_sequence'})
        error('config_reglage1 : population ''%s'' inconnue (''sequence'', ''hors_sequence'' ou rien).', population);
    end
    config.pre.population = population;
end
end
