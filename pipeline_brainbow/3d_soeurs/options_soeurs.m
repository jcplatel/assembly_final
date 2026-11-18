function opts = options_soeurs(config)
%OPTIONS_SOEURS  Réglages de l'étape 3d (probabilité d'être sœurs, puis
%   similarité d'activité selon cette probabilité). Valeurs par défaut,
%   remplacées par config.soeurs quand le champ existe.

opts = struct();
opts.bornes_um      = [0 25 50 75 100 150 200 300 400 600 Inf];  % tranches de distance pour pi
opts.pi_max_um      = 200;          % pi = 0 au-delà (choix de l'utilisateur, 7 nov. : taille maximale d'un clone)
opts.angle_um       = [25 200];     % paires utilisées pour trouver la direction des bandes de couleur
opts.carte_pas_um   = 50;           % carte 2D alignée sur les bandes (le long, en travers) : pas
opts.carte_max_um   = 500;          % carte 2D : jusqu'à
opts.min_observees  = 5;            % pi d'une couleur dans une tranche : si moins de 5 paires de cette couleur, celui de toutes les couleurs
opts.mesures        = {'sttc', 'sce', 'coassemblee'};   % STTC au repos ; co-participation aux SCE ; même assemblée (K retenu)
opts.nPerm          = 500;          % mélanges des couleurs par session
opts.categories_pi  = [0 0.05 0.1 0.2 Inf];      % catégories de pi pour la figure dose-réponse
opts.nPerm_coeur    = 1000;         % cœur clonal de chaque assemblée : tirages par p
opts.k_voisins      = 5;            % p « voisines » : chaque membre coloré remplacé par une de ses 5 plus proches cellules colorées
opts.min_colorees   = 3;            % assemblée testable : au moins 3 membres colorés
opts.fdr_q          = 0.10;         % FDR du verdict (comme les assemblées biaisées)
if isfield(config, 'soeurs') && isstruct(config.soeurs)
    noms = fieldnames(config.soeurs);
    for k = 1:numel(noms)
        opts.(noms{k}) = config.soeurs.(noms{k});
    end
end
end
