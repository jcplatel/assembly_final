function opts = options_hierarchique(config)
%OPTIONS_HIERARCHIQUE  Réglages de l'étape 3c (clustering hiérarchique des
%   cellules, puis couleur de chaque groupe de l'arbre). Valeurs par défaut,
%   remplacées par config.hier quand le champ existe.

opts = struct();
opts.arbres         = {'sce', 'sttc', 'sttc_100ms', 'dff', 'fca'};
                                             % 'sce' = co-participation aux SCE (phi) ;
                                             % 'sttc' = STTC résiduel au repos (fenêtre +/- 200 ms) ;
                                             % 'sttc_100ms' = pareil, +/- 100 ms ;
                                             % 'dff' = Pearson du dF/F au repos, résiduel ;
                                             % 'fca' = Functional Clustering Algorithm (Feldt et al. 2009)
opts.permutations   = {'globale', 'blocs'};  % hasard : couleurs mélangées dans tout le champ (principal) ;
                                             % dans des carrés de taille_bloc_um (contrôle spatial)
opts.taille_bloc_um = 100;
opts.liaison        = 'average';             % lien moyen (comme le consensus)
opts.min_sce        = 3;                     % cellules dans au moins 3 SCE (mêmes cellules pour tous les arbres)
opts.alpha          = 0.05;                  % groupe « avec une couleur en excès » : min des p hypergéométriques < alpha
opts.niveaux        = [5 10 20 40];          % coupes de l'arbre : groupes de 5, 10, 20, 40 cellules en moyenne
opts.min_colorees   = 3;                     % groupe testé : au moins 3 cellules colorées
opts.fdr_q          = 0.10;                  % FDR des groupes d'un niveau (comme les assemblées biaisées)
opts.nPerm          = 1000;
opts.fca_sigma_s    = 10;                    % FCA : écart-type du décalage des surrogates (s ; 10 s dans l'article)
opts.fca_nsurr      = 500;                   % FCA : surrogates par comparaison (5000 à 10 000 dans l'article)
if isfield(config, 'paires') && isfield(config.paires, 'taille_bloc_um')
    opts.taille_bloc_um = config.paires.taille_bloc_um;
end
if isfield(config, 'hier') && isstruct(config.hier)
    noms = fieldnames(config.hier);
    for k = 1:numel(noms)
        opts.(noms{k}) = config.hier.(noms{k});
    end
end
end
