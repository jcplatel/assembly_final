function opts = options_sequences(config)
%OPTIONS_SEQUENCES  Réglages de l'étape 5 (séquences pendant la course).
%   Valeurs de sequence_2026_01_27 (code de l'utilisateur), remplacées par
%   config.seq quand le champ existe (un ancien config_brainbow.m sans
%   config.seq marche donc aussi).

opts = struct();
opts.lissage_vitesse        = 20;    % trames : smoothdata(speed, 'gaussian', 20)
opts.seuil_vitesse          = 2;     % cm/s : course quand la vitesse lissée dépasse ce seuil
opts.max_lag                = 100;   % trames : décalage maximal de la corrélation croisée
opts.fraction_apprentissage = 0.5;   % courses tirées au hasard pour ranger les cellules ; les autres vérifient
opts.psth_avant             = 100;   % trames avant le début de course (figure)
opts.psth_apres             = 300;   % trames après le début de course (figure)
opts.nPerm_couleur          = 2000;  % mélanges des couleurs (test même couleur / couleurs différentes)
opts.min_paires_meme        = 5;     % paires minimum de même couleur ET de couleurs différentes dans la séquence
opts.min_courses            = 5;     % au moins 5 courses (choix de l'utilisateur, 6 nov.) : 3 pour ranger, 2 pour vérifier
% 19 nov. (choix de l'utilisateur) : validation croisée répétée (validation_kfold)
opts.validation             = 'moitie';  % 'moitie' = une coupure 50/50 (avant) ; 'kfold' = 5-fold x 20
opts.kfold_k                = 5;     % paquets de courses (80 % pour ranger, 20 % pour vérifier)
opts.kfold_repetitions      = 20;    % mélanges : 5 x 20 = 100 coupures
opts.kfold_min_courses      = 10;    % en dessous : 2 paquets (50/50), 50 mélanges
opts.kfold_nPerm            = 1000;  % permutations pour la p du rho moyen
opts.kfold_seuil_valide     = 0.5;   % cellule validée : retenue ET retrouvée dans >= 50 % des coupures
% 20 nov. (choix de l'utilisateur) : choix des cellules de la séquence
opts.selection              = 'otsu';   % 'otsu' = avant (courses collées, seuil d'Otsu) ; 'faux_debuts' = test par cellule
opts.max_lag_s              = 40;    % 'faux_debuts' : décalage maximal (s), signal continu
opts.n_faux_debuts          = 500;   % 'faux_debuts' : tirages de faux débuts (pendant la course)
opts.fdr_sequence           = 0.05;  % 'faux_debuts' : FDR sur les cellules de la session
% 20 nov. (choix de l'utilisateur) : cellules robustes = validées ET position stable
opts.cellules_tests         = 'validees';   % 'validees' = avant ; 'robustes' = validées ET position stable
opts.stabilite_nPerm        = 1000;  % mélanges des positions de test, par coupure
opts.stabilite_fdr          = 0.05;  % FDR sur les cellules validées de la session
% 20 nov. (choix de l'utilisateur) : signal utilisé
opts.signal                 = 'calcium';   % 'calcium' = F de data_recording (avant) ; 'cascade' = probabilité de spike
opts.cascade_modele         = 'GC8_EXC_10Hz_smoothing150ms_high_noise';   % 'cascade' : fichier cascade_<modèle>.mat
if isfield(config, 'seq') && isstruct(config.seq)
    noms = fieldnames(config.seq);
    for k = 1:numel(noms)
        opts.(noms{k}) = config.seq.(noms{k});
    end
end
end
