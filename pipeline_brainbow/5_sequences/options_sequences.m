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
if isfield(config, 'seq') && isstruct(config.seq)
    noms = fieldnames(config.seq);
    for k = 1:numel(noms)
        opts.(noms{k}) = config.seq.(noms{k});
    end
end
end
