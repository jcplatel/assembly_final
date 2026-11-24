# Pipeline brainbow — lignes de commande (au 11 nov. 2026)

À taper dans la fenêtre de commande MATLAB, dossier `pipeline_brainbow` dans le chemin.

## 0. Préparer `config` (à mettre avant n'importe quelle commande ci-dessous)

```matlab
config = config_brainbow();                                   % réglages de config_brainbow.m
config.PathSave = 'E:\Data\Aurelie\analysis\Final\assembly\nocues\resultatsFDRSCE';   % dossier de l'étape 1
config.analyse  = 'analyse_2026_11_08consensus4-12';          % sous-dossier de analyses\ ('' = analyse principale)
config.seulement = {'444175'};                                % facultatif : un animal ou une session
```

Sans `config`, chaque commande prend `config_brainbow.m` tel quel.

## 1. Interface et pipeline principal

```matlab
interface_brainbow                    % interface graphique (étapes 1 à 4, 3b, 5, résultats)
main_brainbow                         % étapes 1 à 4 sans interface
main_brainbow(config)
```

## 2. Analyses complémentaires (pas lancées par l'interface)

```matlab
main_paires(config)                   % 3b : même couleur / couleurs différentes au repos (STTC, dF/F, SCE)
main_soeurs(config)                   % 3d : probabilité d'être sœurs (pi), cœurs clonaux
main_hierarchique(config)             % 3c : clustering hiérarchique, groupes de couleur
main_sequences(config)                % 5  : séquences pendant la course (aussi case « Step 5 » de l'interface)
main_clonal_tous_K(config)            % biais clonal à chaque K (p initial et p_minp)
main_paires_assemblees(config)        % paires de même couleur dans les assemblées (« mesurer au lieu de compter »)
```

## 3. Contrôles et diagnostics (ne modifient rien)

```matlab
verifier_couleurs(config)             % fichiers colorcell : taille, valeurs, alignement, composition
diagnostic_doublons(config)           % ROI en double, ROI trop petites / trop grandes
diagnostic_SCE_assemblees(config)     % petits SCE et activations d'assemblées (K retenu)
diagnostic_SCE_assemblees(config, 10) %   ... à K = 10
visionneuse_composite('E:\Data\Aurelie\data\nocues\444175\221121_plane1\')   % image composite + masques (aussi bouton de l'interface)
diagnostic_kmeans(dossier_session, K) % réplicats du k-means
diagnostic_membres(dossier_K)         % faux membres de RACEORTHO
```

## 4. Comparer des réglages

```matlab
comparer_analyses('E:\Data\Aurelie\analysis\Final\assembly\nocues')   % tableau de toutes les analyses déjà faites
test_seuils_SCE(config)                                 % seuil et nombre de SCE par session, sans refaire les assemblées
test_seuils_SCE([], [90 95 99], [10 14 20], 5)          %   percentiles, seuils fixes, minimum
comparer_seuils_SCE(true, {}, [95 99], 4:16)            % refait l'étape 1 sur quelques sessions et compare (long)
comparer_seuils_SCE(false)                              %   compare seulement
comparer_reglages('444119_220919_plane0_sans_indice')   % une session, plusieurs réglages
retrouver_ancien('444119_220919_plane0_sans_indice')    % ancien pipeline, un changement à la fois
```

## 5. Table des sessions

```matlab
creer_table_sessions(config, 'sans_indice')
ajouter_sessions(config)
ajouter_sessions(config, {'E:\...\data_1.mat'}, 'avec_indice')
```

## 6. Figures seules (sans recalcul)

```matlab
figure_paires_resultats(dossier, 100, 'sce')   % 3b, rayon 100 µm, mesure 'sttc' / 'dff' / 'sce'
figures_controles_champ(dossier)
figure_diagnostic_SCE(dossier)
figure_clonal_tous_K(dossier)
figure_paires_assemblees(dossier)
```

## 7. Lignes de config modifiées récemment

```matlab
config.pre.sce_regle          = 'fdr';   % seuil des SCE selon la part de hasard
config.pre.sce_fdr            = 0.5;
config.pre.sce_decalage_fdr_s = Inf;     % hasard : décalage circulaire
config.pre.sce_minimum        = 5;
config.agreg.min_SCE_asm      = 10;      % activations minimum d'une assemblée valide (étape 4)
config.agreg.min_cellules_asm = 10;
config.choixK.regle           = 'consensus';
config.choixK.consensus.seuil = 0.6;
config.liste_K                = 4:12;
config.paires.mesures         = {'sttc', 'dff', 'sce'};
config.soeurs.mesures         = {'sttc', 'sce', 'coassemblee'};
```

Changer `config.pre` ou `config.asm` refait l'étape 1 ; `config.agreg` : seulement l'étape 4.

## 8. Vérifier quel fichier MATLAB utilise

```matlab
which ROI_extraction2 -all            % le premier de la liste est celui lancé
```
