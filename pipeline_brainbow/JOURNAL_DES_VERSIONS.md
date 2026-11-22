# Journal des versions du pipeline brainbow

Chaque envoi du code a un numéro (`version_pipeline`), dans le nom du zip
(`pipeline_brainbow_<version>.zip`) et dans `parametres_analyse.txt` de chaque
analyse lancée avec. Pour revenir en arrière : reprendre le zip de la version
voulue. Les versions avant le 11 nov. au soir n'avaient pas de numéro (zip
toujours au même nom) ; elles sont résumées ici de mémoire du projet.

## 2026.11.22d (22 nov.)
- Étape 4, demandé par l'utilisateur (« fais un code pour quantifier ça et une
  figure au moment de l'agrégation ») : `4_agregation/profil_assemblees.m`,
  appelé par agregation_stats_clonales3 après les variantes ; sorties dans
  `agregation\profil_assemblees\` et deux lignes dans le résumé.
  1. **Profil des assemblées biaisées** (exploratoire) : assemblées testables
     du test principal, biaisées (p_minp < seuilstat) contre les autres, sur
     des caractéristiques sans couleur (taille, compacité, étalement,
     activations, SCE de la session, membres actifs par activation, % du
     champ, cellules du champ, K où l'assemblée est vue et co-appartenance du
     consensus, fréquence des membres au repos, part de gros événements, plan)
     + cellules colorées (puissance). Wilcoxon rang-somme, seuil de
     Bonferroni indiqué. `profil_assemblees_biaisees.csv`,
     `profil_assemblees_detail.csv`, `biased_assembly_profile.png`.
  2. **Taille des événements** (CASCADE en PA discrets seulement) :
     événement = images consécutives avec PA, au repos ; petit (type ripple)
     = <= 3 PA et <= 300 ms. Distribution pour tous les événements et pour
     les participations aux SCE, et part de gros par animal. PA comme à
     l'étape 1 (images mauvaises et gros événements filtrés exclus).
     `taille_evenements_par_session.csv`, `event_size.png`.
  Vérifié sur resultats_cascade75_mpd6 : mêmes chiffres que les scripts
  séparés du 22 nov. (92,8 % de petits événements ; 87,3 % des participations
  aux SCE ; aucune caractéristique ne distingue les biaisées, p >= 0,14).

## 2026.11.22c (22 nov.)
- Détection CASCADE (PA discrets), demandé par l'utilisateur : option pour
  **retirer les gros événements** (une cellule tire 0 à 2 PA par ripple ;
  une bouffée riche en PA ou longue n'est probablement pas liée à un ripple).
  Événement = suite d'images consécutives avec >= 1 PA (session entière),
  retiré en entier si son nombre de PA dépasse `opts.cascade_max_pa_evenement`
  ou si sa durée dépasse `opts.cascade_max_duree_evenement_s` (en images :
  round(durée x fréquence) ; 0,3 s = 3 images, 0,2 s = 2 images).
  Champs facultatifs (absents = rien n'est retiré, comme avant) ; ils ne sont
  pas dans config_brainbow.m, donc les étapes 1 déjà faites restent valables.
  Nombre d'événements et de PA retirés affiché et gardé dans
  `pretraitement.mat` (`diagnostic.bouffees.gros_evenements`).
  Vérifié sur 443669_220614_plane0 : sans option, même Raster qu'avant ;
  3 PA / 0,3 s : 9 % des événements et 40 % des PA retirés.

## 2026.11.22b (22 nov.)
- Interface, tableau des analyses (« quand je choisis une analyse ça ne la
  charge pas ») : après « Choose » (ou double-clic), la session affichée à
  droite est rechargée dans la nouvelle analyse (la même session si elle y
  est, sinon la première) ; avant, la partie droite restait celle de
  l'analyse précédente tant qu'on ne recliquait pas sur une session.
  « Choose » prend aussi la ligne surlignée, et le double-clic sa propre ligne.

## 2026.11.22a (22 nov.)
- **Tableau des analyses** (demandé par l'utilisateur : « les noms des
  analyses deviennent incompréhensibles ») : `utilitaires
egistre_analyses.m`
  inventorie toutes les analyses de `...ssembly
ocues` et `...\cues`
  (principale et `analyses\<nom>` de chaque dossier d'étape 1) et donne à
  chacune un **code stable** (A001, A002... ; gardé dans
  `...ssembly
egistre_analyses.mat`, jamais renuméroté ; tableau aussi
  dans `registre_analyses.csv`). Colonnes : code, condition, date, dossier
  d'étape 1, analyse, version (v1/v2), réglages résumés, sessions faites à
  chaque étape, sessions ok, fin du dernier lancement (fini / erreurs /
  arrêté), SCE (total), assemblées, testables, % p < 0,05 (non corrigée),
  % p_minp < 0,05, nombre à q < 0,05 et q < 0,10 (test champ / globale).
  Cache : seules les analyses dont les fichiers ont changé sont relues.
- Interface : bouton **« Analyses table... »** (fenêtre avec ce tableau,
  triable ; « Choose » ou double-clic affiche l'analyse) ; le code est devant
  le nom dans le menu « Analysis » et dans le titre de la fenêtre ;
  l'analyse affichée est **gardée pour la prochaine ouverture** (setpref
  'brainbow_interface').
- Interface, correction : la version des assemblées suit l'analyse choisie.
  Avant, avec `version_assemblees = 'v2'` dans config_brainbow, une analyse
  v1 (principale, ou nom sans v2) était cherchée dans `analyses2` ou
  `<nom>_v2` : « K choice not done », couleurs cachées (ex.
  `resultats_cascade75_mpd6`). Maintenant : nom en v2 -> v2 ; autre nom -> v1 ;
  principale -> v1 si elle a déjà ses étapes 2-4, sinon la version de
  config_brainbow (lancement comme avant). « New... » ajoute `_v2` au nom
  quand config_brainbow est en v2 (même dossier que celui où le pipeline écrit).

## 2026.11.21b (21 nov.)
- **mid-p et Tarone par défaut** dans le test de l'étape 3
  (`stats_clonales_champ3`), demandés par l'utilisateur après le screening
  CASCADE : `config.clonal.midp = true`, `config.clonal.tarone = true`
  (avant : absents = ancien calcul ; mettre `false` pour le retrouver).
  - mid-p : permutations ex aequo avec l'observé comptées à moitié,
    p = (n_< + 0,5 (n_= + 1)) / (nPerm + 1) ; `p_min_atteignable`,
    `taux_nul` et `p0_sup` suivent la même règle (sinon l'excès est gonflé).
  - Tarone : pour chaque assemblée, une couleur dont le meilleur p possible
    (toutes ses cellules dans l'assemblée) est >= seuilstat est retirée du
    minimum. Aussi dans la variante 'territoire'.
  - Vérifié sur selection34 (analyse `selection34_midp_tarone_v2`) :
    6 / 32 assemblées à p_minp < 0,05 (attendu 4,5 %), 4 après FDR.
  - Pas encore dans le screening (`screening\table_minp.m`,
    `screening_agregation.m`) : il garde l'ancien calcul.

## 2026.11.21a (21 nov.)
- Deux nouvelles options du prétraitement, pour le screening des réglages
  CASCADE. Sans ces champs dans `config.pre`, rien ne change (défauts = comme
  avant), et les étapes 1 déjà faites restent à jour.
  - `config.pre.cascade_min_pa` (détection `cascade_discrets`, défaut 1) :
    une image compte comme événement si elle a au moins ce nombre de PA
    (2 : seulement les images avec >= 2 PA, 8 à 20 % des événements).
  - `config.pre.hyperactives_seuil_hz` (défaut 0,2) et
    `config.pre.hyperactives_min` (défaut 200 transients) : seuil des
    cellules hyperactives = max(durée du repos x seuil_hz, min). Avec
    CASCADE, 0,2 Hz et 200 retirent jusqu'à 2/3 des cellules ; 0,5 Hz sans
    minimum retire 1 à 13 % des cellules (environ le 95e percentile).

## 2026.11.20c (20 nov.)
- Détection `cascade_discrets`, nouvelle option
  `config.pre.cascade_reechantillonnage_Hz` (défaut : absente, rien ne
  change) : pour un modèle CASCADE entraîné à une autre fréquence (ex.
  `GC8m_EXC_15Hz_smoothing50ms_high_noise`, entraîné sur GCaMP8m), DFF0 est
  interpolé (linéaire) sur une grille à cette fréquence, CASCADE et les PA
  discrets y sont calculés, puis chaque PA est ramené à l'image d'origine la
  plus proche dans le temps (somme par image). Caches suffixés `_<Hz>Hz`
  dans le dossier suite2p. L'interpolation n'ajoute pas d'information
  temporelle. L'interface ne trouve pas ces caches : son raster d'activité
  montre alors le `Raster` de l'étape 1 (repos seulement).

## 2026.11.20b (20 nov.)
- `interface_brainbow`, **correction** du raster d'activité avec la
  détection CASCADE (`opts.detection = 'cascade_discrets'`) : il était vide,
  car il refaisait findpeaks avec `th_detection`, qui vaut NaN pour CASCADE.
  Il montre maintenant les images avec au moins 1 PA, lues dans le cache
  `cascade_discrets_<modèle>.mat` du dossier suite2p (session entière,
  course comprise) ; sans ce cache, le `Raster` de l'étape 1 (repos
  seulement). Le titre dit d'où viennent les points. Rien ne change pour
  la détection actuelle, ni pour les calculs du pipeline.

## 2026.11.20a (20 nov.)
- Étape 5, **correction** de la p du rho de chaque session (k-fold) : le
  hasard échange maintenant les cellules entre elles une seule fois pour
  toutes les coupures. Avant, les décalages étaient mélangés
  indépendamment dans chaque coupure ; comme les 100 coupures réutilisent
  les mêmes courses, ce hasard était trop étroit et la p trop petite. Le
  rho lui-même et le test par animal (Wilcoxon) ne changent pas ; le
  nombre de sessions à p < 0,05 peut baisser.
- Étape 5, nouveau : assemblées et séquence (`5_sequences/sequence_assemblees`),
  sur les mêmes cellules et le même ordre que le test de couleur (cellules
  validées) : (a) chaque assemblée du K retenu contient-elle plus de
  cellules de la séquence que le hasard (hypergéométrique, FDR 5 % par
  session) ; (b) les cellules d'une même assemblée sont-elles proches dans
  la séquence (même test que la couleur, assemblées mélangées). Colonnes
  `asm_*` dans `sequences_sessions.csv`, section dans `resume_sequences.txt`.
- Testés, pas retenus (options gardées, valeurs par défaut inchangées) :
  - `config.seq.selection = 'faux_debuts'` (avant et toujours : `'otsu'`) :
    signal continu, décalage jusqu'à 40 s, test par cellule contre des faux
    débuts tirés pendant la course, FDR. Presque aucune cellule ne passe
    (0 à 18 par session). Fonctions `decalage_cellules_test`,
    `choisir_cellules_sequence`.
  - `config.seq.cellules_tests = 'robustes'` (gardé : `'validees'`) :
    cellules validées ET de position stable entre apprentissage et test
    (FDR). Une fois le hasard corrigé, aucune cellule ne passe : la
    séquence est une tendance d'ensemble, pas un ordre fiable cellule par
    cellule. L'écart de position (`ecart_position`, `p_stabilite`) reste
    calculé dans `R.kfold`.
- `validation_kfold` : coupures en parallèle (parfor), décalages de chaque
  coupure gardés pour les permutations.
- Étape 5 à relancer (cues et nocues). Étapes 1 à 4 inchangées.

## 2026.11.19c (19 nov.)
- Étape 5 (séquences), choix de l'utilisateur : validation croisée répétée
  (`config.seq.validation = 'kfold'`, avant `'moitie'`, ancienne valeur en
  commentaire dans `config_brainbow.m`). 5-fold x 20 = 100 coupures 80 / 20 ;
  sessions de moins de 10 courses : 50 / 50 x 50. rho d'une session =
  moyenne des 100 coupures ; p par permutation (décalages de test mélangés
  entre les cellules, 1000 fois). Nouvelle fonction
  `5_sequences/validation_kfold`.
- Test de couleur dans la séquence : seulement les cellules validées
  (retenues sur l'apprentissage ET retrouvées dans le test dans >= 50 % des
  coupures), rangées par leur décalage moyen sur les courses de test (avant :
  toutes les cellules de la séquence sur toutes les courses).
- La coupure 50/50 unique reste pour la figure de chaque session
  (`rho_moitie`) ; le panneau C affiche aussi le résultat de la validation.
- Étape 5 à relancer (cues et nocues). Étapes 1 à 4 inchangées.

## 2026.11.19b (19 nov.)
- Interface : menu de condition à côté de « Overall setting » (*no cues* /
  *cues*). Il choisit les sessions de la table (`sans_indice` ou
  `avec_indice`) et met le dossier de l'étape 1 sur le dossier par défaut de
  la condition (`...\nocues\resultats_reglage1` ou
  `...\cues\resultats_reglage1`) ; les analyses de ce dossier apparaissent
  dans « Analysis ». Les lancements depuis l'interface utilisent cette
  condition. `config_brainbow.m` n'est pas modifié. Aucun résultat ne change.

## 2026.11.19a (19 nov.)
- `config.pre.detection = 'cascade_discrets'` (nouveau, facultatif ; champ
  absent = détection actuelle, rien ne change) : raster = potentiels
  d'action discrets de CASCADE (`infer_discrete_spikes`, aucun réglage,
  graine fixée `cascade_graine`), 1 à chaque image avec au moins 1 PA, au
  repos seulement. Pas de regroupement en événements. Nouvelles fonctions
  `1_fonctionnel/transient_detection_cascade_discrets`, `cascade_discrets.py`
  (en parallèle, `cascade_processus`, 24 par défaut). PA gardés dans le
  dossier suite2p (`cascade_discrets_<modele>.mat`).
- Test demandé : 34 sessions nocues, étapes 1-2, cellules hyperactives non
  retirées (`retirer_hyperactives = false` pour ce test seulement).

## 2026.11.18b (18 nov., soir)
- `config.pre.detection = 'cascade'` (nouveau, facultatif ; champ absent =
  détection actuelle, rien ne change) : événements tirés de la déconvolution
  CASCADE (Rupprecht et al. 2021), modèle `GC8_EXC_10Hz_smoothing150ms_high_noise`
  (`config.pre.cascade_modele`), sur DFF0 non lissé. Événement = suite
  d'images > `cascade_seuil_image` (0,1 potentiel / image) ; gardé si sa
  taille (potentiels estimés) >= `cascade_min_spikes` (1 par défaut ; 2-3 =
  bouffées seulement) et si son maximum est au repos.
  Nouvelles fonctions `1_fonctionnel/transient_detection_cascade`,
  `cascade_spikes` (appelle `cascade_lancer.py`). CASCADE installé dans
  `C:\Users\jcplatel\Desktop\matlab\Cascade` (environnement Python
  `venv_cascade`, TensorFlow 2.3). spike_prob gardé dans le dossier de la
  session (`cascade_<modele>.mat`), recalculé si la trace change.
- Hors pipeline (même jour) : masques suite2p corrigés pour 446178 et 446198
  (`color_maxwellV2/corriger_masques`, `charger_masques`) et fluorescences
  recalculées avec ces masques (`assembly_final/reextraire_fluo_masques`,
  ancien fichier gardé en `data_recording_avant_masquescorr.mat`).

## 2026.11.18a (18 nov.)
- `config.pre.detection` (nouveau, facultatif ; champ absent = détection
  actuelle, rien ne change) : `'arnaud'` = détection de preprocessing_Arnaud
  (code d'Arnaud Malvache, Malvache et al. 2016)
  (médiane locale + 3 IQR sur +/- 4 s, bouffées exclues, 1 s entre deux
  transients) ; `'bouffees'` = notre détection + bouffées exclues + 1 s.
  Nouvelle fonction `1_fonctionnel/transient_detection_bouffees`, appelée par
  `preprocessing_7` ; `diagnostic.bouffees`. Testé sur 4 sessions nocues :
  moins d'assemblées, non retenu.
- Nouveau dossier `6_place_cells` : `place_cells_session` (détection de
  place_cells_2024_08_28.mlx, corrigée : bins, corrélation « circulaire »
  remplacée par distance circulaire des champs et corrélation des courbes,
  tests de couleur par permutation, aussi par blocs spatiaux),
  `agregation_place_cells`, `projection_assemblees` (assemblées contre place
  cells et séquence de course).
- Nouvelles figures (`figures/`) : `comparer_seuils_transients`,
  `etendue_spatiale_SCE`, `ressemblance_SCE` (diagnostics des SCE, sans couleur).
- Sessions cues ajoutées à `sessions_brainbow.xlsx` (condition `avec_indice`) ;
  les analyses nocues (`conditions = {'sans_indice'}`) ne changent pas.

## 2026.11.15c (15 nov.)
- `figure_screening_comptes` refait en récapitulatif : par défaut, chaque
  réglage du prétraitement est montré avec SA meilleure combinaison
  d'assemblées (`'meilleure'`), dont le texte est écrit sur la ligne ;
  `'premiere'` = l'ancien défaut ; ou une ligne du résumé. Plusieurs
  screenings à la fois (`{dossier1, dossier2}`). Deuxième figure par animal
  (`figure_screening_comptes_animaux.png`) : SCE, assemblées testables,
  biaisées et excès. `screening_comptes.csv` : colonnes par animal en plus.
  Ne change aucun résultat.

## 2026.11.15b (15 nov.)
- `config_brainbow.m` rempli avec le réglage 1 du screening fdr (demandé par
  l'utilisateur) : `sce_regle` 'fdr' (avant 'fixe'), `sce_fdr` 0.2 (0.05),
  `sce_decalage_fdr_s` Inf (1), `synchronous_frames` 1 (2),
  `retirer_hyperactives` true (false), `liste_K` 4:14 (4:24),
  `asm.seuil_clusters` 'tous' ('rang'), `asm.FDR_Q` 0.10 (0.20),
  `version_assemblees` 'v2' ('v1'), `agreg.min_cellules_asm` 15 (10),
  `agreg.min_SCE_asm` 5 (20), `PathSave` ...\resultats_reglage1 (resultats).
  Les anciennes valeurs sont notées en commentaire sur chaque ligne.
  Pour revenir à l'ancien config : version 2026.11.15a (tag git).

## 2026.11.15a (15 nov.)
- `config.pre.population = 'sequence'` ou `'hors_sequence'` (nouveau,
  facultatif ; champ absent = toutes les cellules, rien ne change) : les SCE,
  leur seuil et les assemblées sont cherchés parmi les cellules de la
  séquence de course seulement, ou parmi les autres. Cellules de la séquence :
  `cellules_sequence` (nouveau, dans `1_fonctionnel`), même calcul que
  l'étape 5 (`sequence_tout`), recopié sans changement. Le raster des cellules
  écartées est mis à zéro dans `preprocessing_7` ; elles sont traitées comme
  non colorées ensuite (`couleurs_sans_doublons`, `screening_agregation`) :
  les couleurs sont mélangées parmi les cellules gardées seulement.
  `diagnostic.population` dans `pretraitement.mat`.
- `screening_reglages` : refuse de mélanger plusieurs populations dans un même
  screening (un screening par population).
- Sans `config.pre.population`, les résultats sont identiques à 2026.11.14a.
- `config_reglage1` (nouveau) : `config_brainbow` + le réglage 1 du screening fdr
  (fdr 0,2, synchronous_frames 1, hyperactives retirées, tous 95, FDR_Q 0,10,
  v2, exces, >= 15 cellules, >= 5 SCE) ; `config_reglage1('sequence')` ajoute
  la population.

## 2026.11.14a (14 nov.)
- `figure_screening_comptes` (nouveau, dans `figures`) : pour chaque réglage
  du prétraitement du screening, nombre de transients, de SCE (et attendus
  par hasard), d'assemblées (toutes, testables) et d'assemblées biaisées (et
  attendues par hasard). Assemblées : une combinaison de référence (la ligne
  du résumé, copiée) et la moyenne de toutes les combinaisons. Écrit
  `figure_screening_comptes.png` et `screening_comptes.csv`. Ne recalcule rien.

## 2026.11.13b (13 nov.)
- `quel_colorcell` (nouveau) : pour chaque session, le fichier de couleurs
  que le pipeline lit (même règle que `charger_couleurs`), sa date, sa taille
  contre iscell et les autres fichiers colorcell du dossier.
  `quel_colorcell(config, true)` : toutes les lignes de la table, même
  `inclure = 0`. Ne change aucun résultat.

## 2026.11.13a (13 nov.)
- `screening_reglages` (nouveau) : essaie beaucoup de réglages sur quelques
  animaux (446178, 444119, 444175 par défaut) et classe les combinaisons par
  excès d'assemblées biaisées. Prétraitement : un réglage changé à la fois,
  étape 1 refaite (long). Assemblées : clusters retenus, règle des membres,
  version 1 / 2, chaque K, 'exces', consensus, critères de l'étape 4, tout
  recalculé sans refaire le k-means. Mesure aussi ce que le tri donne par
  hasard (couleurs mélangées, screening refait : `p_selection`).
  Fichiers : `screening/` (`grille_screening`, `screening_session`,
  `screening_agregation`, `table_minp`, `resume_screening`),
  `figures/figure_screening`.
- `screening_reglages(config, struct('etapes', {{'plan'}}))` : ne calcule
  rien, affiche les sessions, les réglages et le nombre d'étapes 1 à faire.
- `consensus_assemblees` : accepte des assemblées données directement
  (`opts.assemblees_K`), pour le screening. Rien ne change pour le pipeline.
- Aucun changement dans les étapes 1 à 5 ni dans les résultats existants.

## 2026.11.12a (12 nov.)
- Contrôles du champ, regroupement spatial des couleurs (figure
  colour_spatial_clustering, panneaux A-C, et résumé) : une seule mesure dans
  toute la figure, le comptage des panneaux D et E. Paires de même couleur à
  moins de `config.paires.rayon_principal_um` (100 µm), observées / attendues
  après mélange des couleurs ; > 1 = regroupées. Avant : indice de distance
  (< 1 = regroupées), gardé dans le CSV et une ligne du résumé.
  `regroupement_couleurs` a un 5e argument (rayon). Étape 4 à relancer.

## 2026.11.11d (11 nov., 23 h)
- `verifier_couleurs`, contrôle 6 : l'image Brainbow (compositecalcium.tif)
  est-elle posée sur les masques suite2p ? Cherche le décalage de l'image
  qui met le plus de lumière dans les masques (`alignement_composite`) ;
  signale aussi les bords vides (grand décalage au recalage) et un canal
  calcium décalé. Le contrôle 4 (voisines de même couleur) ne voyait pas un
  recalage raté (443669_220615_plane0).
- La visionneuse et le contrôle ouvrent d'abord compositecalcium.tif
  (`trouver_composite`).

## 2026.11.11c (11 nov., soir)
- Version 2 des assemblées (`config.version_assemblees = 'v2'`) : une cellule
  reste membre de tous les clusters où elle est significative
  (`RACEORTHO_V2`, `creer_assemblees_v2`, `fichier_assemblees`, `est_v2`).
  Analyse à part, nom en `_v2`. Version 1 inchangée.
- Consensus : une paire comptée une fois par K (sans effet en version 1).
- Numéro de version (`version_pipeline`) et ce journal.

## 2026.11.11b (11 nov.)
- `comparer_analyses` : tableau de toutes les analyses déjà faites (relit
  aussi les anciennes sans bilan).
- `verifier_couleurs` : contrôle des fichiers colorcell.
- `visionneuse_composite` + bouton « Composite viewer... » de l'interface
  (canaux rouge, vert, bleu, calcium ; masques numérotés parmi iscell).
- Étape 4 : FDR de sensibilité sur les assemblées à `p_min_atteignable` <
  `config.agreg.fdr_filtre_pmin`.
- `main_paires_assemblees` : paires de même couleur dans les assemblées,
  pondérées par la rareté.

## 2026.11.09 – 10
- Interface : image du milieu = couleurs de tout le champ ; camemberts
  (champ, assemblée affichée) ; raster plus étroit.
- Relecture complète du code ; fonction interne des camemberts renommée.

## 2026.11.08
- Seuil des SCE : `sce_minimum` séparé du seuil fixe.
- Co-participation aux SCE ('sce') en 3b et 3d ; `main_clonal_tous_K` ;
  `diagnostic_SCE_assemblees` + figure ; 3e image de l'interface (une
  assemblée à la fois, p et p_minp à chaque K).

## Avant
Voir `LISEZMOI.md` (sections datées) et le document du projet
`etat_pipeline_brainbow.md`.
