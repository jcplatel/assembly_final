# Journal des versions du pipeline brainbow

Chaque envoi du code a un numéro (`version_pipeline`), dans le nom du zip
(`pipeline_brainbow_<version>.zip`) et dans `parametres_analyse.txt` de chaque
analyse lancée avec. Pour revenir en arrière : reprendre le zip de la version
voulue. Les versions avant le 11 nov. au soir n'avaient pas de numéro (zip
toujours au même nom) ; elles sont résumées ici de mémoire du projet.

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
