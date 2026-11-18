# Pipeline brainbow — version MATLAB

Assemblées (Malvache et al. 2016, modifié) puis analyse clonale. Les noms de tes
variables sont conservés (`Tr1b`, `Raster`, `Race`, `TRace`, `IDX2`, `sCl`,
`NClOK`, `assemblyortho`, `PCl_final`, `colorcell`, `p_minp`, ...).

## Démarrage

1. **Fonctions externes.** Garde dans ton chemin MATLAB tes fonctions
   `color_substraction` et `lap_calculator` (inchangées). Le reste est dans
   ce dossier ; `main_brainbow` l'ajoute lui-même en tête du chemin.
2. **Test sur données simulées** (quelques minutes) : `test_pipeline_synthetique`.
   Attendu : Jaccard proche de 1 entre assemblées implantées et retrouvées,
   `p_minp` petit pour les 2 assemblées clonales de chaque animal.
   Je n'ai pas MATLAB ici : la version MATLAB n'a pas tourné. La version
   Python (même logique) a tourné sur ces mêmes données simulées. S'il y a
   une erreur au premier lancement, envoie-la-moi.
3. **Table des sessions** (une fois) :
   `config = config_brainbow(); creer_table_sessions(config, 'sans_indice')`.
   Puis vérifier l'Excel : identifiants, même `id_champ` pour un même champ de
   vue revu plusieurs jours, `inclure = 0` pour une session à écarter.
4. **Chemins** dans `config_brainbow.m`, puis `main_brainbow`.
5. **Une seule session ou un seul animal** : `config = config_brainbow();
   config.seulement = {'444119_220919_plane0_sans_indice'};` (ou `{'444119'}`)
   puis `main_brainbow(config)` ; l'agrégation est alors sautée.
6. **Nouvelles sessions** : `ajouter_sessions(config_brainbow())` (fichiers de
   ta liste matfile pas encore dans la table) ou
   `ajouter_sessions(config_brainbow(), {'E:\...\data.mat'})`, vérifier les
   lignes ajoutées, puis `main_brainbow` : seules les nouvelles sessions
   passent par l'étape 1.

## Structure

```
interface_brainbow.m       interface : lancer (tout / sélection), suivre, regarder rasters, GCaMP, couleurs
                           (valeurs par défaut = config_brainbow.m, relu à chaque ouverture ;
                           seuil des SCE : percentile, fixe ou FDR, la case à côté suit la règle)
main_brainbow.m            lance les étapes 1 à 4 (et 3b si config.faire_etape3b = true,
                           case « Step 3b: pairs » de l'interface)
main_paires.m              étape 3b seule : synchronie au repos des voisines de même couleur
                           contre couleurs différentes (sans les assemblées)
config_brainbow.m          TOUS les paramètres (à figer avant les résultats clonaux)
                           config_brainbow('valide') ou ('ancien_rang') : réglage d'ensemble
                           (aussi dans l'interface) ; résultats dans resultats\ ou
                           resultats_ancien_rang\
creer_table_sessions.m     crée la table des sessions depuis l'ancienne liste matfile
ajouter_sessions.m         ajoute de nouvelles sessions à la table (lignes existantes intactes)
comparer_reglages.m        refait UNE session avec plusieurs réglages (seuil des SCE,
                           anciens modèles nuls) pour comparer les assemblées
retrouver_ancien.m         refait UNE session à UN K avec TOUS les anciens réglages, puis
                           remet les nouveaux UN PAR UN : montre quel changement fait
                           perdre (ou gagner) des assemblées
1_fonctionnel/diagnostic_membres.m
                           sur un dossier k<K> : faux membres et fausses assemblées que
                           chaque réglage de RACEORTHO (FDR global / par cellule,
                           pondéré / uniforme) trouve dans des données mélangées
1_fonctionnel/             ÉTAPE 1, par session, tous les K, AUCUNE couleur
2_choix_K/                 ÉTAPE 2, choix du K sans couleur
3_clonal/                  ÉTAPE 3, SEULE étape qui lit les couleurs
4_agregation/              ÉTAPE 4, toutes les sessions et tous les animaux
figures/  utilitaires/  tests/
```

Résultats, dans `config.PathSave` :

```
<animal>_<session>_<champ>_<condition>/
    pretraitement.mat          étape 1 (gardé en cache)
    k9/assemblees.mat ...      étape 1, un dossier par K, + analysis.xlsx,
                               clusters.png, rasterall*.png, graphSCE.png
    choix_K.mat                étape 2
    clonal/stats_clonales.mat  étape 3 (+ figures couleur)
agregation/<variante>/         étape 4 (dont biased_assemblies_per_animal.png : un rond par
                               assemblée biaisée, de la couleur en excès)
agregation/controles_champ/    étape 4, contrôles du champ sans les assemblées (5 nov.) :
                               regroupement_couleurs.csv, detection_couleurs.csv,
                               colour_spatial_clustering.png, detection_by_colour.png
agregation/analysis_brainbow_ALL.xlsx   une ligne par session (K retenu)
analysisALL.xlsx               une ligne par session et par K
analysis_brainbow_best.xlsx, AnalysePerformanceK_taille.png, best_k_par_session.mat
log_analyse_*.txt
```
Les fichiers Excel de résumé sont réécrits entièrement à chaque lancement :
plus de lignes en double.
`analysisALL.xlsx` est mis à jour après CHAQUE session de l'étape 1 (une ligne
par session et par K), comme l'ancien exportdata ; s'il est ouvert dans
Excel, il est mis à jour après la session suivante.

Les étapes 2, 3 et 4 sont refaites à chaque lancement (rapides). L'étape 1 est
gardée en cache : elle n'est refaite que si son fichier manque ou si ses
paramètres ont changé (`config.forcer_pretraitement` / `forcer_assemblees`
pour la forcer). L'agrégation refuse de tourner tant qu'une session n'est pas
finie ; une session exclue pour une raison prévue (pas d'assemblée...) compte
comme finie et est listée.

## Correspondance avec tes anciens fichiers

| Ancien | Nouveau |
|---|---|
| brainbow_analysis_2026_09_03 | main_brainbow + config_brainbow + table des sessions |
| preprocessing_6 | preprocessing_7 |
| motion_correction_substraction | motion_correction_substraction2 |
| exportdata | exportdata2 (analysis.xlsx par K, analysisALL.xlsx à l'étape 2) |
| rastercolor | figure_raster_assemblees (étape 1 : couleurs d'assemblées, pas brainbow) |
| graphSCE | figure_graphSCE (étape 1) |
| distance_calculation | distance_calculation2 (étape 3) |
| export_data_brainbow | export_data_brainbow2 (analysis_brainbow_ALL.xlsx à l'étape 4) |
| transient_detection_shot_theta | transient_detection_shot_theta2 |
| SCE_detection4 | SCE_detection5 |
| baseline_calculation_glissante2, CovarMjc | inchangés |
| SCE_clustering (script) | SCE_clustering2 (fonction) + figure_clusters |
| kmeansoptnew / kmeansoptrndnew | kmeansoptnew2 / kmeansoptrndnew2 + silhouette_clusters + curveball |
| RACEORTHO_FDR (script) | RACEORTHO_FDR2 (fonction) + tirages_ponderes + bh_fdr_mask |
| analyse_assembly_results5parfor_taille | analyse_assembly_results6 |
| brainbowassemblies2026_05_07 | charger_geometrie + charger_couleurs + figures_clonales (l'ancien bloc de stats est supprimé) |
| stats_clonales_champ2 | stats_clonales_champ3 (fonction) + test_biais_detection |
| new_test_brainbow_2026_10_26 | étape 3 de main_brainbow (etape3_clonal) |
| agregation_stats_clonales2 | agregation_stats_clonales3 |

## Ce qui a changé (détail dans l'en-tête de chaque fichier)

**Structure**
- Tout est en fonctions : plus de variable restée d'un K ou d'une souris.
- Les couleurs ne sont lues qu'à l'étape 3 : corrige le bug de `results.mat`
  (sauvegardé avant `brainbowassemblies`, il contenait les couleurs de la session
  précédente pour le premier K de chaque session).
- Contrôle exact d'alignement : `numel(colorcell) == nombre de cellules`.
- Un fichier en échec est toujours supprimé avant le calcul : jamais d'ancien
  résultat agrégé par erreur.
- Graines fixées partout : résultats reproductibles, y compris dans les `parfor`.

**Détection**
- Bruit estimé à l'intérieur des épisodes de repos ; bad frames élargies par le
  filtre ; `speed` mise à la longueur de `F`.
- Seuil des SCE : percentile 99 (réglable), plus de plafond à 15,
  `opts.sce_n_cells_threshold` n'est plus écrasé ; nombre de faux SCE attendus
  dans `diagnostic.n_SCE_hasard` (cellules décalées à l'intérieur du repos).
- **8 nov.** : `sce_n_cells_threshold` ne sert plus qu'à la règle `'fixe'`.
  Avant, il servait aussi de minimum aux règles `'percentile'` et `'fdr'`
  (choix de Claude, jamais voulu par l'utilisateur) : avec 14 dans le config,
  15 sessions sur 26 étaient bloquées à 14 au 95e percentile. Minimum des
  règles calculées : `config.pre.sce_minimum` (5).
- Seuil des SCE, deuxième règle au choix (`config.pre.sce_regle = 'fdr'`) :
  plus petit seuil où les SCE attendus par hasard restent sous 5 % des SCE
  détectés. Pourquoi : les shuffles de la règle `'percentile'` (+/- 10
  trames) gardent les vrais SCE étalés sur 21 trames, donc le seuil monte
  avec le nombre de grands SCE de la session. Simulation (607 cellules,
  150 SCE de 40-120 cellules + 150 de 15-30) : `'percentile'` → seuil 26,
  254 SCE ; `'fdr'` +/- 1 s → 24, 264 SCE ; `'fdr'` circulaire → 19, 298 SCE ;
  sans aucun grand SCE, le seuil `'percentile'` serait 13-15.
  `diagnostic_seuil_SCE(dossier_session)` trace ces courbes pour une session
  réelle déjà prétraitée : à regarder avant de choisir.
  Sur une vraie session (411582_230320, ~10 Hz, 607 cellules) : à 10 Hz,
  +/- 1 s = +/- 10 trames, donc `'fdr'` à +/- 1 s donne 49, plus haut que la
  règle actuelle (40). Seul le hasard circulaire descend à 26, mais il ignore
  les fluctuations lentes communes, qui semblent présentes (le hasard à
  +/- 10 trames est 5 fois plus haut que ce qu'expliqueraient les vrais SCE
  étalés). Sous 25 cellules, le hasard domine quelle que soit la règle
  (fenêtre de 300 ms : ~15 cellules co-actives par hasard). Règle
  `'percentile'` conservée ; `diagnostic_seuil_SCE(dossier, 1)` teste une
  fenêtre plus courte.

**Assemblées**
- Seuil des clusters (`config.asm.seuil_clusters`) : `'rang'` par défaut, le
  i-e meilleur cluster des données contre le i-e meilleur cluster des
  surrogates (95e percentile), arrêt au premier échec. `'max'` (tous contre
  le meilleur cluster de chaque surrogate) et `'tous'` (ancien seuil) restent
  disponibles. Simulation (489 cellules, 305 SCE, 4 assemblées comme celles
  de l'ancien raster de 444119) : `'max'` retrouve 1,9 assemblée sur 4 à K = 9
  et 1,2 à K = 13 ; `'rang'` 2,5 et 2,0 ; `'tous'` 2,4 et 2,1. Sans
  assemblée, aucune règle n'en crée (0 sur 16 analyses chacune).
- Percentile de ce seuil : `config.asm.percentile_seuil` (95 par défaut,
  case « percentile » de l'interface), pour les trois règles. Garder 95 :
  `'max'` contrôle déjà l'erreur pour tous les clusters à la fois, et un 99e
  percentile sur 100 surrogates repose sur le 1er ou 2e plus haut (instable ;
  il en faudrait au moins 500).
- `FDR_Q` (case « FDR cellules » de l'interface) sert à l'étape 1 seulement, dans
  RACEORTHO : quelles cellules appartiennent à une assemblée, et quels SCE
  sont des activations significatives. Les tests clonaux ont leurs propres
  seuils (`clonal.seuilstat`, `agreg.fdr_seuil`).
- `config.asm.FDR_Q_activations` (facultatif) : FDR des activations seule ;
  absent = `FDR_Q` pour les deux. Sur 444119 (diagnostic_membres, cellules
  mélangées par curveball) : FDR par cellule 5 % = ~4 faux membres par
  cluster et 0,6 fausse assemblée par mélange (1,3 et 17,6 faux membres dans
  le 1er cluster avec l'ancien modèle uniforme) ; FDR global 10 % ou 20 % =
  aucun faux membre ni fausse assemblée.
- Mêmes options k-means pour données et surrogates ; `MinSCE` (les singletons
  avaient une silhouette de 1).
- Modèle nul `curveball` (garde l'activité de chaque cellule et la taille de
  chaque SCE) et `pondere` dans RACEORTHO.
- FDR global, minimum de participations, taille avant les cellules
  « promiscuous », cluster d'origine conservé (`cluster_origine`), `NClOK = 1`
  traité, au moins 3 activations significatives.
- Une seule définition de la silhouette.

**Choix de K, règle par défaut `'exces'` (30 oct.)** : pour chaque K, somme
sur les clusters significatifs de (silhouette du cluster - seuil des
surrogates du même rang) ; on garde le K où cette silhouette en excès du
hasard est la plus grande (à égalité, le plus petit K). L'ancien score
(`config.choixK.regle = 'score'`) choisissait le plus petit K dans 21
sessions sur 38 : recall, pureté et spécificité augmentent mécaniquement
quand K diminue, et rien ne les comparait au hasard. Simulations (489
cellules, 370 SCE, activité hétérogène ; 4 assemblées comme l'ancien raster
de 444119, ou 6 assemblées, ou aucune) : assemblées retrouvées (Jaccard >= 0,5)
'exces' 2,2/4 et 4,8/6, ancien score 1,6/4 et 3,2/6, silhouette seule 2,0/4
et 4,0/6, plus grand nombre d'assemblées 2,4/4 et 5,0/6 mais deux fois plus
de fragments ; aucune règle ne crée d'assemblée sans assemblée implantée.
`config.choixK.K_fixe` impose le même K partout (sensibilité).

**Ancien score de K** : critères, poids et normalisation inchangés ; spécificité et
pureté calculées sur `Race` (et non plus sur les trames de `RasterRace`) ;
plus aucune couleur lue.

**Clonal**
- Test principal sur la référence `'assemblees'` (couleurs mélangées entre
  membres d'assemblées), `'champ'` et la permutation par blocs spatiaux en
  sensibilité.
- `taux_nul` : hasard exact attendu pour chaque assemblée, utilisé par l'agrégation.
- Biais de détection : bruit et F0 comparés entre couleurs (`p_detection`, `p_F0`).

**Scripts annexes (rastercolor, graphSCE, exportdata, distance_calculation,
export_data_brainbow, motion_correction_substraction)**
- `motion_correction_substraction2` : le seuil « -3σ » valait en fait -1,8σ
  (`std` des seules déviations négatives = 0,6σ pour un bruit gaussien), soit
  ~3,6 % des trames marquées, ~10 % après élargissement à ±1 trame. Le seuil
  est maintenant à -3 vrais σ : attends-toi à beaucoup moins de bad frames.
  Les lignes « bad frames sans mouvement » (non utilisées, et `speed'` pouvait
  créer une matrice Nz x Nz) sont retirées ; la longueur de `corrXY` est vérifiée.
- `distance_calculation2` : le test « chaque couleur est-elle regroupée dans
  l'espace ? » prenait les positions des premières cellules du champ au lieu
  des cellules colorées (masque logique plus court que `NCell`) : les
  `pvalue_distance` d'avant sont fausses. Toutes les distances sont en µm
  (avant : assemblées en pixels, champ en µm) et `distance_assembly` est
  recalculé à chaque fois (avant, des valeurs d'un K précédent pouvaient
  rester). Les p de regroupement de chaque couleur sont dans
  `analysis_brainbow_ALL.xlsx` (`p_regroupement_R`, ...).
- `exportdata2` / `export_data_brainbow2` : plus d'ajout ligne à ligne
  (doublons à chaque relance, et une ligne par assemblée pour CHAQUE K dans
  `analysis_assemblies_pvalues.mat`). Les colonnes couleur d'avant
  (`successrate`, `color_assembly`, p hypergéo non corrigées) sont remplacées
  par le test principal : assemblées testables, biaisées, attendu par
  hasard, `p_champ`, biais de détection.
- `figure_raster_assemblees` (rastercolor) : au-delà de 17 assemblées, les
  couleurs manquantes étaient du noir ; figures fermées (elles
  s'accumulaient d'un K à l'autre).
- `figure_graphSCE` : figure fermée ; sinon identique.

## Contrôles du champ et nouvelles figures (5 nov.)

L'étape 3 écrit maintenant une ligne par cellule (`cellules` dans
`stats_clonales.mat` : couleur, position en µm, bruit, F0, transients au
repos, participation aux SCE, membre d'une assemblée, hyperactive). **Il faut
relancer l'étape 3** pour en profiter ; sinon l'étape 4 le signale et continue.

L'étape 4 (`4_agregation/controles_champ.m`) fait alors, sans les assemblées :

1. **Couleurs regroupées dans l'espace ?** (`3_clonal/regroupement_couleurs.m`)
   Par champ : distance moyenne des paires de même couleur ÷ paires de couleurs
   différentes, divisée par le même rapport après mélange des couleurs entre
   cellules colorées (`nshuffle_distance` mélanges). Indice 1 = hasard, < 1 =
   même couleur plus proche. Toutes couleurs, puis chaque couleur (couleur c
   contre les autres). Chaque session est son propre champ de vue
   (sessions moyennées par animal) ; Wilcoxon sur les animaux. Vérifié contre un calcul direct ; sur
   champs au hasard, 7 % (toutes couleurs) et 4,7 % (par couleur) de p < 0,05
   pour 5 % attendus (100 champs).
2. **La couleur change-t-elle la détection ?** Pour chaque couleur : moyenne
   de la couleur / moyenne des cellules colorées du champ, pour le bruit, F0,
   les transients au repos, la participation aux SCE, la part de membres ; et
   les cellules hyperactives (observées / attendues). Wilcoxon sur les
   animaux, non corrigé (pour un contrôle, c'est le choix prudent).

3. **Selon la distance** (panneaux D et E de `colour_spatial_clustering.png`,
   calculés par l'étape 3b et `agregation_paires`, demandé le 5 nov.) :
   paires de même couleur observées / attendues si les couleurs étaient
   mélangées dans le champ. D : disques cumulés. E : anneaux disjoints
   (0-50, 50-100, 100-200, 200-300, 300-500 µm). **Test** (`regroupement_anneaux`
   dans `agregation_paires`) : attendu = paires colorées de l'anneau x part
   des paires de même couleur dans tout le champ ; y = log((observé + 0,5) /
   (attendu + 0,5)) ; LMM `y ~ 1 + log2(distance / 100 µm) + (1|animal) +
   (1|session)`, anneaux 0-500 µm, pondéré par les paires attendues, ddl de
   Satterthwaite. Ordonnée = excès à 100 µm ; pente = changement quand la
   distance double. Toutes couleurs (test), puis chaque couleur (descriptif).
   Simulation sans regroupement (250 cellules, 400 champs) : biais du log
   < 0,5 % pour toutes couleurs ; jusqu'à -12 % à 0-50 µm pour la couleur la
   plus rare (estimation par couleur un peu prudente aux courtes distances).
   Dans `resume_paires.txt` et le titre du panneau E. Il suffit de relancer
   `agregation_paires` puis `figures_controles_champ`.

Les lignes correspondantes sont ajoutées à `resume_resultats.txt`. Le bouton
« Results figures » de l'interface montre ces figures (« field controls ») et
`biased_assemblies_per_animal.png` (une par variante). Les figures se
redessinent sans calcul : `figures_controles_champ(dossier)`,
`figure_biaisees_par_animal(dossier_variante)`.

**Compacité des assemblées (5 nov., sans couleur).** `4_agregation/compacite_assemblees.m`
(appelé une fois par l'étape 4) : pour chaque assemblée, l'étape 3 calcule
déjà son étalement (distance moyenne entre ses cellules) et celui de 1000
groupes de même taille tirés au hasard parmi les cellules actives. Rapport
observé / hasard (< 1 = plus compacte) ; écart au hasard (rapport - 1, pas
le log, biaisé par Jensen) moyenné par session puis animal ; Wilcoxon
unilatéral sur les animaux et LMM (ddl de Satterthwaite). Test par tirages au
hasard, pas hypergéométrique. `agregation\controles_champ\compacite_assemblees.csv`
et `assembly_compactness.png` (`figure_compacite`). Pas besoin de refaire
l'étape 3 : seulement l'étape 4.

**Champ de vue = session (5 nov.).** Un même plan (plane0...) revu un autre
jour n'est en général pas à la même profondeur : ce sont d'autres cellules.
`id_champ` n'est donc plus qu'une étiquette. Partout, chaque session est son
propre champ de vue : contrôles du champ, et effet aléatoire « champ » des
modèles mixtes de l'étape 4 (GLMM A2, LMM A5, figures), qui regroupait avant
les plane0 de jours différents. Refaire l'étape 4 (pas l'étape 3).

Le diagnostic de l'étape 1 garde aussi le nombre de trames à plus de 150
cellules mises à zéro (`diagnostic.n_trames_artefact`, nouvelles étapes 1
seulement).

## Étape 3b : synchronie des paires selon la couleur (`main_paires`)

Question : à distance égale, deux cellules de même couleur sont-elles plus
synchrones au repos que deux cellules de couleurs différentes ? N'utilise
que l'étape 1 (pas les assemblées) : `main_paires` après l'étape 1.

1. Activité seule, sans couleur (`similarite_paires`, `residus_paires`) :
   STTC des pics de transients au repos (mesure principale, +/- 0,2 s =
   +/- 2 trames à 10 Hz, comme un SCE ; refait à 0,1 et 0,5 s en
   sensibilité) et corrélation du ΔF/F au repos ; puis on retire l'effet de la distance et le niveau propre de
   chaque cellule.
2. Couleurs (`test_couleur_paires`) : delta = moyenne des paires de même
   couleur - moyenne des paires de couleurs différentes, par tranche de
   distance (0-20, 20-40, 40-80, 80-160, > 160 µm), par couleur et toutes
   couleurs ; p par permutation des couleurs (globale ; par blocs spatiaux
   et par niveaux de couplage au réseau en sensibilité). ROI qui se
   chevauchent exclues.
   Contrôle « couplage » : si les cellules d'une couleur sont chacune plus
   synchrones avec tout le réseau, les paires de cette couleur paraissent
   plus synchrones entre elles. Permuter les couleurs entre cellules de même
   niveau de couplage (20 niveaux) retire ce piège sur simulation (faux
   positifs de 70-100 % à 0-17 %) sans toucher à un vrai effet entre sœurs,
   même dispersées dans le champ. La correction multiplicative (a_i a_j)
   testée d'abord ne suffisait pas. C'est donc la permutation principale.
   Limite connue : quand une ou deux couleurs sont fortement clonales, les
   AUTRES couleurs ont environ 10 % de faux positifs au lieu de 5 % (leurs
   cellules ne peuvent pas occuper les places des clones : les couleurs ne
   sont plus échangeables). Un effet nettement plus faible que celui des
   couleurs les plus fortes est donc à prendre avec prudence.
3. Agrégation (`agregation_paires`) : l'animal est l'unité ; Wilcoxon sur
   les animaux et z combiné (Stouffer) ; q de Benjamini-Hochberg.
   « enrichissement » = paires de même couleur observées / attendues à
   chaque distance : > 1 si les couleurs (donc les clones) sont regroupées.

Simulation (500 cellules, 40 clones de 8 cellules regroupés, ~10 Hz) :
sans effet clonal, pas d'excès de faux positifs ; avec une entrée commune à
30 % des sœurs, l'effet est vu dans tous les champs jusqu'à 100 µm (STTC).
Piège vérifié : si une couleur est plus active (cellules R 2 fois plus
actives), la corrélation du ΔF/F donne de faux positifs (jusqu'à 5 champs
sur 6), le STTC non. D'où le STTC en mesure principale ; vérifier aussi
« Activité au repos différente selon la couleur » dans l'agrégation.

**Depuis le 5 nov. (relancer `main_paires`)** :
- distances en **disques cumulés** (`config.paires.distances = 'disques'`) :
  toutes les paires à moins de R µm, R = 50, 100, 200, 300, 500 µm et
  toutes ; **test principal à 100 µm** (`rayon_principal_um`, choix de
  l'utilisateur le 5 nov. ; d'abord 200 µm, taille des clones Emx), les
  autres rayons servent de profil. `'tranches'` reste possible.
- test principal = STTC, permutation `'couplage'`, toutes couleurs, paires à
  moins de 100 µm : Wilcoxon unilatéral sur les animaux (delta > 0) et p de
  permutation des sessions combinés. Écrit en phrase dans
  `paires_agregation\resume_paires.txt` (aussi affiché par « Results summary »
  de l'interface), avec le profil selon la distance, chaque couleur, la
  sensibilité (autres permutations et mesures) et le regroupement spatial.
- **un test sur toutes les distances** (idée de l'utilisateur, ajoutée après
  le test à 200 µm, donc à présenter comme complément) : chaque session donne
  aussi des **anneaux disjoints** (0-50, 50-100, 100-200, 200-300, 300-500,
  > 500 µm, obtenus par différence des disques, sans calcul en plus) ; modèle
  linéaire mixte sur tous les anneaux de toutes les sessions :
  `delta ~ 1 + (1|animal) + (1|session)` (delta = même couleur - couleurs
  différentes, sans correction du hasard), pondéré par le nombre de paires,
  ddl de Satterthwaite (simulation, 8 animaux : 3,7 % de faux positifs, 77 %
  de puissance pour un effet de 0,0004 ; le test par défaut donnait 7 à 12 %
  de faux positifs). Aussi : pente selon log2(distance), Wilcoxon sur les
  animaux (moyenne de toutes les distances). Dans `resume_paires.txt` et le
  titre du panneau B.
- **valeurs brutes (5 nov.)** : `paires_brutes_session` calcule pour chaque
  session la moyenne BRUTE (sans résidus : ni distance ni niveau des
  cellules retirés) des paires de même couleur et de couleurs différentes,
  pour chaque mesure et distance (rapide, sans permutation ; gardé dans
  `<session>\paires\paires_brutes.mat`). `agregation_paires` les regroupe
  (Wilcoxon sur les animaux, `paires_brutes_resultats.csv`, une ligne du
  résumé) ; panneau A1b de la figure et 2e ligne de `figure_paires_meme_diff`.
- tests sur les animaux et modèle mixte : différence même couleur - couleurs
  différentes telle quelle (choix de l'utilisateur, 5 nov.) ; les
  permutations ne servent qu'aux p de chaque session et à leur combinaison.
- chaque session garde aussi les moyennes des paires de même couleur et de
  couleurs différentes (`moyenne_meme`, `moyenne_diff`) et la similarité
  brute moyenne (delta exprimé en %).
- autre rayon ou autre mesure, sans rien recalculer (exploration ; le test
  principal reste celui de `config.paires`) :
  `figure_paires_resultats(dossier, 100)` (paires à moins de 100 µm),
  `figure_paires_resultats(dossier, 100, 'dff')` (Pearson du dF/F),
  `figure_paires_meme_diff(dossier, 100)` : « même couleur contre couleurs
  différentes » en STTC et en Pearson côte à côte (`paires_meme_diff_100um.png`).
- nouvelle figure `paires_agregation.png` (`figures/figure_paires_resultats`,
  redessinable sans calcul) : A. même couleur contre couleurs différentes,
  une ligne par animal (paires à moins de 200 µm) ; B. différence selon le
  rayon (un trait par animal, moyenne ± sem, p) ; C. chaque couleur ; D.
  paires de même couleur observées / attendues selon la distance.

## Étape 3c : clustering hiérarchique (`main_hierarchique`, 7 nov., exploratoire)

Question : en partant de petits groupes de cellules qui s'activent ensemble
et en montant vers des groupes de plus en plus gros, jusqu'à quelle taille
trouve-t-on des groupes avec une couleur en excès ? Lancer `main_hierarchique`
(pas de case dans l'interface). N'utilise ni les assemblées ni le K.

- **Cellules** : actives dans au moins 3 SCE (`config.hier.min_sce`), les
  mêmes pour les deux arbres.
- **Cinq arbres, sans couleur** (choix de l'utilisateur : comparés) :
  - `sce` : corrélation des lignes de Race (phi : corrige le niveau
    d'activité de chaque cellule), distance = 1 − phi, lien moyen ;
  - `sttc` : résidus du STTC au repos de l'étape 3b (distance et niveau de
    chaque cellule retirés), distance = max − résidu ; paire non mesurée =
    résidu 0 ; lien moyen ;
  - `sttc_100ms` (7 nov.) : pareil avec une fenêtre de coïncidence de
    ± 100 ms (± 1 trame à 10 Hz) au lieu de ± 200 ms ;
  - `dff` (7 nov.) : corrélation de Pearson du dF/F au repos, résidus comme à
    l'étape 3b (distance, niveau de chaque cellule et couplage au signal
    commun retirés), distance = max − résidu, lien moyen ;
  - `fca` (ajouté le 7 nov., comparé aux autres) : Functional Clustering
    Algorithm de Feldt et al. 2009 (Phys Rev E 80:056104), utilisé par Feldt
    Muldoon et al. 2013 (PNAS) (`fca_arbre.m`). Événements au repos ;
    similarité = distance moyenne au plus proche événement (AMD, normalisée) ;
    surrogates : chaque événement décalé (σ = 10 s comme l'article,
    `config.hier.fca_sigma_s`) ; on fusionne la paire la plus significative
    (au-delà de 95 % des surrogates), les deux trains deviennent un seul ;
    coupe quand plus rien n'est significatif ; clusters = groupes d'au moins
    3 cellules formés avant la coupe. Écarts à l'article pour le temps de
    calcul : 500 surrogates (5000-10 000 dans l'article), seule la ligne du
    nouveau groupe est recalculée, chaque sens de l'AMD décale le train de
    départ. L'arbre est continué jusqu'à la racine pour la comparaison par
    taille ; les clusters seuls sont aussi testés (`fca_amas`).
    **Simulation (cellules indépendantes, sans aucun groupe)** : 40 cellules →
    19 fusions « significatives » et 5 clusters ; 150 cellules → 114 fusions
    et 12 clusters (jusqu'à 25 cellules). L'article n'applique aucune
    correction pour le grand nombre de paires testées, et avec 300 à 500
    cellules la meilleure paire passe presque toujours le seuil de 95 % par
    hasard. Le test des couleurs reste valable (couleurs mélangées sur les
    mêmes clusters), mais les clusters FCA ne sont pas tous de vrais groupes.
    Pas de correction ajoutée (à décider avec l'utilisateur).
- **Coupes de l'arbre (7 nov., demandé par l'utilisateur : peu de niveaux
  pour ne pas multiplier les tests ; remplace les classes de taille)** :
  chaque arbre est coupé en groupes séparés de ~5, ~10, ~20 et ~40 cellules
  en moyenne (nombre de groupes = cellules ÷ niveau, `cluster(Z,'maxclust',k)`),
  fixés avant de voir les résultats. Les arbres étant peu structurés, une coupe
  donne quelques gros groupes et beaucoup de cellules seules ; seuls les
  groupes d'au moins 3 cellules colorées sont testés (`test_couleur_groupes`).
  - Par groupe : couleur la plus en excès (min des 6 p hypergéométriques).
  - **Un test par coupe** : nombre de groupes avec p < 0,05, contre la même
    chose sur la même coupe avec les couleurs mélangées (1000 fois) dans le
    champ (principal) et par blocs de 100 µm (contrôle). Le mélange tient
    compte des 6 couleurs et du nombre de groupes.
  - **Quels groupes** : p_minP de chaque groupe (choix de la meilleure couleur
    corrigé par les mêmes mélanges, comme les assemblées), puis FDR 10 %
    (Benjamini-Hochberg) sur les groupes de la coupe.
  - Clusters FCA : testés de la même façon, comme une coupe de plus (niveau 0).
- **Animaux** : somme de leurs sessions, % de groupes testés avec une couleur
  en excès, observé et attendu ; Wilcoxon unilatéral sur les animaux par
  coupe (4 tests par arbre), sans correction ; nombre de groupes
  significatifs après FDR. Simulation (comptes de Poisson, 8 animaux) : ce
  Wilcoxon est prudent quand peu de groupes sont attendus par animal.
- Sorties : `<session>\hierarchique\` (`hierarchique.mat`, `arbres.png` :
  l'arbre, la couleur de chaque cellule sous les feuilles, les 4 coupes en
  pointillé, un rond sur chaque groupe significatif après FDR, et observés /
  attendus par coupe) ;
  `<PathSave>\hierarchique_agregation\` (`hierarchique_sessions.csv`,
  `hierarchique_animaux.csv`, `resume_hierarchique.txt`, `hierarchique.png`).

## Étape 3d : probabilité d'être sœurs (`main_soeurs`, 7 nov.)

Question : **la similarité d'activité d'une paire augmente-t-elle avec sa
probabilité d'être sœurs ?** (proposé et accepté le 7 nov. pour sortir des
tests par assemblée ou par groupe, qui multiplient les corrections.)

Pourquoi : la couleur n'est pas le clone. Parmi les paires de même couleur,
seule une partie sont de vraies sœurs (beaucoup moins pour une couleur
fréquente, presque toutes pour une couleur rare à courte distance). Toute
comparaison « même couleur contre couleur différente » est donc diluée.

1. **pi** (`carte_soeurs`, couleurs et positions seulement, toutes les
   sessions) : pour la couleur c et la tranche de distance b, observé /
   attendu des paires de même couleur ; pi = 1 − 1 / rapport (0 si ≤ 1).
   **Attendu** (7 nov., 18 h 15, accepté par l'utilisateur) = paires de
   couleurs différentes de la tranche × (paires de même couleur c / paires
   de couleurs différentes au-delà de 200 µm), session par session. Les
   paires de couleurs différentes ne sont jamais sœurs et se répartissent
   comme les paires de même couleur non sœurs ; au-delà de 200 µm il n'y a
   pas de sœurs. L'ancien attendu (part des paires de même couleur dans tout
   le champ) contenait les sœurs elles-mêmes et sous-estimait pi, surtout
   pour les couleurs rares (question de l'utilisateur : « une couleur rare
   comme le bleu, à moins de 200 µm, quasi 100 % sœurs »). Simulation
   (clones 200 × 500 µm, champs de 600 µm) :

   | sœurs parmi les paires de même couleur | 0-25 µm | 50-75 | 100-150 | 150-200 |
   |---|---|---|---|---|
   | couleur rare (1 clone par champ), vrai | 100 % | 100 % | 100 % | 100 % |
   | ancien calcul | 86 % | 80 % | 61 % | 26 % |
   | nouveau calcul | 98 % | 97 % | 94 % | 87 % |
   | couleur fréquente, vrai | 55 % | 48 % | 30 % | 16 % |
   | ancien calcul | 31 % | 27 % | 11 % | 3 % |
   | nouveau calcul | 52 % | 45 % | 25 % | 13 % |

   Repli sur pi de toutes les couleurs seulement si moins de 5 paires de
   cette couleur sont **observées** dans la tranche (`min_observees`) ; un
   attendu < 5 est normal pour une couleur rare (l'ancienne règle « attendu
   < 5 » aurait remplacé 98 % par ~40 %). Hypothèses : un clone a une seule
   couleur (25 % d'erreurs de classement simulées : couleur fréquente 32 %
   pour 51 % vrai) ; pas de sœurs au-delà de 200 µm (clones de 500 µm de
   long : pi un peu trop bas) ; tout l'excès est clonal (un gradient
   technique de couleur le gonflerait). Les proportions de couleur et le
   niveau du hasard sont propres à chaque session ; pi est ensuite commun à
   toutes les sessions.
   **pi = 0 au-delà de 200 µm** (`pi_max_um`, choix de l'utilisateur le
   7 nov. : taille maximale d'un clone). Sur 444119, les paires G à plus de
   400 µm étaient en excès (pi 23-61 %), ce qui ne peut pas venir de sœurs ;
   la figure montre encore ces valeurs, en zone grise.
   **Direction des bandes de couleur** (7 nov. : les bandes n'ont pas la même
   direction d'une session à l'autre) : dans chaque session, direction où
   les paires de même couleur à 25-200 µm sont en excès (orientation de
   chaque paire, forme du champ retirée). Pour la carte « le long / en
   travers des bandes » et pour le contrôle, les paires de la couleur c sont
   alignées sur l'angle trouvé **sans aucune paire touchant la couleur c**
   (d'autres clones) : sinon on s'aligne sur son propre bruit et on voit un
   faux allongement (+9 % en simulation avec des clones ronds). Simulation
   (40 champs) : clones allongés à 30° → angle retrouvé à 7° près (médiane),
   le long 1,13 / en travers 0,99, 95 % des champs le long > en travers ;
   clones ronds → le long = en travers (−0,01 ± 0,01), 50 % des champs. Un p
   par mélange des couleurs a été écarté (70 % de « direction nette » avec
   des clones ronds). Figure `angles_bandes_<animal>.png` pour vérifier
   l'angle à l'œil, tableau `angles_bandes.csv`. La carte est descriptive :
   pi reste calculé sur la simple distance.
2. **Similarité d'activité** (`etape_soeurs_session`, activité seulement,
   résidus comme l'étape 3b : distance et niveau de chaque cellule retirés) :
   STTC au repos ; même assemblée au K retenu (1/0).
3. **Test** (`test_soeurs_session`) : par session, moindres carrés sur toutes
   les paires de cellules colorées : résidu ~ 1 + pi + 6 indicateurs « même
   couleur c ». Les indicateurs absorbent tout effet de couleur indépendant
   de la distance (optique…) : la pente beta ne vient que de la comparaison
   des paires de même couleur proches (probables sœurs) aux lointaines (pas
   sœurs). beta = différence de similarité entre une paire sûrement sœur et
   une non sœur. p de session : couleurs mélangées dans le champ (500).
4. **Animaux** (`agregation_soeurs`) : moyenne des beta de leurs sessions ;
   **un test par mesure** : Wilcoxon unilatéral sur les animaux. Effet en
   % = beta / similarité brute moyenne à ≤ 100 µm. Dose-réponse (figure) :
   résidu moyen des paires de couleurs différentes, de même couleur avec
   pi = 0, puis par catégorie de pi.
5. **Par couleur** (7 nov., demande de l'utilisateur : le STTC peut dépendre
   de la couleur ; **descriptif**, 6 tests sans correction) : même pente
   limitée aux paires qui contiennent une cellule de la couleur c (résidu ~
   1 + pi + [même couleur c]) ; Wilcoxon sur les animaux ; lignes de couleur
   dans la dose-réponse (« diff. colour » = paires entre une cellule c et une
   autre couleur ; points avec < 20 paires non montrés) ;
   `soeurs_sessions_couleurs.csv`. Simulation (effet sœurs seulement pour B,
   décalage constant pour toutes les paires R-R) : pente B 0,037 pour 0,03,
   autres couleurs ~0 ; le décalage de R donne un niveau, pas une pente.

**Simulation** (clones de ~7 cellules étalés sur ~60 µm, 6 couleurs de
fréquences inégales, 6 sessions × 8 répétitions) :

| effet sœurs | artefact « même couleur » | beta (avec indicateurs) | beta sans indicateurs |
|---|---|---|---|
| 0 | 0 | −0,03 ± 0,04 | 0,00 |
| 0 | +0,1 partout | 0,00 ± 0,05 | **0,44 (faux effet)** |
| 0,3 | 0 | 0,37 ± 0,05 | 0,37 |
| 0,3 | +0,1 partout | 0,38 ± 0,07 | 0,84 |

Les indicateurs de couleur sont indispensables. L'effet chez les sœurs est
un peu surestimé (~+20 %) car pi est un peu sous-estimé (0,33 contre 0,37).

Sorties : `<session>\soeurs\soeurs.mat` ; `<PathSave>\soeurs_agregation\`
(`carte_soeurs.mat`, `angles_bandes.csv`, `angles_bandes_<animal>.png`,
`soeurs_sessions.csv`, `resume_soeurs.txt`, `soeurs.png` : A pi selon la
distance, B carte alignée sur les bandes, puis une dose-réponse par mesure).
Avec `config.seulement`, pi n'est estimé que sur ces sessions.

### Cœur clonal de chaque assemblée (fin de `main_soeurs`, 7 nov.)

Demandé pour pouvoir dire **oui / non pour chaque assemblée** (et relier
ensuite le cœur clonal aux séquences, par exemple), en plus du test général.

- **Score** (`coeurs_assemblees`) = somme des pi des paires de membres
  colorés = nombre attendu de paires de sœurs dans l'assemblée. Les paires de
  même couleur lointaines comptent 0 ; les couleurs rares pèsent plus ; un
  petit cœur de sœurs proches compte même dans une grosse assemblée.
- **Testable** : au moins 3 membres colorés (`min_colorees`) et au moins une
  paire de membres colorés à une distance où pi peut être > 0.
- **p_champ** : couleurs mélangées entre les cellules colorées du champ
  (mêmes membres, mêmes positions). **Descriptif seulement** : les clones
  étant compacts, une assemblée compacte contient forcément les sœurs de son
  voisinage, même si elle ne les choisit pas.
- **p_voisines (le verdict)** : 1000 fausses assemblées où chaque membre
  coloré est remplacé par une de ses 5 cellules colorées les plus proches
  (`k_voisins`), avec leurs vraies couleurs. Question : l'assemblée a-t-elle
  plus de sœurs que des cellules prises au même endroit ?
- **Verdict** (`agregation_coeurs`, FDR 10 % sur toutes les assemblées
  testables) : « cœur clonal » (p_voisines significatif), « contient des
  sœurs » (seulement p_champ), « pas de preuve », « pas testable ».
  « Pas de preuve » ne veut pas dire « pas clonale ». Pour relier à une autre
  propriété, utiliser plutôt le score continu `z_voisines`.

**Simulation** (6 champs, 36 assemblées de 20 cellules par cas ; part des
assemblées avec p < 0,05, avant FDR) :

| assemblée simulée | p_champ | p_voisines |
|---|---|---|
| compacte, aucun choix clonal | 0,50 | 0,08 |
| dispersée au hasard | 0,14 | 0,00 |
| 4 sœurs + 16 cellules au hasard | 0,50 | 0,28 |

p_champ n'est pas spécifique (d'où le verdict sur p_voisines, revu le 7 nov.
après cette simulation) ; p_voisines est à peu près juste mais ne voit
qu'environ un cœur sur quatre.

Sorties : `<session>\soeurs\coeurs.mat` ; dans `soeurs_agregation\` :
`coeurs_assemblees.csv` (une ligne par assemblée : score, p, verdict,
cellules et couleurs du cœur), `resume_coeurs.txt`, `coeurs_assemblees.png`
(A score de chaque assemblée par animal, B les deux p, C part des assemblées
avec un cœur clonal par animal). Demande les étapes 2 et 3 (K retenu).

## Interface : 3e image = une assemblée en couleurs Brainbow (8 nov.)

Demande de l'utilisateur (l'ancienne image des couleurs de toutes les cellules
ne servait à rien) : l'image de fond (GCaMP, ou l'image couleur choisie avec
« Colour background image... ») avec UNE assemblée, chaque membre rempli de sa
couleur Brainbow, rien d'entouré ; si p_minp < 0,05, le cadre de l'image prend
la couleur en excès ; boutons ▲ ▼ sur le côté pour passer d'une assemblée à
l'autre (le menu « Assembly » la choisit aussi). Titre sur une ligne : cellules,
cellules colorées, couleur en excès, p (non corrigé) et p_minp, pris dans
l'étape 3 si elle a utilisé ce K, sinon dans `main_clonal_tous_K` ; sinon
« not computed at this K ». Couleurs seulement après l'étape 3, comme avant.
Texte de la session : pour chaque assemblée, p (non corrigé) et p_minp, à
TOUS les K : ceux de l'étape 3 au K retenu, ceux de `main_clonal_tous_K` s'il a
tourné, sinon calculés par l'interface pour la session et le K affichés (même
test que l'étape 3, variante principale, sans les descripteurs lents ;
quelques secondes, gardés en mémoire). Seulement une fois l'étape 3 faite
pour la session (couleurs).

Image du milieu (8 nov., demande de l'utilisateur) : toutes les cellules du
champ, chacune remplie de sa couleur Brainbow (gris : sans couleur, y compris
les doublons et les ROI retirées pour leur taille), rien d'entouré, sur le même
fond que la 3e image ; titre : nombre de cellules colorées. Tant que les
couleurs sont cachées (étape 3 pas faite, ou case « Show brainbow colours »
décochée) : les cellules des assemblées comme avant.

Camemberts (9 nov., demande de l'utilisateur) : colonne à droite de la 3e
image (le raster des SCE est un peu plus étroit pour faire la place). En haut,
les couleurs des cellules colorées du champ (l'urne du test) ; en bas, celles de
l'assemblée de la 3e image (suit les boutons ▲ ▼). Pourcentage écrit dans les
parts d'au moins 10 % ; titre : nombre de cellules colorées. Cellules sans
couleur non comptées. Seulement quand les couleurs sont affichées.

## Biais clonal pour chaque K : `main_clonal_tous_K` (8 nov., exploratoire)

Demande de l'utilisateur : donner le p initial et le p_minp, pour chaque K.
Même test que l'étape 3 (variante principale) pour tous les K de
`config.liste_K` et le consensus (K = 0), sans rien modifier. Par assemblée :
`p_initial` (p hypergéométrique de la couleur la plus enrichie, NON corrigé
pour les 6 couleurs : sous H0, p_initial < 0,05 arrive jusqu'à ~6 x 5 % du
temps) et `p_minp` (corrigé par permutation), chacun avec son taux attendu par
hasard (`taux_nul_initial`, `taux_nul`, ajoutés au test de l'étape 3). Par K :
testables, observé contre attendu, FDR 10 %, Wilcoxon sur les animaux (sans
correction pour le nombre de K). Sorties dans `<PathSave>\clonal_tous_K\`.
**Ne pas choisir le K d'après ces résultats** (meilleur de ~20 essais) ; les
mêmes assemblées reviennent d'un K à l'autre.

Limites de p_minp (quand il peut être trop petit) : (1) les couleurs sont
supposées échangeables dans tout le champ ; or les clones regroupent les
couleurs, et une assemblée compacte contient alors trop d'une couleur sans
aucun choix clonal (simulation du 4 nov. : 36 % de faux positifs pour des
assemblées compactes, 12 % pour des assemblées dispersées ; contrôles
'blocs' et 'territoire') ; (2) toutes les cellules colorées de l'urne sont
supposées pouvoir être membres ; si l'activité ou la détection dépend de la
couleur, la couleur la plus active semble enrichie (contrôle 'activite') ;
(3) une cellule comptée deux fois (doublon de ROI) double sa couleur
(corrigé par le retrait des doublons à l'étape 1) ; (4) choisir le K, la
variante ou le seuil d'après le résultat. Les ex aequo sont comptés contre
nous (p un peu trop grand, pas trop petit).

## Co-participation aux SCE, mesure 'sce' des étapes 3b et 3d (8 nov.)

Question de l'utilisateur : « est-ce que les cellules de même couleur
participent plus aux mêmes SCE ? » (jamais testé directement : seulement la
participation cellule par cellule, l'arbre 'sce' de 3c et « même assemblée »).
`3b_paires/coparticipation_sce.m` : coefficient phi des deux lignes 0/1 de
Race (tient compte du nombre de SCE de chaque cellule), puis mêmes résidus que
le STTC (distance et niveau de chaque cellule retirés). Ajoutée à
`config.paires.mesures = {'sttc', 'dff', 'sce'}` (étape 3b : même couleur
contre couleurs différentes ; figure `paires_sce_<rayon>um.png`) et à
`config.soeurs.mesures = {'sttc', 'sce', 'coassemblee'}` (étape 3d : selon la
probabilité d'être sœurs, couleur par couleur). Phi a une moyenne proche de
0 : l'effet est donné en différence de phi, pas en %.
Premier regard (444175_221121_plane1 seulement, couleurs mélangées) : phi même
− différente = +0,0095 à 0-50 µm (p = 0,10), +0,0065 à 200-400 µm (p = 0,004),
+0,0026 sur toutes les paires (p = 0,04) ; B +0,011 à ≤ 200 µm, R −0,006.

## Doublons : même cellule comptée deux fois (7 nov.)

Trouvé sur 444175_221121_plane1 : ROI 4 et 224 (centres à 0,2 µm, traces au
repos corrélées à 0,94), 2 et 154 (0,93), 14 et 42 (0,80), 28 et 214 (0,75),
1 et 5 (0,70) ; une paire ordinaire est à 0,03. Les analyses de paires
(3b, 3d, cœurs) les excluaient déjà (masques qui se touchent), mais elles
comptaient deux fois dans les SCE et les assemblées : 29 des 406 SCE de
cette session n'atteignent le seuil que grâce à elles.

- `1_fonctionnel/detecter_doublons.m` : traces au repos corrélées à plus de
  `config.pre.doublons_r_min` (0,7) **et** masques qui se touchent ou centres
  à moins d'un diamètre ; la ROI au plus petit masque est retirée.
- `config.pre.retirer_doublons = true` (défaut ; `false` dans le réglage
  'ancien_rang') : à l'étape 1, raster de la ROI retirée mis à zéro avant
  les SCE (elle n'entre plus dans les SCE ni les assemblées) ;
  `diagnostic.doublons` dans `pretraitement.mat`. **Change config.pre : toute
  l'étape 1 est refaite.**
- `utilitaires/couleurs_sans_doublons.m` : aux étapes 3, 3b, 3c, 3d et 5, la
  ROI retirée est traitée comme non colorée (sa couleur compterait sinon deux
  fois). Sans effet tant que l'étape 1 n'a pas été refaite.
- **Tailles des ROI** (7 nov., demande de l'utilisateur) : `tailles_roi` ;
  `config.pre.retirer_tailles = true`, ROI de moins de `roi_min_px` (4) pixels
  ou de plus de `roi_max_facteur` (3) fois la surface médiane de la session
  retirées comme les doublons (avant eux). Règle relative pour les grandes,
  pour ne pas dépendre du zoom ; seuils à ajuster après `diagnostic_doublons`.
- `diagnostic_doublons` : sans rien recalculer, liste les doublons et les ROI
  hors taille de toutes les sessions déjà prétraitées (`doublons_paires.csv`,
  `roi_hors_taille.csv`, `doublons_sessions.csv` avec les tailles min /
  médiane / max dans `PathSave`), avec leurs couleurs et s'ils sont dans une
  même assemblée.

## Petits SCE et assemblées : `diagnostic_SCE_assemblees` (7 nov.)

Question de l'utilisateur : les petits SCE sont en nombre proche du hasard,
faut-il les éliminer ? (« si la participation moyenne est de 30 %, pour avoir
des assemblées de 30 il me faut des SCE de 10 »). Sans rien recalculer : SCE
du hasard (train de chaque cellule décalé dans le repos, même détection que
l'étape 1) ; pour chaque assemblée, membres présents dans ses activations
retenues (PCl_final) ; activations imitables par le hasard = SCE du hasard
avec au moins autant de membres que le 10e percentile des vraies
activations (un maximum). `diagnostic_SCE_assemblees` (K retenu) ou
`diagnostic_SCE_assemblees([], 10)` ; dans `<PathSave>\diagnostic_SCE\` :
`diagnostic_SCE_sessions.csv`, `diagnostic_SCE_assemblees.csv`,
`resume_diagnostic_SCE.txt`, `diagnostic_SCE.png` (figure_diagnostic_SCE ;
aussi sur d'anciens CSV : `figure_diagnostic_SCE(dossier_des_csv)`). Visibles
dans l'interface (« Results figures », « Results summary »), comme les figures
et résumés des étapes 3c, 3d (sœurs, cœurs clonaux) et de main_clonal_tous_K.

**444175_221121_plane1, K = 10** : seuil 9 cellules, 406 SCE (231 attendus
par hasard) ; 186 petits SCE (10-12 cellules), dont 81 sont des activations
d'une assemblée : elles contiennent 3,5 à 6 membres de leur assemblée
(médiane) contre 0 à 1 dans les SCE du hasard de même taille (95e percentile
2 à 3) ; participation des membres 27-50 % ; 241 activations, au plus 2 à 4
imitables par le hasard (1-2 %). Un seuil à 13 cellules (percentile 99)
perdrait environ un tiers des vraies activations. Les petits SCE qui ne sont
l'activation d'aucune assemblée (105) sont probablement le bruit ; ils vont
surtout dans les clusters rejetés. Le décalage des cellules surestime le
hasard quand il y a de vraies co-activations : « autant de petits SCE que le
hasard » ne veut pas dire qu'ils sont tous du hasard.

## Étape 5 : séquences pendant la course (`main_sequences`, 6 nov.)

Case « Step 5: sequences » de l'interface, ou `main_sequences`. D'abord : y
a-t-il des cellules activées en séquence pendant la course ? Ensuite : les
cellules de même couleur sont-elles plus proches dans la séquence ?

- Méthode de l'utilisateur (`sequence_2026_01_27`), sans changement
  (`5_sequences/decalage_cellules.m`) : F de `data_recording` ÷ polynôme de
  degré 3, Savitzky-Golay (3, 5), z-score ; courses = vitesse lissée
  (gaussienne 20 trames) > 2 cm/s ; corrélation croisée de chaque cellule avec
  les débuts de course, sur les courses seulement (± 100 trames) ; cellules de
  la séquence = maximum > seuil d'Otsu et décalage < 0 ; rangées selon le décalage.
- Malvache et al. 2016 : PCA sur les traces lissées (σ = 5 s), composante
  choisie à la main, sa dérivée corrélée à chaque cellule, seuil d'Otsu ;
  validation par mélange de l'activité entre cellules dans chaque course.
- **Vérification (idée déjà dans le code d'origine, `test_spearman`)** : les
  courses sont tirées au hasard en deux moitiés ; la 1re range les cellules, la
  2e recalcule leur décalage ; rho de Spearman (unilatéral) entre les deux. Les
  moitiés n'ont aucune course en commun. Toutes sessions : rho moyen par
  animal, Wilcoxon unilatéral sur les animaux (rho > 0).
- Au moins 5 courses par session (`config.seq.min_courses`, choix de
  l'utilisateur) : 3 pour ranger, 2 pour vérifier ; sinon la session est notée
  « trop peu de courses ».
- Simulation sans séquence (550 sessions) : le p de Spearman d'une session est
  un peu optimiste (≈ 7 % de p < 0,05) ; ça ne touche que le nombre de sessions
  à p < 0,05, pas le test sur les animaux (fait sur rho). Avec séquence : rho ≈ 0,9.
- **Couleurs (test de l'utilisateur, 6 nov.)** (`distance_sequence_couleur`) :
  ordre de TOUTES les courses ; position = rang / (cellules de la séquence - 1) ;
  distance d'une paire = écart de position (fraction de la séquence), comme la
  « distance intra-assemblée » du code d'origine. Paires de cellules colorées
  de la séquence, même couleur contre couleurs différentes, sans correction ;
  p de session : couleurs mélangées entre les cellules colorées de la séquence
  (2000) ; animaux : moyenne de leurs sessions, Wilcoxon unilatéral (même
  couleur plus proche). Au moins 5 paires de même couleur et 5 de couleurs
  différentes (`config.seq.min_paires_meme`). Contrôle (pas une correction) :
  rho de Spearman entre distance dans la séquence et distance dans l'espace,
  parce que les couleurs sont regroupées dans l'espace. Figure
  `sequence_colour.png` (A même / différentes par animal ; B chaque session ;
  C contrôle spatial).
- Par session : `<session>\sequences\sequences.mat`, `sequence_session.png`
  (A : moitié qui range ; B : autre moitié, même ordre ; C : décalages des deux
  moitiés). Toutes sessions : `sequences_agregation\` (`sequences_sessions.csv`,
  `resume_sequences.txt`, `sequences.png`). Réglages : `config.seq`.

## Ce que les données simulées ont montré

- **Calibration.** Sans effet clonal, 4,2 % des assemblées testables sortent à
  p < 0,05 pour 3,7 % attendus ; `p_champ` < 0,05 dans 5,3 % des champs.
- **Fausses assemblées.** Sans vraies assemblées, mais avec des cellules très
  inégalement actives, le nouveau modèle nul n'en trouve aucune ; l'ancien
  modèle nul (même avec les autres corrections) en trouve dans 6 % des analyses.
- **Choix de K.** Avec la normalisation `'session'`, la règle choisissait
  souvent K = 3, qui fusionne deux assemblées (2 animaux sur 3 dans le test
  MATLAB, 8 sur 10 sur une simulation plus large), parce qu'elle étire à
  pleine échelle des différences minuscules de spécificité et de pureté
  (0,974 contre 0,971). Avec `'brute'`, le bon K (4) est retenu 10 fois sur
  10. La configuration utilise donc `'brute'` (aucune couleur n'intervient
  dans ce choix) ; `'session'` reste disponible. Sur les vraies données,
  vérifie le diagnostic de l'étape 2 (« meilleur K sur la borne basse »).
- **Référence `'assemblees'` et assemblées très clonales.** Quand une partie
  des assemblées est fortement clonale, leurs couleurs gonflent l'urne de
  référence ; les AUTRES assemblées paraissent alors en excès dans les
  couleurs restantes. Faux positifs parmi les assemblées non clonales
  (simulation, p < 0,05) :

  | Scénario | `'assemblees'` | `'champ'` | hasard attendu |
  |---|---|---|---|
  | aucun clone | 4,1 % | 4,4 % | ~4 % |
  | comme le test (2 assemblées sur 4 à 9/15 cellules d'un clone) | **17,7 %** | 5,3 % | ~4 % |
  | 10 assemblées de 20, 5 avec 6 cellules d'un clone | 5,6 % | 4,4 % | ~4 % |
  | 10 assemblées de 20, 3 avec 4 cellules d'un clone | 4,4 % | 3,8 % | ~4 % |

  La conclusion « il existe un biais clonal » reste valide (sans aucun clone,
  pas d'excès), mais la PROPORTION d'assemblées biaisées et leur liste sont
  surestimées si l'effet est fort. `'champ'` est aussi un peu plus puissant
  (96 % contre 88 % dans le scénario du test). Son point faible : un biais de
  détection selon la couleur, que `p_detection` et `p_F0` contrôlent.
- **Assemblée à deux clones.** Une assemblée qui contient deux clones de
  couleurs différentes n'est pas détectée par `p_minp` (0,054), qui cherche UNE
  couleur en excès ; le test global du champ (`p_champ`), lui, la voit.

## Pas encore traité

- Anciens scripts non relus : `raster_rastermap` (étape 1) et
  `figures_biais_clonal` (étape 4). Ils sont lancés avec tes noms de variables
  habituels ; une erreur s'affiche en avertissement sans bloquer. Un script de
  l'étape 1 qui demande `colorcell` échouera : c'est voulu.
- `figures_biais_clonal` est lancé une fois par variante, et son titre ne dit
  pas laquelle : la variable `nom_variante` est disponible pour l'y ajouter.
  Il classe les assemblées à FDR 20 % alors que le reste de l'agrégation
  utilise `fdr_seuil_figure` (10 %) : 9 « biaisées » dans sa figure contre
  7 dans la console pour la même variante.
- La version Python n'a pas encore ces derniers changements (scripts annexes).
- `score_sequence` n'est calculé nulle part : si tu l'ajoutes dans
  `assemblees.mat` (une valeur par assemblée, sans les couleurs), l'agrégation
  l'utilisera.
- Temps de calcul de l'étape 1 : les 100 surrogates à 50 réplicats sont le
  poste principal. `config.asm.kmeans_surrogate` peut être réduit à 50.

## Paires de même couleur dans les assemblées (`main_paires_assemblees`, 10 nov.)

« Mesurer au lieu de compter ». L'étape 3 teste chaque assemblée (p_minp) et
l'étape 4 compte les significatives ; une assemblée seule a peu de cellules
colorées et son test manque de puissance. Ici, une seule mesure par session,
toutes les assemblées ensemble : nombre de paires de cellules de même couleur
dans les assemblées, chaque paire pondérée par 1 / (part de sa couleur parmi
les cellules colorées du champ), comparé au même nombre après mélange des
couleurs. Excès = observé / attendu - 1 ; moyenne par animal ; Wilcoxon
unilatéral sur les animaux.

- Test principal : assemblées valides (critères de l'étape 4 :
  `config.agreg.min_cellules_asm`, `min_SCE_asm`), mélange dans tout le champ.
- Contrôle spatial : mélange par blocs (comme `champ_blocs`).
- Sensibilité : toutes les assemblées.
- Mêmes choix que l'étape 3 : urne = cellules colorées du champ, couleurs de
  moins de `config.clonal.nmin_couleur` cellules ignorées, K retenu à l'étape 2.
- 2000 mélanges (`config.paires_asm.nPerm` pour changer ; facultatif).
- Il faut les étapes 1 et 2. Lancer : `main_paires_assemblees`.
- Sorties dans `<analyse>\paires_assemblees\` : CSV par session et par animal,
  `resume_paires_assemblees.txt`, `paires_assemblees.png` (visibles par
  « Results summary » et « Results figures » de l'interface).
- Limites : ne dit pas quelle assemblée est clonale ; trompé comme l'étape 3 si
  la couleur est liée à la qualité du signal.

## FDR sur les assemblées qui pouvaient gagner (étape 4, 10 nov.)

Sensibilité ajoutée à l'étape 4 (« C ter »), le classement principal ne change
pas. Benjamini-Hochberg sur les seules assemblées testables, valides, dont le
meilleur p possible (`p_min_atteignable`) est sous
`config.agreg.fdr_filtre_pmin` (0,001 ; facultatif dans le config). Le tri ne
dépend que de la taille de l'assemblée et des couleurs du champ, pas du
résultat. Colonnes `fdr_filtre` et `q_filtre` dans le tableau des assemblées,
une ligne « Sensibilité » dans `resume_resultats.txt`. Seulement l'étape 4 à
relancer.

## Vérification des fichiers de couleurs (`verifier_couleurs`, 10 nov.)

Ne modifie rien. Pour chaque session : taille de colorcell contre le nombre
de cellules iscell de Fall.mat et le Raster de l'étape 1 (signale un
colorcell calculé sur toutes les ROI), valeurs hors 1-8, accord avec les
autres fichiers de couleurs du dossier, alignement couleurs / ROI (les
voisines doivent partager leur couleur plus que le hasard, et plus du tout
quand on décale colorcell d'une position), sessions qui partagent le même
fichier ou les mêmes couleurs. Composition de chaque session comparée aux
autres sessions de son animal, et composition moyenne de chaque animal.
Lancer : `verifier_couleurs` (ou avec `config.seulement`). Sorties dans
`<dossier de l'étape 1>\verification_couleurs\`.

## Visionneuse de l'image composite (`visionneuse_composite`, 10 nov.)

Bouton « Composite viewer... » de l'interface (session sélectionnée), ou en
ligne de commande `visionneuse_composite('<dossier suite2p>')`. Ouvre une
fenêtre avec le fichier dont le nom contient « composite » (.tif) du dossier
des données (sinon son dossier parent ; « Other file... » pour en choisir un
autre). Image ImageJ à plusieurs canaux : chaque canal affiché ou non, avec
sa couleur (défaut : canal 1 rouge, 2 vert, 3 bleu, 4 gris = calcium,
l'ordre écrit par color_predictor_model) et son gain. Masques GCaMP (cellules iscell = 1) en blanc ou dans
la couleur de colorcell, numéros parmi les cellules iscell = 1 (= ligne de
colorcell et du Raster), « Highlight ROI » pour retrouver une cellule. Ne modifie rien.

## Comparer toutes les analyses déjà faites (`comparer_analyses`, 11 nov.)

`comparer_analyses('E:\...\nocues')` cherche tous les résultats de l'étape 4
sous ce dossier (tous les dossiers d'étape 1, toutes les analyses, toutes les
variantes) et fait un tableau : seuil des SCE, K calculés, règle du K, K
retenus, assemblées testables, % à p_minp < 0,05, attendu, excès, biaisées
après FDR, Wilcoxon, GLMM, score, proportion estimée. Ne recalcule rien.
Sorties : `comparaison_analyses.csv`, `resume_comparaison_analyses.txt`.
Attention : choisir le réglage qui donne le plus d'assemblées biaisées,
c'est choisir d'après le résultat ; le tableau sert à juger la stabilité.

## Version 2 des assemblées : cellules partagées (`RACEORTHO_V2`, 11 nov.)

Version 1 (RACEORTHO_FDR2, héritée de Malvache 2016) : une cellule
significative dans plusieurs clusters de SCE est mise dans le seul où elle
participe le plus. Sur 444175_221121_plane1, K = 10 : 22 cellules sur 126
concernées, 26 appartenances sur 152 retirées, une assemblée (9 membres -> 3)
supprimée. Version 2 : la cellule reste membre de tous ces clusters ; le reste
est identique (test des membres, taille minimale, test des activations ; dans
le modèle nul pondéré, une cellule tirée compte pour toutes ses assemblées).
Même session en version 2 : 8 assemblées (12, 20, 19, 27, 19, 24, 22, 9
cellules) au lieu de 7, 22 cellules dans 2 ou 3 assemblées.

- Lancer : `config.version_assemblees = 'v2';` puis étapes 2 à 4 (interface
  ou `main_brainbow`). Pas d'étape 1 à refaire, pas de k-means.
- La version est une propriété de l'analyse : son nom est `v2`, commence par
  `v2_` ou finit par `_v2` (ajouté tout seul par `dossiers_analyse`). Dans
  l'interface : nommer la nouvelle analyse avec `_v2` à la fin.
- Fichiers : `<étape 1>\<session>\k<K>\assemblees_v2.mat`, créés à la demande
  (`fichier_assemblees`, `creer_assemblees_v2`) à partir de `assemblees.mat`,
  refaits si celui-ci est plus récent. Champs en plus : `assemblee_principale`
  (sert seulement au tri du raster), `n_appartenances`, `version_assemblees`.
- Tous les lecteurs (étapes 2, 3, 3d, 4, consensus, interface,
  `main_clonal_tous_K`, `main_paires_assemblees`, `diagnostic_SCE_assemblees`)
  passent par `fichier_assemblees`. Le consensus compte une paire une seule
  fois par K ; ses assemblées restent disjointes (coupe d'un arbre).
- Silhouettes, clusters de SCE et leur test : inchangés (avant RACEORTHO).
- Raster : cellules triées par leur assemblée principale ; une cellule
  partagée prend la couleur de la dernière assemblée affichée.

Contrôle 6 de `verifier_couleurs` (11 nov., 23 h) : l'image Brainbow
(compositecalcium.tif, dans le dossier des données) est-elle posée sur les
masques suite2p ? `alignement_composite` cherche, de -80 à +80 pixels, le
déplacement de l'image qui met le plus de lumière dans les masques ; bien
posée = meilleur décalage proche de 0. Colonnes `decalage_image_x_px`,
`decalage_image_y_px`, `contraste_masques` (lumière dans les masques /
moyenne de l'image), `bord_vide_px`, et le même contrôle pour le canal
calcium. Simulation : image décalée de (30, -20) px retrouvée exactement.

## Assemblées dans une partie des cellules (`config.pre.population`, 15 nov.)

Idée : il y aurait deux types d'assemblées, celles des cellules activées en
séquence pendant la course et celles des cellules inhibées pendant la course
et très actives au repos ; ces dernières pourraient noyer les premières dans
la détection des SCE.

```matlab
config.pre.population = 'sequence';        % SCE et assemblées parmi les cellules de la séquence de course
config.pre.population = 'hors_sequence';   % ... parmi toutes les autres cellules
% ne pas mettre le champ : toutes les cellules, comme avant
```

- Cellules de la séquence : même calcul que l'étape 5 (corrélation croisée
  avec les débuts de course, seuil d'Otsu, décalage < 0, sur toutes les
  courses), fait dans l'étape 1 à partir de l'activité et de la vitesse
  seulement. Réglages : ceux de `options_sequences` ; pour les changer ici,
  `config.pre.population_seq.<champ>`.
- Le raster des cellules écartées est mis à zéro : le seuil des SCE, les SCE
  et les assemblées sont cherchés parmi les cellules gardées.
- Couleurs : les cellules écartées sont traitées comme non colorées dans
  toutes les étapes suivantes. Le test compare donc une assemblée à la
  composition en couleurs de SA population, et les couleurs sont mélangées
  parmi les cellules gardées. (Dans l'interface, elles apparaissent sans
  couleur pour cette analyse.)
- Changer `config.pre.population` refait l'étape 1 : prendre un nouveau
  `config.PathSave` par population.
- Une session avec moins de 5 courses s'arrête avec un message : la retirer
  (`inclure = 0`) pour cette analyse.
- Screening : un screening par population (`config.pre.population` dans la
  base, un `opts.dossier` chacun).
- Le test de couleur de l'étape 5 (séquences) est à lancer sur une analyse
  faite avec toutes les cellules.

## Screening des réglages (`screening_reglages`, 13 nov.)

But : chercher, sur quelques animaux (446178, 444119, 444175 par défaut), les
réglages qui donnent le plus d'assemblées biaisées en couleur ; puis tester
ces réglages, fixés, sur les autres animaux et sur les sessions avec indice.

```matlab
config = config_brainbow();           % la « base » : tes réglages actuels
screening_reglages(config, struct('etapes', {{'plan'}}))   % ne calcule rien : sessions, réglages, nombre d'étapes 1 à faire
screening_reglages(config)            % tout (long), reprend où il s'est arrêté
screening_reglages(config, struct('etapes', {{'agregation'}}))             % refaire seulement le classement
screening_reglages(config, struct('facteurs', {{'sce_seuil', 'threshold_peak'}}))   % seulement ces réglages du prétraitement
```

Trois étapes :

1. **Étape 1 de chaque réglage du prétraitement** (`config.pre`), un seul
   changé à la fois par rapport à la base : seuil des SCE (8 règles), seuil
   des transients, largeur des SCE, distances minimales, lissage, garde-fou,
   cellules hyperactives. 24 réglages + la base. C'est la partie longue :
   chaque réglage refait l'étape 1 des sessions du screening. Dossiers
   d'étape 1 ordinaires dans `screening\etape1\<réglage>\`. Si l'étape 1 de
   `config.PathSave` a exactement les réglages de base, elle est lue telle
   quelle (rien n'y est écrit).
2. **Variantes d'assemblées**, sans refaire le k-means (il est relu) :
   6 règles de clusters (rang 90 / 95 / 99, max 95 / 99, tous 95) x 4 règles
   de membres (FDR_Q 0.05 / 0.10 / 0.20, participations >= 3) x 2 versions
   (v1, v2) x choix du K (chaque K, 'exces', 9 consensus) x 9 critères de
   l'étape 4 (min cellules 6 / 10 / 15, min SCE 5 / 10 / 20). Avec les
   réglages de base, les assemblées sont exactement celles du pipeline
   (mêmes tirages). Aucune couleur n'est lue jusqu'ici.
3. **Agrégation** : seule partie qui lit les couleurs. Pour chaque
   combinaison : assemblées testables, biaisées (p_minp < 0,05, test
   champ_globale), attendues par hasard, et l'excès (par animal, puis
   moyenne), comme le test A1 de l'étape 4.

Pour ajouter ou retirer des réglages : `screening\grille_screening.m`.

### Le tri par le résultat

Choisir le meilleur de milliers de réglages donne un bon score même sans
aucun effet. Le screening le mesure : les couleurs sont mélangées dans chaque
champ (500 fois, le même mélange pour tous les réglages), le screening est
refait à chaque fois, et on garde le meilleur excès de chaque mélange. C'est
la max-statistique (comme le min-p, mais sur les réglages).

- `p_seul` : p d'une combinaison comme si elle était la seule essayée.
- `p_selection` : part des mélanges dont le meilleur réglage fait au moins
  aussi bien. C'est la p corrigée pour tous les réglages essayés.

Les réglages retenus doivent ensuite être testés sur des animaux qui n'ont pas
servi au screening. Dans ce test, les animaux du screening sont à part.

### Sorties (`screening\resultats\`)

- `resume_screening.txt` : à lire en premier. Ce que donne le hasard après
  le tri, la meilleure combinaison, les 10 premiers résultats, l'effet moyen
  de chaque réglage, les réglages conseillés, les lignes de config à copier.
- `figure_screening.png` : A, toutes les combinaisons contre le meilleur de
  chaque mélange ; B à F, effet moyen de chaque niveau de chaque réglage
  (point noir), hasard (barre grise), meilleure combinaison (triangle).
- `screening_top.csv` (50 premières), `screening_combinaisons.csv` (toutes),
  `screening_facteurs.csv`, `screening_pretraitements.csv`, `screening.mat`.

L'effet moyen d'un niveau (moyenne sur des centaines de combinaisons) est
plus fiable que la meilleure combinaison seule.

### Durée et place sur le disque

L'étape 1 est refaite pour chaque réglage du prétraitement et chaque session
des animaux du screening : 24 réglages x (nombre de sessions) étapes 1 de
session. Le mode `'plan'` donne ce nombre ; la durée d'une étape 1 de session
est dans les `log_analyse_*.txt` (« Step 1 - ... done »). Chaque étape 1 de
session prend aussi la place d'un dossier de session ordinaire. Pour
commencer par les réglages les plus importants :
`struct('facteurs', {{'sce_seuil', 'threshold_peak', 'synchronous_frames'}})`.
On peut arrêter et relancer : ce qui est fini n'est pas refait, et
`{'agregation'}` classe à tout moment les réglages déjà finis.

### Vérifications (13 nov., données de test)

- Avec les réglages de base, les assemblées du screening sont identiques à
  celles du pipeline : version 1, version 2 et consensus (54 comparaisons
  sur 54).
- p_minp de la table contre le mélange complet des couleurs : mêmes valeurs.
- Couleurs mélangées (aucun effet), 200 screenings : `p_selection` < 0,05
  dans 14 (7 % ; 5 % attendus, écart dans le bruit d'un test fait avec 150
  à 200 mélanges au lieu de 500) ; excès moyen des combinaisons : +0,2 point.
- Biais implanté dans 2 assemblées de 2 animaux sur 3 : retrouvé
  (`p_selection` = 0,02).

### Limites

- Test champ_globale seulement ; les contrôles blocs et territoire sont
  faits par le pipeline normal sur le réglage retenu.
- p_minp calculé par une table (loi de la statistique selon le nombre de
  cellules colorées, 20 000 tirages par champ) : même test que l'étape 3,
  petites différences de tirage.
- Version 2 et consensus : mêmes assemblées que le pipeline à de rares
  différences de tirage près (activations).
- Prétraitement : un réglage à la fois. Le résumé propose la ligne pour
  essayer ensemble les meilleurs niveaux (`opts.pre_supplementaires`).

