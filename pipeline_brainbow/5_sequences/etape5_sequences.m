function statut = etape5_sequences(session, config)
%% =========================================================================
%  ÉTAPE 5 : y a-t-il des cellules activées en séquence pendant la course ?
%  -------------------------------------------------------------------------
%  Première partie, SANS couleur (6 nov.) : on vérifie d'abord que les
%  séquences existent et sont reproductibles ; la structure de couleur
%  viendra ensuite.
%  Méthode de l'utilisateur (sequence_2026_01_27), sans changement :
%    - F de data_recording, divisé par un polynôme de degré 3, filtre de
%      Savitzky-Golay (3, 5), z-score par cellule ;
%    - courses : vitesse lissée (gaussienne, 20 trames) > 2 cm/s ;
%    - pour chaque cellule : corrélation croisée avec les débuts de course,
%      sur les courses seulement, maximum et son décalage (decalage_cellules) ;
%    - cellules de la séquence : maximum > seuil d'Otsu et décalage < 0 ;
%      rangées selon leur décalage.
%  Vérification (idée déjà dans le code d'origine, la ligne test_spearman) :
%    les courses sont tirées au hasard en deux moitiés. La 1re moitié range
%    les cellules ; sur la 2e moitié, on recalcule leur décalage et on
%    compare (corrélation de Spearman entre les deux décalages des cellules
%    de la séquence). Une vraie séquence garde le même ordre : rho > 0.
%    Les deux moitiés sont indépendantes (aucune course en commun).
%  Aussi sur TOUTES les courses (decalage_tout, sequence_tout, ordre_tout).
%  Couleurs (6 nov., test de l'utilisateur) : sur l'ordre de toutes les
%  courses, distance dans la séquence des paires de même couleur contre
%  couleurs différentes (distance_sequence_couleur).
%  19 nov. : opts.validation = 'kfold' (validation_kfold) : rho et p_rho
%  viennent de 100 coupures (5-fold x 20 ; 50/50 x 50 sous 10 courses), p
%  par permutation ; test de couleur sur les cellules validées (>= 50 % des
%  coupures), rangées par leur décalage moyen sur les courses de test.
%  La coupure 50/50 unique reste pour la figure (rho_moitie, p_rho_moitie).
%  20 nov. (choix de l'utilisateur) : opts.selection = 'faux_debuts'
%  (choisir_cellules_sequence -> decalage_cellules_test) : même corrélation
%  croisée, mais sur le signal continu (courses non collées), décalage
%  jusqu'à opts.max_lag_s (40 s), et seuil par cellule contre des faux
%  débuts tirés pendant la course (opts.n_faux_debuts = 500), FDR 5 %, à la
%  place d'Otsu. 'otsu' = méthode d'avant.
%  20 nov. : opts.cellules_tests = 'robustes' : tests de couleur et
%  d'assemblées sur les cellules validées ET de position stable
%  (validation_kfold, ordre_robuste) ; 'validees' = avant.
%  20 nov. : opts.signal = 'cascade' : même analyse sur la probabilité de
%  spike de CASCADE au lieu du calcium ; résultats dans
%  <session>\sequences_cascade (dossier_sequences).
%  Assemblées (sequence_assemblees) : assemblées du K retenu à l'étape 2,
%  enrichies en cellules de la séquence ? regroupées dans la séquence ?
%  (mêmes cellules et même ordre que le test de couleur).
%  Écrit <session>\sequences\sequences.mat et sequence_session.png.
%  =========================================================================

opts = options_sequences(config);
dossier = fullfile(session.dossier, dossier_sequences(opts));   % 'sequences' (calcium) ou 'sequences_cascade'
if ~isfolder(dossier), mkdir(dossier); end
fichier = fullfile(dossier, 'sequences.mat');
if isfile(fichier), delete(fichier); end

%% 1. Activité et vitesse (data_recording, comme sequence_2026_01_27)
S = load(session.fichier_data, 'data_recording');
F = double(S.data_recording.F);
speed = double(S.data_recording.speed(:)');
fs = S.data_recording.sampling_rate;
clear S
[NCell, Nz] = size(F);
if strcmp(opts.signal, 'cascade')
    % 20 nov. (choix de l'utilisateur) : probabilité de spike de CASCADE
    % (spike_prob, PA estimés par image, déjà lissés sur 150 ms), calculée
    % à l'étape 1 (cascade_<modèle>.mat dans le dossier suite2p). Pas de
    % polynôme ni de Savitzky-Golay ; NaN des bords (fenêtre CASCADE) mis à
    % 0 ; z-score par cellule comme pour le calcium.
    fc = fullfile(session.dossier_suite2p, ['cascade_' opts.cascade_modele '.mat']);
    if ~isfile(fc)
        error('CASCADE pas calculé pour cette session (%s absent)', fc);
    end
    C = load(fc, 'spike_prob');
    DF = double(C.spike_prob);
    clear C
    if ~isequal(size(DF), [NCell, Nz])
        error('CASCADE : taille %dx%d, data_recording %dx%d', size(DF, 1), size(DF, 2), NCell, Nz);
    end
    DF(~isfinite(DF)) = 0;
else
    DF = F;
    ws = warning('off', 'all');
    for i = 1:NCell
        p0 = polyfit(1:Nz, DF(i, :), 3);
        DF(i, :) = DF(i, :) ./ polyval(p0, 1:Nz);
    end
    warning(ws);
    DF = sgolayfilt(DF', 3, 5)';
end
DF = normalize(DF, 2, 'zscore');
DF(~isfinite(DF)) = 0;      % cellule sans aucune activité (écart-type nul)
clear F

% vitesse : même longueur que la fluorescence (comme preprocessing_7)
if numel(speed) >= Nz
    speed = speed(1:Nz);
else
    speed(end+1:Nz) = NaN;
end

%% 2. Courses
vitesse_lissee = smoothdata(speed, 'gaussian', opts.lissage_vitesse);
en_course = double(vitesse_lissee > opts.seuil_vitesse);
changements = diff([0, en_course, 0]);
% comme le code d'origine : début = trame juste avant la course, fin =
% dernière trame de la course ; une course qui commence à la 1re trame est
% retirée
debuts = find(changements == 1) - 1;
fins = find(changements == -1) - 1;
garder = debuts >= 1;
debuts = debuts(garder);
fins = fins(garder);
n_courses = numel(debuts);
trames_course = find(en_course);       % candidats pour les faux débuts
if strcmp(opts.selection, 'faux_debuts')
    opts.max_lag_trames = round(opts.max_lag_s * fs);
else
    opts.max_lag_trames = opts.max_lag;
end

R = struct();
R.nom = session.nom;
R.animal = session.id_animal;
R.session = session.id_session;
R.opts = opts;
R.sampling_rate = fs;
R.NCell = NCell;
R.n_courses = n_courses;
R.debuts = debuts;
R.fins = fins;
R.duree_course_s = sum(fins - debuts + 1) / fs;

if n_courses < opts.min_courses
    statut = sprintf('%d course(s) : trop peu (minimum %d), session non analysée', n_courses, opts.min_courses);
    R.statut = statut;
    save(fichier, '-struct', 'R');
    return
end

%% 3. Moitié des courses pour ranger, l'autre pour vérifier
fixer_flux_aleatoire(session.graine, 5, 1);
melange = randperm(n_courses);
n_app = max(1, min(n_courses - 1, round(n_courses * opts.fraction_apprentissage)));
courses_app = sort(melange(1:n_app));
courses_test = sort(melange(n_app + 1:end));

% (seuil_* = seuil d'Otsu en mode 'otsu', p de chaque cellule en mode 'faux_debuts')
fixer_flux_aleatoire(session.graine, 10, 1);
[R.corr_app, R.decalage_app, R.seuil_app, R.sequence_app, R.ordre_app] = ...
    choisir_cellules_sequence(DF, debuts, fins, courses_app, opts, trames_course);
[R.corr_test, R.decalage_test, R.seuil_test, R.sequence_test] = ...
    choisir_cellules_sequence(DF, debuts, fins, courses_test, opts, trames_course);
[R.corr_tout, R.decalage_tout, R.seuil_tout, R.sequence_tout, R.ordre_tout] = ...
    choisir_cellules_sequence(DF, debuts, fins, 1:n_courses, opts, trames_course);
R.courses_app = courses_app;
R.courses_test = courses_test;

%% 4. Reproductibilité : même ordre sur des courses indépendantes ?
seq = R.sequence_app;
R.n_sequence = sum(seq);
R.n_sequence_tout = sum(R.sequence_tout);
R.n_deux = sum(seq & R.sequence_test);                       % retrouvées sur l'autre moitié
R.n_deux_hasard = sum(seq) * sum(R.sequence_test) / NCell;   % si les deux sélections étaient indépendantes
R.rho = NaN; R.p_rho = NaN;
if sum(seq) >= 5
    [R.rho, R.p_rho] = corr(R.decalage_app(seq)', R.decalage_test(seq)', 'Type', 'Spearman', ...
        'Tail', 'right');
end

% 19 nov. (demandé par l'utilisateur) : validation croisée répétée. La
% coupure 50/50 ci-dessus reste pour la figure de la session (rho_moitie).
R.validation = 'moitie (1 coupure 50/50)';
R.rho_moitie = R.rho; R.p_rho_moitie = R.p_rho;
ordre_couleur = R.ordre_tout;
if strcmp(opts.validation, 'kfold')
    V = validation_kfold(DF, debuts, fins, opts, session.graine, trames_course);
    R.kfold = V;
    R.validation = V.validation;
    R.rho = V.rho; R.p_rho = V.p_rho;
    R.n_sequence = V.n_sequence;            % moyenne sur les coupures (courses d'apprentissage)
    R.n_deux = V.n_deux; R.n_deux_hasard = V.n_deux_hasard;
    R.n_validees = V.n_validees;
    R.n_robustes = V.n_robustes;
    ordre_couleur = V.ordre_valide;         % test de couleur : cellules validées, ordre sur les courses de test
    if strcmp(opts.cellules_tests, 'robustes')
        ordre_couleur = V.ordre_robuste;    % 20 nov. : validées ET position stable
    end
end
R.cellules_tests = opts.cellules_tests;

%% 5. Activité moyenne autour du début de course (figure)
psth = @(courses) moyenne_autour(DF, debuts(courses), opts.psth_avant, opts.psth_apres);
[R.psth_app, R.n_psth_app] = psth(courses_app);
[R.psth_test, R.n_psth_test] = psth(courses_test);
R.psth_app = single(R.psth_app);
R.psth_test = single(R.psth_test);

%% 6. Couleurs : même couleur plus proche dans la séquence ?
path = avec_separateur(session.dossier_suite2p);
[outline_gcampx, ~, gcamp_centroid] = charger_geometrie(path);
colorcell = charger_couleurs(path, session.fichier_colorcell);
colorcell = couleurs_sans_doublons(colorcell, session);    % doublons de l'étape 1 : non colorés
verifier_alignement(NCell, colorcell, outline_gcampx, session.nom);
um = session.um_par_px;
if isnan(um), um = config.um_par_px_defaut; end
R.colorcell = colorcell(:);
R.ordre_couleur = ordre_couleur;
Cc = distance_sequence_couleur(ordre_couleur, colorcell(:), gcamp_centroid * um, opts.nPerm_couleur, ...
    opts.min_paires_meme, session.graine);
noms_c = fieldnames(Cc);
for k = 1:numel(noms_c)
    R.(['coul_' noms_c{k}]) = Cc.(noms_c{k});
end

%% 7. Assemblées (étape 2) : suivent-elles la séquence ? (20 nov.)
R.asm_statut = 'pas d''assemblées';
try
    choix = load(fullfile(session.dossier, 'choix_K.mat'));
    if strcmp(choix.statut_choix, 'ok')
        fichier_asm = fichier_assemblees(session.dossier, session.dossier_source, choix.best_K);
        As = load(fichier_asm, 'assemblyortho');
        R.asm_K = choix.best_K;
        R.assemblees = As.assemblyortho;
        Ca = sequence_assemblees(ordre_couleur, NCell, As.assemblyortho, opts.nPerm_couleur, ...
            opts.min_paires_meme, session.graine);
        noms_a = fieldnames(Ca);
        for k = 1:numel(noms_a)
            R.(['asm_' noms_a{k}]) = Ca.(noms_a{k});
        end
        R.asm_statut = 'ok';
    else
        R.asm_statut = choix.statut_choix;
    end
catch err
    R.asm_statut = ['erreur : ' err.message];
end

statut = sprintf(['%d courses, %.0f/%d cellules dans la séquence (%.0f %%), ordre sur les courses de test ' ...
    '(%s) : rho = %.2f, p = %.3g ; même couleur : %+.1f %% de distance dans la séquence, p = %.3g'], ...
    n_courses, R.n_sequence, NCell, 100 * R.n_sequence / NCell, R.validation, R.rho, R.p_rho, Cc.pct, Cc.p);
if strcmp(R.asm_statut, 'ok')
    statut = sprintf('%s ; assemblées : %d/%d enrichies, même assemblée %+.1f %%, p = %.3g', statut, ...
        R.asm_n_enrichies, R.asm_n_assemblees, R.asm_pct, R.asm_p);
end
if isfield(R, 'n_validees')
    statut = sprintf('%s ; %d cellules validées, %d robustes', statut, R.n_validees, R.n_robustes);
end
R.statut = statut;
save(fichier, '-struct', 'R');

if ~isfield(config, 'faire_figures') || config.faire_figures
    try
        figure_sequence_session(dossier);
    catch err
        warning('figure_sequence_session : %s', err.message);
    end
end
end


function [M, n] = moyenne_autour(DF, debuts, avant, apres)
% Activité moyenne des cellules de debut-avant à debut+apres (courses
% assez loin des bords seulement)
Nz = size(DF, 2);
debuts = debuts(debuts - avant >= 1 & debuts + apres <= Nz);
n = numel(debuts);
M = nan(size(DF, 1), avant + apres + 1);
if n == 0, return; end
M = zeros(size(DF, 1), avant + apres + 1);
for k = 1:n
    M = M + DF(:, debuts(k) - avant:debuts(k) + apres);
end
M = M / n;
end
