function A = agregation_paires(sessions, config)
%AGREGATION_PAIRES  Regroupe l'étape 3b de toutes les sessions. L'animal est
%   l'unité statistique.
%   Pour chaque mesure x permutation x couleur x tranche de distance :
%     - par animal : delta moyen de ses sessions, z combiné (Stouffer) des
%       p de permutation de ses sessions
%     - entre animaux : delta moyen +/- sem, Wilcoxon (delta > 0, animaux),
%       z combiné des animaux (Stouffer) et son p unilatéral, q de
%       Benjamini-Hochberg sur les couleurs x tranches de chaque mesure
%   Écrit dans <PathSave>\paires_agregation\ : paires_toutes_sessions.csv,
%   paires_par_animal.csv, paires_resultats.csv, paires_agregation.mat,
%   resume_paires.txt (le test principal en une phrase, 5 nov.) et
%   paires_agregation.png (figure_paires_resultats).
%   Test principal (5 nov.) : mesure principale (STTC), permutation
%   principale ('couplage'), toutes couleurs, paires à moins de
%   config.paires.rayon_principal_um (disques) ; Wilcoxon unilatéral sur les
%   animaux (delta > 0) et p de permutation des sessions combinés (Stouffer).

opts = config.paires;
dossier = fullfile(config.PathSave, 'paires_agregation');
if ~isfolder(dossier), mkdir(dossier); end

lignes = {};
diag_taux = cell(0, 2);
diag_couplage = cell(0, 2);
couplage_couleurs = zeros(0, 6);
for s = 1:height(sessions)
    f = fullfile(sessions.dossier{s}, 'paires', 'paires.mat');
    if ~isfile(f), continue; end
    P = load(f);
    lignes{end+1, 1} = P.res; %#ok<AGROW>
    diag_taux(end+1, :) = {sessions.nom{s}, P.p_taux_couleur}; %#ok<AGROW>
    if isfield(P, 'p_couplage_couleur')
        diag_couplage(end+1, :) = {sessions.nom{s}, P.p_couplage_couleur}; %#ok<AGROW>
        couplage_couleurs(end+1, :) = P.couplage_par_couleur; %#ok<AGROW>
    end
end
A = table();
if isempty(lignes)
    fprintf('Agrégation des paires : aucune session avec des résultats.\n');
    return
end
R = vertcat(lignes{:});
writetable(R, fullfile(dossier, 'paires_toutes_sessions.csv'));
n_sessions = numel(lignes);
n_animaux = numel(unique(R.animal));
% Anneaux disjoints (sessions refaites depuis le 5 nov.) : seulement pour le
% modèle mixte sur toutes les distances ; le reste utilise disques ou tranches
Rann = R([], :);
if ismember('forme', R.Properties.VariableNames)
    est_anneau = strcmp(R.forme, 'anneau');
    Rann = R(est_anneau, :);
    R = R(~est_anneau, :);
end
% delta = moyenne des paires de même couleur - moyenne des paires de couleurs
% différentes, SANS retirer ce que donnent les permutations (choix de
% l'utilisateur, 5 nov.) : les tests sur les animaux et le modèle mixte
% comparent directement même couleur et couleurs différentes. Les
% permutations ne servent qu'aux p de chaque session (et à leur combinaison).
% (delta_hasard reste enregistré dans paires_toutes_sessions.csv.)

%% Par animal
Rt = R(R.testable, :);
[G, mes, perm, coul, dmin, dmax, anim] = findgroups(Rt.mesure, Rt.permutation, Rt.couleur, ...
    Rt.dist_min_um, Rt.dist_max_um, Rt.animal);
Ani = table(mes, perm, coul, dmin, dmax, anim, ...
    'VariableNames', {'mesure', 'permutation', 'couleur', 'dist_min_um', 'dist_max_um', 'animal'});
Ani.n_sessions     = splitapply(@numel, Rt.delta, G);
Ani.delta          = splitapply(@mean, Rt.delta, G);
% z de chaque session tiré de son p de PERMUTATION (et non de (delta - moyenne) / sd :
% le hasard est asymétrique quand des assemblées existent, et ce z-là exagère)
p_borne = min(max(Rt.p, 1 / (opts.nPerm + 1)), 1 - 1 / (opts.nPerm + 1));
Rt.z_p = norminv(1 - p_borne);
Ani.z              = splitapply(@(z) sum(z) / sqrt(numel(z)), Rt.z_p, G);
Ani.n_paires_meme  = splitapply(@sum, Rt.n_paires_meme, G);
Ani.enrichissement = splitapply(@(e) mean(e, 'omitnan'), Rt.enrichissement, G);
% moyennes des paires de même couleur / de couleurs différentes et similarité
% brute moyenne (sessions de l'étape 3b depuis le 5 nov.)
for v = {'moyenne_meme', 'moyenne_diff', 'similarite_moyenne'}
    if ismember(v{1}, Rt.Properties.VariableNames)
        Ani.(v{1}) = splitapply(@(x) mean(x, 'omitnan'), Rt.(v{1}), G);
    else
        Ani.(v{1}) = nan(height(Ani), 1);
    end
end
writetable(Ani, fullfile(dossier, 'paires_par_animal.csv'));

%% Entre animaux
[G2, mes, perm, coul, dmin, dmax] = findgroups(Ani.mesure, Ani.permutation, Ani.couleur, ...
    Ani.dist_min_um, Ani.dist_max_um);
A = table(mes, perm, coul, dmin, dmax, ...
    'VariableNames', {'mesure', 'permutation', 'couleur', 'dist_min_um', 'dist_max_um'});
A.n_animaux      = splitapply(@numel, Ani.delta, G2);
A.delta_moyen    = splitapply(@mean, Ani.delta, G2);
A.sem            = splitapply(@(x) std(x) / sqrt(numel(x)), Ani.delta, G2);
A.p_wilcoxon     = splitapply(@(x) signrank(x, 0, 'tail', 'right'), Ani.delta, G2);
A.z_combine      = splitapply(@(z) sum(z) / sqrt(numel(z)), Ani.z, G2);
A.p_combine      = 1 - normcdf(A.z_combine);
A.n_paires_meme  = splitapply(@sum, Ani.n_paires_meme, G2);
A.enrichissement = splitapply(@mean, Ani.enrichissement, G2);
A.moyenne_meme   = splitapply(@(x) mean(x, 'omitnan'), Ani.moyenne_meme, G2);
A.moyenne_diff   = splitapply(@(x) mean(x, 'omitnan'), Ani.moyenne_diff, G2);
A.similarite_moyenne = splitapply(@(x) mean(x, 'omitnan'), Ani.similarite_moyenne, G2);
A.delta_pct      = 100 * A.delta_moyen ./ A.similarite_moyenne;   % delta en % de la similarité moyenne
A.delta_pct(strcmp(A.mesure, 'sce')) = NaN;   % phi de co-participation : moyenne proche de 0, pas de %
A.q = nan(height(A), 1);
[Gm, ~, ~] = findgroups(A.mesure, A.permutation);
for g = 1:max(Gm)
    i = find(Gm == g);
    A.q(i) = bh_q(A.p_combine(i));
end
% Enrichissement spatial : toujours par rapport à une répartition des couleurs
% au hasard dans le champ (permutation globale), quelle que soit la ligne
% (5 nov.) la clé contient les DEUX bornes : avec les disques, dist_min vaut 0
% partout, et tous les rayons recevaient l'enrichissement du premier (50 µm)
cle = strcat(A.mesure, '|', A.couleur, '|', string(A.dist_min_um), '|', string(A.dist_max_um));
glob = strcmp(A.permutation, 'globale');
if any(glob)
    [trouve, j] = ismember(cle, cle(glob));
    enr_glob = A.enrichissement(glob);
    A.enrichissement(trouve) = enr_glob(j(trouve));
    A.enrichissement(~trouve) = NaN;
end
A = sortrows(A, {'mesure', 'permutation', 'couleur', 'dist_min_um'});
writetable(A, fullfile(dossier, 'paires_resultats.csv'));

%% Affichage
fprintf('\n=========== PAIRES : même couleur contre couleurs différentes ===========\n');
fprintf(['%d sessions, %d animaux. Unité : l''animal. p_combine : p de permutation des sessions ' ...
    'combinés (Stouffer, unilatéral).\n'], ...
    n_sessions, n_animaux);
p_taux = cell2mat(diag_taux(:, 2));
fprintf('Activité au repos différente selon la couleur (Kruskal-Wallis p < 0.05) : %d / %d sessions\n', ...
    sum(p_taux < 0.05), sum(~isnan(p_taux)));
if ~isempty(diag_couplage)
    p_coup = cell2mat(diag_couplage(:, 2));
    fprintf('Couplage au réseau différent selon la couleur (Kruskal-Wallis p < 0.05) : %d / %d sessions\n', ...
        sum(p_coup < 0.05), sum(~isnan(p_coup)));
    fprintf('Couplage moyen par couleur (z dans chaque session, moyenne des sessions) :\n');
    disp(array2table(mean(couplage_couleurs, 1, 'omitnan'), ...
        'VariableNames', {'R', 'G', 'B', 'RG', 'RB', 'GB'}))
    fprintf(['(> 0 : cellules de cette couleur plus synchrones avec l''ensemble du réseau. ' ...
        'Si l''excès des paires de même couleur disparaît avec la permutation ''couplage'', ' ...
        'il vient de ce couplage individuel, pas d''un lien propre aux paires.)\n']);
end
mesures = liste_mesures_paires(opts);
for im = 1:numel(mesures)
    for ip = 1:numel(opts.permutations)
        sel = strcmp(A.mesure, mesures{im}) & strcmp(A.permutation, opts.permutations{ip});
        if im == 1 && ip == 1, etiquette = 'PRINCIPAL'; else, etiquette = 'sensibilité'; end
        fprintf('\n--- %s, permutation %s (%s) ---\n', mesures{im}, opts.permutations{ip}, etiquette);
        disp(A(sel & strcmp(A.couleur, 'toutes'), {'dist_min_um', 'dist_max_um', 'n_animaux', ...
            'delta_moyen', 'sem', 'p_wilcoxon', 'p_combine', 'q', 'enrichissement'}))
        if im == 1 && ip == 1
            fprintf('Par couleur :\n');
            disp(A(sel & ~strcmp(A.couleur, 'toutes'), {'couleur', 'dist_min_um', 'dist_max_um', ...
                'n_animaux', 'delta_moyen', 'p_wilcoxon', 'p_combine', 'q', 'enrichissement'}))
        end
    end
end
fprintf(['enrichissement : paires de même couleur observées / attendues à cette distance ' ...
    '(> 1 : couleurs regroupées dans l''espace, ce qu''on attend des clones).\n']);

%% Test principal, résumé, sauvegarde et figure (5 nov.)
[bornes, cumulatif, etiquettes] = distances_paires(opts);
if cumulatif
    rayons = bornes(2:end);
    R0 = opts.rayon_principal_um;
    if ~any(rayons == R0)
        [~, i0] = min(abs(rayons - R0));
        R0 = rayons(i0);
    end
    principal_dist = [0, R0];
else
    principal_dist = bornes(1:2);       % tranches : les paires les plus proches
end
if ~exist('p_coup', 'var'), p_coup = []; end
MD = modele_distances(Rann, opts);
[AniB, AB] = paires_brutes(sessions, config, opts);
[MR, AniR, AR] = regroupement_anneaux(Rann, opts);
lignes_resume = resume_paires(A, opts, principal_dist, cumulatif, n_sessions, n_animaux, p_taux, p_coup, MD);
lignes_resume = [lignes_resume, ligne_brute(AB, opts, principal_dist), lignes_regroupement(MR)];
fid = fopen(fullfile(dossier, 'resume_paires.txt'), 'w', 'n', 'UTF-8');
if fid >= 0
    fprintf(fid, '%s\n', lignes_resume{:});
    fclose(fid);
end
fprintf('\n%s\n', lignes_resume{:});
date_paires = datestr(now, 'yyyy-mm-dd HH:MM'); %#ok<NASGU>
Rsess = R;    % une ligne par session (figure : un point par session)
save(fullfile(dossier, 'paires_agregation.mat'), 'A', 'Ani', 'Rsess', 'Rann', 'MD', 'AniB', 'AB', 'MR', 'AniR', 'AR', ...
    'opts', 'principal_dist', ...
    'cumulatif', 'etiquettes', 'bornes', 'n_sessions', 'n_animaux', 'lignes_resume', 'date_paires');
try
    figure_paires_resultats(dossier);
catch err
    warning('figure_paires_resultats non exécutée : %s', err.message);
end
if any(strcmp(opts.mesures, 'sce'))
    try
        figure_paires_resultats(dossier, [], 'sce');     % co-participation aux SCE (8 nov.) : paires_sce_<rayon>um.png
    catch err
        warning('figure_paires_resultats (sce) non exécutée : %s', err.message);
    end
end
end


function [AniB, AB] = paires_brutes(sessions, config, opts)
% (5 nov., demande de l'utilisateur) Similarité BRUTE, sans résidus, des paires
% de même couleur et de couleurs différentes (paires_brutes_session : gardé
% dans <session>\paires\paires_brutes.mat, recalculé si les distances ont
% changé). Par animal : moyenne de ses sessions ; entre animaux : Wilcoxon
% unilatéral (même > différente).
AniB = table(); AB = table();
[bornes, ~] = distances_paires(opts);
dist_cfg = 'tranches';
if isfield(opts, 'distances'), dist_cfg = opts.distances; end
lb = {};
for s = 1:height(sessions)
    if ~isfile(fullfile(sessions.dossier{s}, 'paires', 'paires.mat')), continue; end
    fb = fullfile(sessions.dossier{s}, 'paires', 'paires_brutes.mat');
    try
        Tb = [];
        if isfile(fb)
            B = load(fb);
            if isfield(B, 'bornes') && isequal(B.bornes, bornes) && strcmp(B.distances, dist_cfg)
                Tb = B.T;
            end
        end
        if isempty(Tb)
            Tb = paires_brutes_session(table2struct(sessions(s, :)), config);
        end
        lb{end+1, 1} = Tb; %#ok<AGROW>
    catch err
        warning('Similarité brute non calculée pour %s : %s', sessions.nom{s}, err.message);
    end
end
if isempty(lb), return; end
Rb = vertcat(lb{:});
Rb.delta = Rb.brut_meme - Rb.brut_diff;
Rb = Rb(isfinite(Rb.delta), :);
if isempty(Rb), return; end
[G, mes, dmin, dmax, anim] = findgroups(Rb.mesure, Rb.dist_min_um, Rb.dist_max_um, Rb.animal);
AniB = table(mes, dmin, dmax, anim, 'VariableNames', {'mesure', 'dist_min_um', 'dist_max_um', 'animal'});
AniB.brut_meme = splitapply(@mean, Rb.brut_meme, G);
AniB.brut_diff = splitapply(@mean, Rb.brut_diff, G);
AniB.delta     = splitapply(@mean, Rb.delta, G);
[G2, mes, dmin, dmax] = findgroups(AniB.mesure, AniB.dist_min_um, AniB.dist_max_um);
AB = table(mes, dmin, dmax, 'VariableNames', {'mesure', 'dist_min_um', 'dist_max_um'});
AB.n_animaux   = splitapply(@numel, AniB.delta, G2);
AB.brut_meme   = splitapply(@mean, AniB.brut_meme, G2);
AB.brut_diff   = splitapply(@mean, AniB.brut_diff, G2);
AB.delta_moyen = splitapply(@mean, AniB.delta, G2);
AB.sem         = splitapply(@(x) std(x) / sqrt(numel(x)), AniB.delta, G2);
AB.delta_pct   = 100 * AB.delta_moyen ./ AB.brut_diff;     % % au-dessus des paires de couleurs différentes
AB.delta_pct(strcmp(AB.mesure, 'sce')) = NaN;
AB.n_pos       = splitapply(@(x) sum(x > 0), AniB.delta, G2);
AB.p_wilcoxon  = splitapply(@(x) signrank(x, 0, 'tail', 'right'), AniB.delta, G2);
writetable(AB, fullfile(config.PathSave, 'paires_agregation', 'paires_brutes_resultats.csv'));
end


function [MR, AniR, AR] = regroupement_anneaux(Rann, opts)
% (5 nov., demandé par l'utilisateur) Test du regroupement spatial des couleurs
% selon la distance (panneau D), sur les ANNEAUX disjoints de chaque session.
%   observé = paires de même couleur dans l'anneau ;
%   attendu = paires (colorées) de l'anneau x probabilité qu'une paire soit de
%             même couleur si les couleurs étaient mélangées dans le champ
%             (= paires de même couleur / toutes les paires, sur tous les
%             anneaux de la session) ;
%   y = log((observé + 0,5) / (attendu + 0,5)) : 0 = comme le hasard.
%   LMM : y ~ 1 + log2(distance / 100 µm) + (1|animal) + (1|session),
%   anneaux de 0 à 500 µm (l'anneau ouvert > 500 µm n'a pas de distance
%   centrale), pondéré par le nombre de paires attendues, ddl de Satterthwaite.
%   Ordonnée à l'origine = excès à 100 µm ; pente = changement quand la
%   distance double. Toutes couleurs, puis chaque couleur (descriptif).
%   AniR / AR : log(observé / attendu) par animal (moyenne de ses sessions) et
%   entre animaux, pour la figure (tous les anneaux, > 500 µm compris).
MR = table(); AniR = table(); AR = table();
if isempty(Rann) || ~ismember('n_paires_meme', Rann.Properties.VariableNames), return; end
D = Rann(strcmp(Rann.mesure, opts.mesures{1}) & strcmp(Rann.permutation, opts.permutations{1}), :);
if isempty(D), return; end
D.session_u = strcat(D.animal, '_', D.champ, '_', D.session, '_', D.condition);
Dt = D(strcmp(D.couleur, 'toutes'), :);
Dt.N = Dt.n_paires_meme + Dt.n_paires_diff;                % toutes les paires colorées de l'anneau
[gs, sess] = findgroups(Dt.session_u);
tot_N = splitapply(@sum, Dt.N, gs);
cle_t = strcat(Dt.session_u, '|', string(Dt.dist_min_um), '|', string(Dt.dist_max_um));
noms = {'toutes', 'R', 'G', 'B', 'RG', 'RB', 'GB'};
morceaux = {};
for k = 1:numel(noms)
    Dk = D(strcmp(D.couleur, noms{k}), :);
    if isempty(Dk), continue; end
    cle_k = strcat(Dk.session_u, '|', string(Dk.dist_min_um), '|', string(Dk.dist_max_um));
    [ok, it] = ismember(cle_k, cle_t);
    Dk = Dk(ok, :); it = it(ok);
    [gk, sk] = findgroups(Dk.session_u);
    meme_tot = splitapply(@sum, Dk.n_paires_meme, gk);
    [~, is] = ismember(sk, sess);
    P = meme_tot ./ tot_N(is);                              % probabilité « même couleur » au hasard
    attendu = Dt.N(it) .* P(gk);
    morceaux{end+1, 1} = table(Dk.animal, Dk.session_u, repmat(noms(k), height(Dk), 1), Dk.dist_min_um, ...
        Dk.dist_max_um, Dk.n_paires_meme, attendu, 'VariableNames', ...
        {'animal', 'session', 'couleur', 'dist_min_um', 'dist_max_um', 'observe', 'attendu'}); %#ok<AGROW>
end
if isempty(morceaux), return; end
L = vertcat(morceaux{:});
L = L(L.attendu > 0, :);
L.lr = log((L.observe + 0.5) ./ (L.attendu + 0.5));

% Figure : par animal (moyenne de ses sessions), puis entre animaux
[G, coul, dmin, dmax, anim] = findgroups(L.couleur, L.dist_min_um, L.dist_max_um, L.animal);
AniR = table(coul, dmin, dmax, anim, splitapply(@mean, L.lr, G), ...
    'VariableNames', {'couleur', 'dist_min_um', 'dist_max_um', 'animal', 'lr'});
[G2, coul, dmin, dmax] = findgroups(AniR.couleur, AniR.dist_min_um, AniR.dist_max_um);
AR = table(coul, dmin, dmax, splitapply(@mean, AniR.lr, G2), ...
    splitapply(@(x) std(x) / sqrt(numel(x)), AniR.lr, G2), splitapply(@numel, AniR.lr, G2), ...
    'VariableNames', {'couleur', 'dist_min_um', 'dist_max_um', 'lr_moyen', 'lr_sem', 'n_animaux'});

% Modèle mixte, toutes couleurs puis chaque couleur
lignes = cell(0, 9);
for k = 1:numel(noms)
    Lk = L(strcmp(L.couleur, noms{k}) & isfinite(L.dist_max_um), :);
    if height(Lk) < 6 || numel(unique(Lk.animal)) < 3, continue; end
    Lk.x = log2(((Lk.dist_min_um + Lk.dist_max_um) / 2) / 100);
    Lk.w = Lk.attendu / mean(Lk.attendu);
    Lk.animal_c = categorical(Lk.animal);
    Lk.session_c = categorical(Lk.session);
    try
        m = fitlme(Lk, 'lr ~ 1 + x + (1|animal_c) + (1|session_c)', 'Weights', Lk.w);
        [~, ~, c] = fixedEffects(m, 'DFMethod', 'satterthwaite');
        lignes(end+1, :) = {noms{k}, c.Estimate(1), c.Lower(1), c.Upper(1), c.pValue(1), ...
            c.Estimate(2), c.pValue(2), numel(unique(Lk.session)), numel(unique(Lk.animal))}; %#ok<AGROW>
    catch err
        warning('Modèle du regroupement (%s) impossible : %s', noms{k}, err.message);
    end
end
if ~isempty(lignes)
    MR = cell2table(lignes, 'VariableNames', {'couleur', 'est_100um', 'ic_bas', 'ic_haut', 'p_100um', ...
        'pente', 'p_pente', 'n_sessions', 'n_animaux'});
end
end


function L = lignes_regroupement(MR)
% Lignes du résumé : regroupement spatial des couleurs selon la distance
L = {};
if isempty(MR), return; end
L{end+1} = '';
L{end+1} = ['REGROUPEMENT SPATIAL DES COULEURS SELON LA DISTANCE (anneaux 0-500 µm ; LMM log(observé / ' ...
    'attendu) ~ log2(distance), animal et session aléatoires, pondéré par les paires attendues, ' ...
    'ddl de Satterthwaite) :'];
r = MR(strcmp(MR.couleur, 'toutes'), :);
if ~isempty(r)
    L{end+1} = sprintf(['   « Toutes couleurs : à 100 µm, %s %% de paires de même couleur par rapport au ' ...
        'hasard [IC 95 %% %s à %s], p = %s ; l''excès est multiplié par %s quand la distance double ' ...
        '(p = %s). »'], nb3(100 * (exp(r.est_100um) - 1), '%+.0f'), nb3(100 * (exp(r.ic_bas) - 1), '%+.0f'), ...
        nb3(100 * (exp(r.ic_haut) - 1), '%+.0f'), pv3(r.p_100um), nb3(exp(r.pente), '%.2f'), pv3(r.p_pente));
end
Rc = MR(~strcmp(MR.couleur, 'toutes'), :);
if ~isempty(Rc)
    t = arrayfun(@(i) sprintf('%s %s %% (p = %s)', Rc.couleur{i}, nb3(100 * (exp(Rc.est_100um(i)) - 1), '%+.0f'), ...
        pv3(Rc.p_100um(i))), 1:height(Rc), 'UniformOutput', false);
    L{end+1} = ['   Par couleur, à 100 µm (descriptif, 6 tests) : ' strjoin(t, ' ; ')];
end
end


function L = ligne_brute(AB, opts, principal_dist)
% Une ligne du résumé : le test principal avec la mesure brute (sans résidus)
L = {};
if isempty(AB), return; end
r = AB(strcmp(AB.mesure, opts.mesures{1}) & AB.dist_min_um == principal_dist(1) & ...
    AB.dist_max_um == principal_dist(2), :);
if isempty(r), return; end
L{end+1} = '';
L{end+1} = sprintf(['SANS RÉSIDUS (%s brut : distance et niveau de chaque cellule NON retirés), mêmes paires : ' ...
    'même couleur %s, couleurs différentes %s, différence %s (%s %%), %d/%d animaux > 0, Wilcoxon unilatéral ' ...
    'sur les animaux : p = %s.'], opts.mesures{1}, nb3(r.brut_meme(1), '%.4f'), nb3(r.brut_diff(1), '%.4f'), ...
    nb3(r.delta_moyen(1), '%+.4f'), nb3(r.delta_pct(1), '%+.1f'), r.n_pos(1), r.n_animaux(1), pv3(r.p_wilcoxon(1)));
end


function MD = modele_distances(Rann, opts)
% (5 nov.) Un seul test sur TOUTES les distances : modèle linéaire mixte sur
% les anneaux disjoints de chaque session (0-50, 50-100, ... µm), mesure et
% permutation principales, toutes couleurs.
%   y = delta (différence même - différente couleur, sans correction du hasard)
%   y ~ 1 + (1|animal) + (1|session), pondéré par le nombre de paires
%       (poids = 1 / (1/n_même + 1/n_différentes)) : l'ordonnée à l'origine
%       est la différence moyenne sur toutes les distances.
%   Pente : y ~ 1 + log2(distance) + ... (l'effet baisse-t-il avec la distance ?)
%   Wilcoxon sur les animaux : moyenne pondérée de y de chaque animal.
MD = struct('ok', false, 'n_lignes', 0, 'n_sessions', 0, 'n_animaux', 0, 'anneaux', '', ...
    'est', NaN, 'ic', [NaN NaN], 'p', NaN, 'pct', NaN, ...
    'pente', NaN, 'pente_p', NaN, 'p_wilcoxon', NaN, 'n_pos', NaN, 'message', '');
if isempty(Rann)
    MD.message = 'pas d''anneaux : refaire l''étape 3b (main_paires) avec config.paires.distances = ''disques''';
    return
end
D = Rann(strcmp(Rann.mesure, opts.mesures{1}) & strcmp(Rann.permutation, opts.permutations{1}) & ...
    strcmp(Rann.couleur, 'toutes') & Rann.testable, :);
D = D(isfinite(D.delta), :);
if height(D) < 5
    MD.message = 'trop peu d''anneaux testables';
    return
end
D.y = D.delta;
D.w = 1 ./ (1 ./ D.n_paires_meme + 1 ./ D.n_paires_diff);
D.w = D.w / mean(D.w);
D.animal_c  = categorical(D.animal);
D.session_c = categorical(strcat(D.animal, '_', D.champ, '_', D.session, '_', D.condition));
centre = (D.dist_min_um + D.dist_max_um) / 2;
D.x = log2(centre / 100);                     % 0 = 100 µm ; +1 = distance doublée
MD.n_lignes   = height(D);
MD.n_sessions = numel(unique(D.session_c));
MD.n_animaux  = numel(unique(D.animal_c));
u = unique([D.dist_min_um, D.dist_max_um], 'rows');
MD.anneaux = strjoin(arrayfun(@(i) strrep(sprintf('%g-%g', u(i, 1), u(i, 2)), '-Inf', '+'), ...
    1:size(u, 1), 'UniformOutput', false), ', ');
try
    % degrés de liberté de Satterthwaite : avec 8 animaux, le test par défaut
    % (résidus) donne trop de faux positifs (simulation : 7 à 12 % pour 5 %),
    % Satterthwaite (~ nb d'animaux - 1) environ 4 %
    m1 = fitlme(D, 'y ~ 1 + (1|animal_c) + (1|session_c)', 'Weights', D.w);
    [~, ~, c1] = fixedEffects(m1, 'DFMethod', 'satterthwaite');
    MD.est = c1.Estimate(1); MD.ic = [c1.Lower(1), c1.Upper(1)]; MD.p = c1.pValue(1);
    sim = sum(D.w .* D.similarite_moyenne, 'omitnan') / sum(D.w(isfinite(D.similarite_moyenne)));
    MD.pct = 100 * MD.est / sim;
    Df = D(isfinite(D.x), :);
    m2 = fitlme(Df, 'y ~ 1 + x + (1|animal_c) + (1|session_c)', 'Weights', Df.w);
    [~, ~, c2] = fixedEffects(m2, 'DFMethod', 'satterthwaite');
    MD.pente = c2.Estimate(2); MD.pente_p = c2.pValue(2);
    MD.ok = true;
catch err
    MD.message = err.message;
end
[ga, ~] = findgroups(D.animal);
ya = splitapply(@(y, w) sum(y .* w) / sum(w), D.y, D.w, ga);
MD.n_pos = sum(ya > 0);
if numel(ya) >= 2, MD.p_wilcoxon = signrank(ya, 0, 'tail', 'right'); end
MD.y_animal = ya;
end


function L = resume_paires(A, opts, principal_dist, cumulatif, n_sessions, n_animaux, p_taux, p_coup, MD)
% Les phrases du résultat de l'étape 3b, avec les tests
if nargin < 9, MD = struct('ok', false, 'message', ''); end
mes0 = opts.mesures{1};
perm0 = opts.permutations{1};
L = {};
L{end+1} = '=================================================================';
L{end+1} = 'PAIRES (étape 3b) : les cellules de même couleur sont-elles plus synchrones au repos ?';
L{end+1} = '=================================================================';
L{end+1} = sprintf(['Mesure : %s au repos, corrigé de la distance et du niveau de chaque cellule ; ' ...
    'permutation : %s. %d sessions, %d animaux ; l''animal est l''unité.'], mes0, perm0, n_sessions, n_animaux);
est = @(m, p, c, d) A(strcmp(A.mesure, m) & strcmp(A.permutation, p) & strcmp(A.couleur, c) & ...
    A.dist_min_um == d(1) & A.dist_max_um == d(2), :);
r = est(mes0, perm0, 'toutes', principal_dist);
if cumulatif
    texte_d = sprintf('paires à moins de %g µm', principal_dist(2));
else
    texte_d = sprintf('paires à %g-%g µm', principal_dist(1), principal_dist(2));
end
L{end+1} = '';
if isempty(r)
    L{end+1} = sprintf('TEST PRINCIPAL (%s) : pas de résultat (étape 3b à refaire avec ces distances).', texte_d);
else
    if r.delta_moyen(1) >= 0, sens = 'plus'; else, sens = 'moins'; end
    L{end+1} = sprintf('TEST PRINCIPAL (%s, toutes couleurs) :', texte_d);
    L{end+1} = sprintf(['   « Les paires de même couleur sont %s synchrones que les paires de couleurs ' ...
        'différentes : différence de %s %s ± %s (moyenne des %d animaux ± sem ; %s %% de la valeur ' ...
        'moyenne de ces paires) ; Wilcoxon unilatéral sur les animaux : p = %s ; permutations des ' ...
        'sessions combinées : p = %s. »'], sens, mes0, ...
        nb3(r.delta_moyen(1), '%+.4f'), nb3(r.sem(1), '%.4f'), ...
        r.n_animaux(1), nb3(r.delta_pct(1), '%+.1f'), pv3(r.p_wilcoxon(1)), pv3(r.p_combine(1)));
    L{end+1} = sprintf('   (Wilcoxon : avec %d animaux, la plus petite p possible est %s.)', ...
        r.n_animaux(1), nb3(1 / 2^r.n_animaux(1), '%.4f'));
end
% Toutes les distances : modèle mixte sur les anneaux
L{end+1} = '';
if MD.ok
    L{end+1} = sprintf(['TOUTES LES DISTANCES (complément ; anneaux disjoints %s µm ; LMM, animal et ' ...
        'session en effets aléatoires, pondéré par le nombre de paires, ddl de Satterthwaite) :'], MD.anneaux);
    L{end+1} = sprintf(['   « Différence même - différente couleur : %s [IC 95 %% %s à %s] ' ...
        '(%s %% du %s moyen), p = %s (bilatéral) ; Wilcoxon unilatéral sur les animaux ' ...
        '(moyenne de toutes les distances) : p = %s (%d/%d animaux > 0). »'], ...
        nb3(MD.est, '%+.4f'), nb3(MD.ic(1), '%+.4f'), nb3(MD.ic(2), '%+.4f'), nb3(MD.pct, '%+.1f'), ...
        opts.mesures{1}, pv3(MD.p), pv3(MD.p_wilcoxon), MD.n_pos, MD.n_animaux);
    L{end+1} = sprintf(['   L''effet change-t-il avec la distance ? pente %s par doublement de la distance, ' ...
        'p = %s (anneaux de distance finie).'], nb3(MD.pente, '%+.4f'), pv3(MD.pente_p));
elseif isfield(MD, 'message') && ~isempty(MD.message)
    L{end+1} = ['TOUTES LES DISTANCES : non calculé (' MD.message ').'];
end

% Profil selon la distance
L{end+1} = '';
L{end+1} = 'Profil selon la distance (toutes couleurs ; description, pas de conclusion propre) :';
Ad = A(strcmp(A.mesure, mes0) & strcmp(A.permutation, perm0) & strcmp(A.couleur, 'toutes'), :);
Ad = sortrows(Ad, {'dist_min_um', 'dist_max_um'});
for i = 1:height(Ad)
    if cumulatif
        td = sprintf('<= %g µm', Ad.dist_max_um(i));
    else
        td = sprintf('%g-%g µm', Ad.dist_min_um(i), Ad.dist_max_um(i));
    end
    td = strrep(td, '<= Inf µm', 'toutes  ');
    L{end+1} = sprintf('   %-10s delta %s (%s %%) ; Wilcoxon p = %s ; combiné p = %s ; paires de même couleur %d', ...
        td, nb3(Ad.delta_moyen(i), '%+.4f'), nb3(Ad.delta_pct(i), '%+.1f'), pv3(Ad.p_wilcoxon(i)), ...
        pv3(Ad.p_combine(i)), round(Ad.n_paires_meme(i))); %#ok<AGROW>
end
% Par couleur
noms = {'R', 'G', 'B', 'RG', 'RB', 'GB'};
t = {};
for c = 1:numel(noms)
    rc = est(mes0, perm0, noms{c}, principal_dist);
    if isempty(rc), continue; end
    t{end+1} = sprintf('%s %s (p = %s)', noms{c}, nb3(rc.delta_moyen(1), '%+.4f'), pv3(rc.p_wilcoxon(1))); %#ok<AGROW>
end
if ~isempty(t)
    L{end+1} = sprintf('Par couleur (%s, Wilcoxon sur les animaux) : %s', texte_d, strjoin(t, ' ; '));
end
% Sensibilité
L{end+1} = sprintf('Sensibilité (%s, toutes couleurs) :', texte_d);
mesures = liste_mesures_paires(opts);
for im = 1:numel(mesures)
    for ip = 1:numel(opts.permutations)
        if im == 1 && ip == 1, continue; end
        rs = est(mesures{im}, opts.permutations{ip}, 'toutes', principal_dist);
        if isempty(rs), continue; end
        L{end+1} = sprintf('   %s / %s : delta %s ; Wilcoxon p = %s ; combiné p = %s', mesures{im}, ...
            opts.permutations{ip}, nb3(rs.delta_moyen(1), '%+.4f'), pv3(rs.p_wilcoxon(1)), pv3(rs.p_combine(1))); %#ok<AGROW>
    end
end
% Regroupement spatial des paires de même couleur
Ag = A(strcmp(A.mesure, mes0) & strcmp(A.permutation, 'globale') & strcmp(A.couleur, 'toutes'), :);
if ~isempty(Ag)
    Ag = sortrows(Ag, {'dist_min_um', 'dist_max_um'});
    t = arrayfun(@(i) sprintf('%g µm x%s', Ag.dist_max_um(i), nb3(Ag.enrichissement(i), '%.2f')), ...
        1:height(Ag), 'UniformOutput', false);
    L{end+1} = ['Paires de même couleur observées / attendues (regroupement spatial des couleurs) : ' strjoin(t, ' ; ')];
end
if ~isempty(p_taux)
    L{end+1} = sprintf('Activité au repos différente selon la couleur (Kruskal-Wallis p < 0,05) : %d / %d sessions', ...
        sum(p_taux < 0.05), sum(~isnan(p_taux)));
end
if ~isempty(p_coup)
    L{end+1} = sprintf('Couplage au réseau différent selon la couleur (Kruskal-Wallis p < 0,05) : %d / %d sessions', ...
        sum(p_coup < 0.05), sum(~isnan(p_coup)));
end
end


function t = nb3(x, fmt)
if isnan(x), t = 'NA'; else, t = strrep(sprintf(fmt, x), '.', ','); end
end


function t = pv3(p)
if isnan(p)
    t = 'NA';
elseif p < 0.001
    t = '< 0,001';
else
    t = strrep(sprintf('%.3f', p), '.', ',');
end
end


function q = bh_q(p)
% q-values de Benjamini-Hochberg (NaN ignorés)
q = nan(size(p));
ok = find(~isnan(p));
m = numel(ok);
if m == 0, return; end
[ps, ordre] = sort(p(ok));
qs = ps .* m ./ (1:m)';
qs = min(1, cummin(qs, 'reverse'));
q(ok(ordre)) = qs;
end
