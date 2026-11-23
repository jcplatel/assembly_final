function S = screening_agregation(racine, G, sessions, config, opts)
%SCREENING_AGREGATION  Score clonal de chaque combinaison de réglages, et
%   ce que le hasard donnerait après le même tri (13 nov., screening_reglages).
%
%   SEULE partie du screening qui lit les couleurs.
%
%   Pour chaque combinaison (prétraitement x clusters x membres x version x
%   choix du K x critères de l'étape 4), sur les sessions du screening :
%     assemblées testables, biaisées (p_minp < seuil, test 'champ' /
%     'globale' de l'étape 3), attendues par hasard (taux_nul), et
%     EXCÈS = % biaisées - % attendu, calculé par animal puis moyenné
%     (comme le test A1 de l'étape 4). Aussi : excès toutes assemblées
%     ensemble, assemblées biaisées après FDR, score continu moyen.
%
%   Le tri par le résultat : choisir le meilleur de N réglages donne un bon
%   score même sans aucun effet. Pour le mesurer, les couleurs sont
%   mélangées dans chaque champ (opts.nNull fois, le MÊME mélange pour
%   toutes les combinaisons) et le screening entier est refait sur ces
%   couleurs mélangées : on garde le meilleur excès de chaque mélange
%   (« max-statistique », comme le min-p mais sur les réglages).
%     p_seul      : p de la combinaison, comme si elle était la seule essayée
%     p_selection : part des mélanges dont le MEILLEUR réglage fait au moins
%                   aussi bien = p corrigée pour tous les réglages essayés
%
%   Sorties dans <racine>\resultats\ : screening_combinaisons.csv (toutes),
%   screening_top.csv (les 50 premières), screening_facteurs.csv (effet moyen
%   de chaque réglage), screening_pretraitements.csv, resume_screening.txt,
%   screening.mat, figure_screening.png.

seuil          = config.clonal.seuilstat;
nmin_couleur   = config.clonal.nmin_couleur;
nmin_assemblee = config.clonal.nmin_assemblee;
% 23 nov. (demandé par l'utilisateur) : mid-p et Tarone comme le pipeline
% (config.clonal.midp / tarone ; absents = ancien calcul)
midp   = isfield(config.clonal, 'midp') && config.clonal.midp;
tarone = isfield(config.clonal, 'tarone') && config.clonal.tarone;
fprintf('Screening : p_minp %s, %s\n', ternaire(midp, 'mid-p', 'classique'), ternaire(tarone, 'avec Tarone', 'sans Tarone'));
fdr_seuil      = config.agreg.fdr_seuil;
B              = opts.nNull;
nTirages       = opts.nTirages;
tol            = 1e-10;
nC6            = 6;

nP = numel(G.pre);
nS = numel(sessions);
nR = size(G.clusters, 1);
nM = numel(G.membres);
nV = numel(G.versions);
nCh = numel(G.choix_K);
nmc = numel(G.min_cellules);
nms = numel(G.min_SCE);
dossier_res = fullfile(racine, 'resultats');
if ~isfolder(dossier_res), mkdir(dossier_res); end

%% ---------- 1. Prétraitements complets (toutes les sessions faites) -------
fichiers = cell(nP, nS);
complet = true(1, nP);
for p = 1:nP
    for s = 1:nS
        fichiers{p, s} = fullfile(racine, 'variantes', G.pre(p).nom, [sessions(s).nom '.mat']);
        if ~isfile(fichiers{p, s}), complet(p) = false; end
    end
end
iP = find(complet);
if isempty(iP)
    error('Screening : aucun réglage de prétraitement n''est fini pour toutes les sessions.');
end
fprintf('\n=== SCREENING, agrégation : %d réglage(s) de prétraitement sur %d, %d sessions ===\n', ...
    numel(iP), nP, nS);

%% ---------- 2. Couleurs, ROI retirées, mélanges communs --------------------
[animaux, ~, animal_s] = unique({sessions.id_animal});
nAn = numel(animaux);
couleurs = cell(1, nS);
retirees = cell(nP, nS);
perms = cell(1, nS);
for s = 1:nS
    fc = '';
    if isfield(sessions, 'fichier_colorcell') && ischar(sessions(s).fichier_colorcell)
        fc = sessions(s).fichier_colorcell;
    end
    cc = charger_couleurs(avec_separateur(sessions(s).dossier_suite2p), fc);
    cc = cc(:);
    couleurs{s} = cc;
    retiree_qq = false(numel(cc), 1);
    for p = iP
        retirees{p, s} = roi_retirees(fullfile(G.pre(p).dossier, sessions(s).nom), numel(cc));
        retiree_qq(retirees{p, s}) = true;
    end
    % mélange entre les cellules colorées gardées par TOUS les prétraitements
    % (les rares autres gardent leur couleur) : même mélange pour tous les réglages
    commun = find(cc >= 1 & cc <= nC6 & ~retiree_qq);
    rng(mod(sessions(s).graine, 2^31 - 1000) + 777);
    PERM = repmat(uint8(cc), 1, B);
    for b = 1:B
        PERM(commun, b) = uint8(cc(commun(randperm(numel(commun)))));
    end
    perms{s} = PERM;
end

%% ---------- 3. Chaque prétraitement ----------------------------------------
nComb_p = nR * nM * nV * nCh * nmc * nms;
C = struct();                    % une ligne par combinaison (tous prétraitements)
champs_num = {'ip', 'r', 'm', 'v', 'c', 'imc', 'ims', 'n_sessions', 'n_asm', 'n_test', 'taille_med', ...
    'n_sig', 'attendu', 'exces', 'exces_pool', 'n_an_pos', 'n_fdr', 'score', 'eligible', 'p_seul'};
for k = 1:numel(champs_num), C.(champs_num{k}) = zeros(0, 1); end
C.exces_animal = zeros(0, nAn);
C.nt_animal = zeros(0, nAn);
max_nul = -inf(1, B);            % meilleur excès de chaque mélange (toutes combinaisons éligibles)
noms_fact = {'pretraitement', 'clusters', 'membres', 'version', 'choix_K', 'min_cellules', 'min_SCE'};
n_niv = [nP, nR, nM, nV, nCh, nmc, nms];
F = struct();
for f = 1:numel(noms_fact)
    F(f).nom = noms_fact{f};
    F(f).somme = zeros(n_niv(f), 1);
    F(f).n = zeros(n_niv(f), 1);
    F(f).meilleur = -inf(n_niv(f), 1);
    F(f).somme_nul = zeros(n_niv(f), B);
end
Tpre = struct('nom', {}, 'facteur', {}, 'niveau', {}, 'n_SCE_median', {}, 'n_cellules_median', {});

for ip = 1:numel(iP)
    p = iP(ip);
    t0 = tic;
    % --- 3a. toutes les assemblées de tous les jeux, session par session
    a_sess = []; a_ntot = []; a_nsce = []; a_test = []; a_taux = []; a_p = []; a_score = []; a_sig = [];
    SIGnul = false(0, B);
    Rs = cell(1, nS); deb = cell(1, nS); nA_e = cell(1, nS); decal = zeros(1, nS);
    nSCE_s = nan(1, nS); nCell_s = nan(1, nS);
    for s = 1:nS
        Rr = load(fichiers{p, s}, 'R');
        R = Rr.R;
        Rs{s} = rmfield(R, 'ens');
        cc = couleurs{s};
        if numel(cc) ~= R.NCell
            error('Session %s : colorcell (%d) et Race (%d cellules) ne décrivent pas les mêmes cellules.', ...
                sessions(s).nom, numel(cc), R.NCell);
        end
        nSCE_s(s) = R.NRace; nCell_s(s) = R.NCell;
        cc(retirees{p, s}) = 8;                      % comme couleurs_sans_doublons
        coloree = cc >= 1 & cc <= nC6;
        Kcol = accumarray(cc(coloree), 1, [nC6 1])';
        nA_e{s} = arrayfun(@(e) numel(e.cellules), R.ens);
        deb{s} = cumsum([0, nA_e{s}(1:end - 1)]);
        tot = sum(nA_e{s});
        decal(s) = numel(a_sess);
        if tot == 0, continue; end
        lig = cell(1, tot); col = cell(1, tot);
        n_tot = zeros(tot, 1); n_sce = zeros(tot, 1);
        k = 0;
        for e = 1:numel(R.ens)
            for a = 1:nA_e{s}(e)
                k = k + 1;
                cel = unique(R.ens(e).cellules{a}(:));
                n_tot(k) = numel(cel);
                n_sce(k) = R.ens(e).n_SCE(a);
                cel = cel(coloree(cel));
                lig{k} = repmat(k, numel(cel), 1);
                col{k} = cel;
            end
        end
        lig = vertcat(lig{:}); col = vertcat(col{:});
        if isempty(lig)
            M = sparse(tot, R.NCell);
        else
            M = sparse(lig, col, 1, tot, R.NCell);    % membres colorés de chaque assemblée
        end
        Obs = zeros(tot, nC6);
        for c = 1:nC6
            Obs(:, c) = full(M * double(cc == c));
        end
        n_col = sum(Obs, 2);
        L = table_minp(Kcol, max([n_col; 1]), nTirages, seuil, nmin_couleur, ...
            mod(sessions(s).graine, 2^31 - 1000) + 31 * p, midp, tarone);
        nmax = size(L.PT, 1);
        nn = min(max(n_col, 1), nmax);                % n_col = 0 : ligne factice, annulée plus bas
        T = ones(tot, 1);
        for c = 1:nC6
            T = min(T, L.PT(nn + (c - 1) * nmax + min(Obs(:, c), nmax) * (nmax * nC6)));
        end
        pm = ones(tot, 1); sc = nan(tot, 1);
        bornes_p = [0.5, nTirages - 0.5] / nTirages;
        for k = find(n_col >= 1)'
            ts = L.Ts(:, nn(k));
            n_le = sum(ts <= T(k) * (1 + tol));
            n_lt = sum(ts < T(k) * (1 - tol));
            if midp
                pm(k) = (n_lt + 0.5 * (n_le - n_lt + 1)) / (nTirages + 1);
            else
                pm(k) = (1 + n_le) / (nTirages + 1);
            end
            if L.z_sd(nn(k)) > 0
                p_mid = (n_lt + 0.5 * (n_le - n_lt)) / nTirages;
                zo = sqrt(2) * erfcinv(2 * min(max(p_mid, bornes_p(1)), bornes_p(2)));
                sc(k) = (zo - L.z_moy(nn(k))) / L.z_sd(nn(k));
            end
        end
        testable = n_col >= nmin_assemblee & L.p_min_att(nn) < seuil;
        % mêmes assemblées, couleurs mélangées : biaisée ou non dans chaque mélange
        Tn = ones(tot, B);
        for c = 1:nC6
            Cnt = full(M * double(perms{s} == c));    % tot x B
            Tn = min(Tn, L.PT(repmat(nn, 1, B) + (c - 1) * nmax + min(Cnt, nmax) * (nmax * nC6)));
        end
        sig_nul = (Tn * (1 + tol) < repmat(L.borne(nn), 1, B)) & repmat(n_col >= 1, 1, B);
        a_sess  = [a_sess; repmat(s, tot, 1)]; %#ok<AGROW>
        a_ntot  = [a_ntot; n_tot]; %#ok<AGROW>
        a_nsce  = [a_nsce; n_sce]; %#ok<AGROW>
        a_test  = [a_test; testable]; %#ok<AGROW>
        a_taux  = [a_taux; L.taux_nul(nn)]; %#ok<AGROW>
        a_p     = [a_p; pm]; %#ok<AGROW>
        a_score = [a_score; sc]; %#ok<AGROW>
        a_sig   = [a_sig; pm < seuil]; %#ok<AGROW>
        SIGnul  = [SIGnul; sig_nul]; %#ok<AGROW>
    end
    nAsm = numel(a_sess);
    Tpre(end+1) = struct('nom', G.pre(p).nom, 'facteur', G.pre(p).facteur, 'niveau', G.pre(p).niveau, ...
        'n_SCE_median', median(nSCE_s), 'n_cellules_median', median(nCell_s)); %#ok<AGROW>

    % --- 3b. combinaisons : quelles assemblées comptent pour chacune
    cr = zeros(nComb_p, 1); cm = cr; cv = cr; cch = cr; cimc = cr; cims = cr;
    c_nsess = cr; c_nasm = cr; c_taille = nan(nComb_p, 1); c_nfdr = cr;
    lig = cell(nComb_p, 1); col = cell(nComb_p, 1);
    ic = 0;
    for r = 1:nR
        for m = 1:nM
            for v = 1:nV
                for c = 1:nCh
                    idx = zeros(1, 0);
                    for s = 1:nS
                        e = jeu_choisi(Rs{s}, G.choix_K(c), r, m, v, c);
                        if e > 0
                            idx = [idx, decal(s) + deb{s}(e) + (1:nA_e{s}(e))]; %#ok<AGROW>
                        end
                    end
                    for imc = 1:nmc
                        for ims = 1:nms
                            ic = ic + 1;
                            cr(ic) = r; cm(ic) = m; cv(ic) = v; cch(ic) = c; cimc(ic) = imc; cims(ic) = ims;
                            c_nasm(ic) = numel(idx);
                            if isempty(idx), continue; end
                            g = idx(a_test(idx)' & a_ntot(idx)' >= G.min_cellules(imc) & ...
                                a_nsce(idx)' >= G.min_SCE(ims));
                            if isempty(g), continue; end
                            lig{ic} = repmat(ic, numel(g), 1);
                            col{ic} = g(:);
                            c_nsess(ic) = numel(unique(a_sess(g)));
                            c_taille(ic) = median(a_ntot(g));
                            c_nfdr(ic) = sum(bh_fdr_mask(a_p(g)', fdr_seuil));
                        end
                    end
                end
            end
        end
    end
    lig = vertcat(lig{:}); col = vertcat(col{:});
    if isempty(lig)
        W = sparse(nComb_p, max(nAsm, 1));
    else
        W = sparse(lig, col, 1, nComb_p, max(nAsm, 1));
    end

    % --- 3c. scores observés et mélangés
    a_animal = zeros(nAsm, 1);
    if nAsm > 0, a_animal = reshape(animal_s(a_sess), [], 1); end
    exA = nan(nComb_p, nAn); ntA = zeros(nComb_p, nAn);
    Enul = zeros(nComb_p, B); nb_an = zeros(nComb_p, 1);
    n_test = zeros(nComb_p, 1); n_sig = n_test; att = n_test; som_score = n_test; nb_score = n_test;
    for a = 1:nAn
        sel = find(a_animal == a);
        if isempty(sel), continue; end
        Wa = W(:, sel);
        nt = full(sum(Wa, 2));
        ns = full(Wa * double(a_sig(sel)));
        at = full(Wa * a_taux(sel));
        ok = nt > 0;
        exA(ok, a) = (ns(ok) - at(ok)) ./ nt(ok);
        ntA(:, a) = nt;
        Ns = full(Wa * double(SIGnul(sel, :)));       % nComb x B
        Enul(ok, :) = Enul(ok, :) + (Ns(ok, :) - repmat(at(ok), 1, B)) ./ repmat(nt(ok), 1, B);
        nb_an(ok) = nb_an(ok) + 1;
        n_test = n_test + nt; n_sig = n_sig + ns; att = att + at;
        scs = a_score(sel);
        sans = isnan(scs);
        scs(sans) = 0;
        som_score = som_score + full(Wa * scs);
        nb_score = nb_score + full(Wa * double(~sans));
    end
    exces = moy_nan(exA);
    Enul = Enul ./ repmat(max(nb_an, 1), 1, B);
    eligible = n_test >= opts.min_testables & all(ntA >= opts.min_par_animal, 2);
    p_seul = (1 + sum(Enul >= repmat(exces, 1, B) - 1e-12, 2)) / (B + 1);
    p_seul(isnan(exces)) = NaN;
    if any(eligible)
        max_nul = max(max_nul, max(Enul(eligible, :), [], 1));
    end
    % effet moyen de chaque niveau (combinaisons éligibles)
    niveaux = {repmat(p, nComb_p, 1), cr, cm, cv, cch, cimc, cims};
    for f = 1:numel(noms_fact)
        for q = unique(niveaux{f}(eligible))'
            sel = eligible & niveaux{f} == q;
            F(f).somme(q) = F(f).somme(q) + sum(exces(sel));
            F(f).n(q) = F(f).n(q) + sum(sel);
            F(f).meilleur(q) = max(F(f).meilleur(q), max(exces(sel)));
            F(f).somme_nul(q, :) = F(f).somme_nul(q, :) + sum(Enul(sel, :), 1);
        end
    end
    C.ip = [C.ip; repmat(p, nComb_p, 1)];
    C.r = [C.r; cr]; C.m = [C.m; cm]; C.v = [C.v; cv]; C.c = [C.c; cch];
    C.imc = [C.imc; cimc]; C.ims = [C.ims; cims];
    C.n_sessions = [C.n_sessions; c_nsess]; C.n_asm = [C.n_asm; c_nasm];
    C.n_test = [C.n_test; n_test]; C.taille_med = [C.taille_med; c_taille];
    C.n_sig = [C.n_sig; n_sig]; C.attendu = [C.attendu; att];
    C.exces = [C.exces; exces];
    C.exces_pool = [C.exces_pool; (n_sig - att) ./ max(n_test, 1)];
    C.n_an_pos = [C.n_an_pos; sum(exA > 0, 2)];
    C.n_fdr = [C.n_fdr; c_nfdr];
    C.score = [C.score; som_score ./ max(nb_score, 1)];
    C.eligible = [C.eligible; double(eligible)];
    C.p_seul = [C.p_seul; p_seul];
    C.exces_animal = [C.exces_animal; exA];
    C.nt_animal = [C.nt_animal; ntA];
    fprintf('  %-28s : %5d assemblées, %5d combinaisons éligibles sur %d, meilleur excès %+.1f points (%.0f s)\n', ...
        G.pre(p).nom, nAsm, sum(eligible), nComb_p, 100 * max([exces(eligible); -Inf]), toc(t0));
end

%% ---------- 4. Correction pour la sélection, classement --------------------
nComb = numel(C.exces);
el = C.eligible > 0;
if ~any(el)
    warning(['Screening : aucune combinaison éligible (au moins %d assemblées testables, %d par animal). ' ...
        'Baisser opts.min_testables / opts.min_par_animal.'], opts.min_testables, opts.min_par_animal);
end
C.p_selection = nan(nComb, 1);
for k = find(el)'
    C.p_selection(k) = (1 + sum(max_nul >= C.exces(k) - 1e-12)) / (B + 1);
end
cle = C.exces; cle(~el | isnan(cle)) = -Inf;
[~, ordre] = sortrows([-cle, -C.n_test]);
C.rang = zeros(nComb, 1);
C.rang(ordre) = 1:nComb;
C.rang(~el) = NaN;

for f = 1:numel(F)
    F(f).moyenne = F(f).somme ./ F(f).n;
    moy_nul = F(f).somme_nul ./ repmat(F(f).n, 1, B);
    F(f).nul_bas = prctile(moy_nul, 2.5, 2);
    F(f).nul_haut = prctile(moy_nul, 97.5, 2);
    F(f).p = (1 + sum(moy_nul >= repmat(F(f).moyenne, 1, B) - 1e-12, 2)) / (B + 1);
    F(f).p(F(f).n == 0) = NaN;
    F(f).meilleur(F(f).n == 0) = NaN;
    switch F(f).nom
        case 'pretraitement', F(f).niveaux = {G.pre.nom};
        case 'clusters',      F(f).niveaux = arrayfun(@(i) sprintf('%s %g', G.clusters{i, 1}, G.clusters{i, 2}), ...
                1:nR, 'UniformOutput', false);
        case 'membres',       F(f).niveaux = {G.membres.nom};
        case 'version',       F(f).niveaux = G.versions;
        case 'choix_K',       F(f).niveaux = {G.choix_K.nom};
        case 'min_cellules',  F(f).niveaux = arrayfun(@(x) sprintf('%g', x), G.min_cellules, 'UniformOutput', false);
        case 'min_SCE',       F(f).niveaux = arrayfun(@(x) sprintf('%g', x), G.min_SCE, 'UniformOutput', false);
    end
end

S = struct();
S.version = version_pipeline();
S.date = datestr(now, 'yyyy-mm-dd HH:MM');
S.G = G; S.C = C; S.F = F; S.Tpre = Tpre;
S.animaux = animaux;
S.sessions = {sessions.nom};
S.pre_complets = {G.pre(iP).nom};
S.pre_incomplets = {G.pre(~complet).nom};
S.max_nul = max_nul;
S.opts = rmfield_si(opts, {'sessions', 'grille'});
S.seuil = seuil; S.fdr_seuil = fdr_seuil;
S.config_base = config;
S.midp = midp;
S.tarone = tarone;

%% ---------- 5. Fichiers ------------------------------------------------------
ecrire_combinaisons(fullfile(dossier_res, 'screening_combinaisons.csv'), S, 1:nComb);
top = ordre(1:min(50, sum(el)));
ecrire_combinaisons(fullfile(dossier_res, 'screening_top.csv'), S, top);
ecrire_facteurs(fullfile(dossier_res, 'screening_facteurs.csv'), S);
ecrire_pretraitements(fullfile(dossier_res, 'screening_pretraitements.csv'), S);
lignes = resume_screening(S, ordre);
fid = fopen(fullfile(dossier_res, 'resume_screening.txt'), 'w', 'n', 'UTF-8');
fprintf(fid, '%s\n', lignes{:});
fclose(fid);
fprintf('\n');
fprintf('%s\n', lignes{:});
S.resume = lignes;
save(fullfile(dossier_res, 'screening.mat'), 'S');
try
    figure_screening(S, dossier_res);
catch err
    warning('figure_screening : %s', err.message);
end
fprintf('\nRésultats du screening : %s\n', dossier_res);
end


%% =========================================================================
function e = jeu_choisi(R, ch, r, m, v, c)
% Numéro du jeu d'assemblées d'une session pour un choix du K (0 = aucun).
e = 0;
switch ch.type
    case 'fixe'
        iK = find(R.liste_K == ch.K, 1);
        if ~isempty(iK), e = R.map_K(iK, r, m, v); end
    case 'exces'            % K où la silhouette dépasse le plus le hasard ; à égalité, le plus petit
        es = reshape(R.map_K(:, r, m, v), [], 1);
        ex = R.exces(:, r);
        ex(es == 0 | isnan(ex)) = -Inf;
        [mx, i] = max(ex);
        if isfinite(mx), e = es(i); end
    case 'consensus'
        e = R.map_cons(r, m, v, c);
end
end


function r = roi_retirees(dossier_session, n)
% ROI retirées à l'étape 1 (doublons, tailles) : non colorées, comme dans
% couleurs_sans_doublons.
r = zeros(1, 0);
f = fullfile(dossier_session, 'pretraitement.mat');
if ~isfile(f) || ~ismember('diagnostic', who('-file', f)), return; end
D = load(f, 'diagnostic');
if isfield(D.diagnostic, 'doublons') && isfield(D.diagnostic.doublons, 'retirees')
    r = [r, D.diagnostic.doublons.retirees(:)'];
end
if isfield(D.diagnostic, 'tailles_roi') && isfield(D.diagnostic.tailles_roi, 'retirees')
    r = [r, D.diagnostic.tailles_roi.retirees(:)'];
end
% population (15 nov.) : cellules écartées par config.pre.population
if isfield(D.diagnostic, 'population') && isfield(D.diagnostic.population, 'retirees')
    r = [r, D.diagnostic.population.retirees(:)'];
end
r = unique(r(r >= 1 & r <= n));
end


function m = moy_nan(X)
% moyenne de chaque ligne sans les NaN (NaN si la ligne est vide)
ok = ~isnan(X);
X(~ok) = 0;
m = sum(X, 2) ./ sum(ok, 2);
end


function s = rmfield_si(s, noms)
for k = 1:numel(noms)
    if isfield(s, noms{k}), s = rmfield(s, noms{k}); end
end
end


function [noms, cols] = colonnes_combinaisons(S, idx)
C = S.C; G = S.G;
idx = idx(:);
pct = @(x) round(1000 * x) / 10;
noms = {'rang', 'pretraitement', 'facteur_pre', 'niveau_pre', 'regle_clusters', 'percentile_clusters', ...
    'membres', 'version', 'choix_K', 'min_cellules', 'min_SCE', 'eligible', 'n_sessions', ...
    'n_assemblees', 'n_testables', 'taille_mediane', 'n_biaisees', 'pct_biaisees', 'pct_attendu', ...
    'exces_points', 'exces_toutes_assemblees_points', 'n_animaux_exces_positif', 'n_biaisees_FDR', ...
    'score_moyen', 'p_seul', 'p_selection'};
cols = {C.rang(idx), {G.pre(C.ip(idx)).nom}', {G.pre(C.ip(idx)).facteur}', {G.pre(C.ip(idx)).niveau}', ...
    G.clusters(C.r(idx), 1), cell2mat(G.clusters(C.r(idx), 2)), {G.membres(C.m(idx)).nom}', ...
    reshape(G.versions(C.v(idx)), [], 1), {G.choix_K(C.c(idx)).nom}', reshape(G.min_cellules(C.imc(idx)), [], 1), ...
    reshape(G.min_SCE(C.ims(idx)), [], 1), C.eligible(idx), C.n_sessions(idx), C.n_asm(idx), C.n_test(idx), ...
    C.taille_med(idx), C.n_sig(idx), pct(C.n_sig(idx) ./ max(C.n_test(idx), 1)), ...
    pct(C.attendu(idx) ./ max(C.n_test(idx), 1)), pct(C.exces(idx)), pct(C.exces_pool(idx)), ...
    C.n_an_pos(idx), C.n_fdr(idx), round(1000 * C.score(idx)) / 1000, C.p_seul(idx), C.p_selection(idx)};
for a = 1:numel(S.animaux)
    noms{end+1} = ['exces_points_' S.animaux{a}]; %#ok<AGROW>
    cols{end+1} = pct(C.exces_animal(idx, a)); %#ok<AGROW>
    noms{end+1} = ['n_testables_' S.animaux{a}]; %#ok<AGROW>
    cols{end+1} = C.nt_animal(idx, a); %#ok<AGROW>
end
end


function ecrire_combinaisons(fichier, S, idx)
[noms, cols] = colonnes_combinaisons(S, idx);
ecrire_csv(fichier, noms, cols);
end


function ecrire_facteurs(fichier, S)
fa = {}; ni = {}; n = []; mo = []; me = []; nb = []; nh = []; pp = [];
for f = 1:numel(S.F)
    for q = 1:numel(S.F(f).niveaux)
        fa{end+1, 1} = S.F(f).nom; %#ok<AGROW>
        ni{end+1, 1} = S.F(f).niveaux{q}; %#ok<AGROW>
        n(end+1, 1) = S.F(f).n(q); %#ok<AGROW>
        mo(end+1, 1) = 100 * S.F(f).moyenne(q); %#ok<AGROW>
        me(end+1, 1) = 100 * S.F(f).meilleur(q); %#ok<AGROW>
        nb(end+1, 1) = 100 * S.F(f).nul_bas(q); %#ok<AGROW>
        nh(end+1, 1) = 100 * S.F(f).nul_haut(q); %#ok<AGROW>
        pp(end+1, 1) = S.F(f).p(q); %#ok<AGROW>
    end
end
arr = @(x) round(100 * x) / 100;
ecrire_csv(fichier, {'reglage', 'niveau', 'n_combinaisons_eligibles', 'exces_moyen_points', ...
    'meilleur_exces_points', 'hasard_bas_points', 'hasard_haut_points', 'p_exces_moyen'}, ...
    {fa, ni, n, arr(mo), arr(me), arr(nb), arr(nh), pp});
end


function ecrire_pretraitements(fichier, S)
T = S.Tpre;
f1 = S.F(1);
[~, i] = ismember({T.nom}, f1.niveaux);
arr = @(x) round(100 * x) / 100;
ecrire_csv(fichier, {'pretraitement', 'facteur', 'niveau', 'n_SCE_median', 'n_cellules_median', ...
    'n_combinaisons_eligibles', 'exces_moyen_points', 'meilleur_exces_points', 'hasard_bas_points', ...
    'hasard_haut_points'}, ...
    {{T.nom}', {T.facteur}', {T.niveau}', [T.n_SCE_median]', [T.n_cellules_median]', f1.n(i), ...
    arr(100 * f1.moyenne(i)), arr(100 * f1.meilleur(i)), arr(100 * f1.nul_bas(i)), arr(100 * f1.nul_haut(i))});
end


function ecrire_csv(fichier, noms, cols)
% CSV simple : colonnes numériques (vecteurs) ou texte (cellstr)
n = numel(cols{1});
fid = fopen(fichier, 'w', 'n', 'UTF-8');
if fid < 0, error('Impossible d''écrire %s (fichier ouvert dans Excel ?).', fichier); end
fprintf(fid, '%s\n', strjoin(noms, ','));
txt = cell(n, numel(cols));
for j = 1:numel(cols)
    cj = cols{j};
    if iscell(cj)
        for i = 1:n
            t = char(cj{i});
            if any(t == ',' | t == '"'), t = ['"' strrep(t, '"', '""') '"']; end
            txt{i, j} = t;
        end
    else
        for i = 1:n
            if isnan(cj(i)), txt{i, j} = ''; else, txt{i, j} = sprintf('%.6g', cj(i)); end
        end
    end
end
for i = 1:n
    fprintf(fid, '%s\n', strjoin(txt(i, :), ','));
end
fclose(fid);
end


function t = ternaire(c, a, b)
if c, t = a; else, t = b; end
end
