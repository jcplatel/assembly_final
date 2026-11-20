function V = validation_kfold(DF, debuts, fins, opts, graine, trames_course)
%% =========================================================================
%  VALIDATION_KFOLD (19 nov. 2026, demandée par l'utilisateur)
%  Reproductibilité de la séquence de course, validation croisée répétée :
%    - n_courses >= opts.kfold_min_courses (10) : opts.kfold_k paquets (5),
%      opts.kfold_repetitions mélanges (20) -> 100 coupures 80 / 20 ;
%    - sinon : 2 paquets (50 / 50), 100 / 2 = 50 mélanges -> 100 coupures.
%    Chaque course est testée exactement une fois par mélange.
%  Pour chaque coupure : cellules de la séquence et décalage sur les courses
%  d'apprentissage (choisir_cellules_sequence : Otsu ou test par faux débuts selon
%  opts.selection), décalage sur les courses de test ;
%  rho de Spearman (apprentissage contre test) sur les cellules de la
%  séquence d'apprentissage (au moins 5).
%  Session :
%    rho   = moyenne des rho des coupures ;
%    p_rho = permutation : cellules échangées entre elles (une fois pour toutes
%            les coupures, corrigé le 20 nov. ; avant : mélange indépendant
%            dans chaque coupure, hasard trop étroit) ; décalages de test mélangés
%            entre les cellules de la séquence, moyenne des rho refaite
%            (opts.kfold_nPerm fois) ; unilatéral (rho > 0).
%  Cellules validées : retenues sur l'apprentissage ET dans la séquence du
%  test, dans au moins opts.kfold_seuil_valide (50 %) des coupures. Leur
%  ordre = décalage moyen mesuré sur les courses de test (coupures où elles
%  sont retenues sur l'apprentissage) ; rangées comme decalage_cellules (du
%  plus proche de 0 au plus tardif).
%  20 nov. : coupures calculées en parallèle (parfor), chacune avec son
%  propre flux aléatoire (faux débuts reproductibles).
%  20 nov. (choix de l'utilisateur) : cellules robustes. Dans chaque
%  coupure, position de chaque cellule dans la séquence (rang, 0 à 1) sur
%  l'apprentissage et sur le test ; écart moyen sur les coupures. Hasard :
%  positions de test mélangées entre les cellules de la coupure
%  (opts.stabilite_nPerm fois). Robustes = validées ET écart plus petit que
%  le hasard (FDR opts.stabilite_fdr sur les validées). ordre_robuste :
%  rangées par leur décalage moyen sur les courses de test. Le rho de la
%  session reste calculé sur toutes les cellules (pas sur les robustes, qui
%  sont choisies pour garder leur place).
%  =========================================================================
if nargin < 6, trames_course = []; end
n = numel(debuts);
NCell = size(DF, 1);
if n >= opts.kfold_min_courses
    K = opts.kfold_k; nrep = opts.kfold_repetitions;
else
    K = 2; nrep = round(opts.kfold_k * opts.kfold_repetitions / 2);
end
ncoup = K * nrep;

% coupures tirées d'abord (dans l'ordre, comme avant)
fixer_flux_aleatoire(graine, 7, 1);
coup_test = cell(ncoup, 1); coup_app = cell(ncoup, 1);
c = 0;
for r = 1:nrep
    melange = randperm(n);
    paquet = mod(0:n - 1, K) + 1;
    for k = 1:K
        c = c + 1;
        coup_test{c} = sort(melange(paquet == k));
        coup_app{c} = sort(melange(paquet ~= k));
    end
end

rho = nan(ncoup, 1);
rangs_a = cell(ncoup, 1); rangs_t = cell(ncoup, 1);
n_app = zeros(NCell, 1); n_val = zeros(NCell, 1);
somme_dec_test = zeros(NCell, 1);
n_sel = nan(ncoup, 1); n_deux = nan(ncoup, 1); n_deux_hasard = nan(ncoup, 1);
% 20 nov. : décalages de chaque coupure gardés (tests par permutation plus bas)
DA = nan(NCell, ncoup); DT = nan(NCell, ncoup); SA = false(NCell, ncoup);
parfor c = 1:ncoup
    fixer_flux_aleatoire(graine, 100 + c, 100 + ncoup);     % faux débuts de cette coupure
    [~, d_app, ~, s_app] = choisir_cellules_sequence(DF, debuts, fins, coup_app{c}, opts, trames_course);
    [~, d_test, ~, s_test] = choisir_cellules_sequence(DF, debuts, fins, coup_test{c}, opts, trames_course);
    s_app = s_app(:); s_test = s_test(:); d_test = d_test(:); d_app = d_app(:);
    sel = s_app & isfinite(d_test);        % décalage de test mesurable
    DA(:, c) = d_app; DT(:, c) = d_test; SA(:, c) = sel;
    n_app = n_app + sel;
    n_val = n_val + (s_app & s_test);
    v = zeros(NCell, 1);
    v(sel) = d_test(sel);
    somme_dec_test = somme_dec_test + v;
    n_sel(c) = sum(s_app);
    n_deux(c) = sum(s_app & s_test);
    n_deux_hasard(c) = sum(s_app) * sum(s_test) / NCell;
    if sum(sel) >= 5
        ra = tiedrank(d_app(sel)); rt = tiedrank(d_test(sel));
        if std(ra) > 0 && std(rt) > 0
            rho(c) = corr(ra, rt);
            rangs_a{c} = ra; rangs_t{c} = rt;
        end
    end
end
coupures_ok = find(isfinite(rho))';
V = struct();
V.validation = sprintf('%d-fold x %d', K, nrep);
V.K = K; V.n_repetitions = nrep; V.n_coupures = ncoup;
V.rho_coupures = rho;
V.rho = mean(rho, 'omitnan');
V.n_coupures_rho = sum(isfinite(rho));
V.n_sequence = mean(n_sel);
V.n_deux = mean(n_deux);
V.n_deux_hasard = mean(n_deux_hasard);

% permutation (corrigée le 20 nov.) : les cellules sont échangées entre
% elles UNE fois pour toutes les coupures (la cellule i reçoit les
% décalages de test de la cellule pi(i) dans chaque coupure). Les coupures
% réutilisent les mêmes courses : mélanger indépendamment dans chaque
% coupure (avant) donnait un hasard trop étroit, donc une p trop petite.
% Échange parmi les cellules retenues au moins une fois.
fixer_flux_aleatoire(graine, 8, 1);
V.p_rho = NaN;
pool = find(n_app > 0);
if ~isempty(coupures_ok) && numel(pool) >= 5
    nulle = zeros(opts.kfold_nPerm, 1);
    for p = 1:opts.kfold_nPerm
        DTp = DT;
        DTp(pool, :) = DT(pool(randperm(numel(pool))), :);
        s = 0;
        for c = coupures_ok
            sel = SA(:, c) & isfinite(DTp(:, c));
            if sum(sel) < 5, continue; end
            s = s + corr(tiedrank(DA(sel, c)), tiedrank(DTp(sel, c)));
        end
        nulle(p) = s / numel(coupures_ok);
    end
    V.p_rho = (1 + sum(nulle >= V.rho - 1e-12)) / (opts.kfold_nPerm + 1);
end

% cellules validées et leur ordre (décalage moyen sur les courses de test)
V.frac_validee = n_val / ncoup;
V.validees = V.frac_validee >= opts.kfold_seuil_valide;
V.decalage_test_moyen = nan(NCell, 1);
V.decalage_test_moyen(n_app > 0) = somme_dec_test(n_app > 0) ./ n_app(n_app > 0);
idx = find(V.validees);
[~, o] = sort(V.decalage_test_moyen(idx), 'descend');
V.ordre_valide = idx(o)';
V.n_validees = numel(idx);

% 20 nov. (choix de l'utilisateur) : cellules robustes = validées ET
% position stable. Dans chaque coupure, position de la cellule (rang, 0 =
% première, 1 = dernière) parmi les cellules retenues, sur l'apprentissage
% et sur le test ; écart moyen sur les coupures (~0,33 au hasard). Hasard :
% même échange des cellules que pour p_rho (une fois pour toutes les
% coupures) ; p = part des échanges où l'écart est aussi petit ; FDR sur
% les cellules validées.
[E, nE] = ecarts_position(DA, DT, SA, coupures_ok);
V.ecart_position = E ./ nE;
V.ecart_position(nE == 0) = NaN;
nP = opts.stabilite_nPerm;
nb_inf = zeros(NCell, 1);
fixer_flux_aleatoire(graine, 11, 1);
for p = 1:nP
    DTp = DT;
    DTp(pool, :) = DT(pool(randperm(numel(pool))), :);
    [Ep, nEp] = ecarts_position(DA, DTp, SA, coupures_ok);
    nb_inf = nb_inf + (Ep ./ max(nEp, 1) <= V.ecart_position + 1e-12);
end
V.p_stabilite = (1 + nb_inf) / (nP + 1);
V.p_stabilite(nE == 0) = NaN;
V.robustes = false(NCell, 1);
iv = find(V.validees & nE > 0);
if ~isempty(iv)
    V.robustes(iv) = bh_fdr_mask(V.p_stabilite(iv)', opts.stabilite_fdr)';
end
idx = find(V.robustes);
[~, o] = sort(V.decalage_test_moyen(idx), 'descend');
V.ordre_robuste = idx(o)';
V.n_robustes = numel(idx);
end


function [E, nE] = ecarts_position(DA, DT, SA, coupures)
% Somme sur les coupures de |position apprentissage - position test| de
% chaque cellule (positions = rang parmi les cellules retenues de la
% coupure, de 0 à 1), et nombre de coupures où elle compte.
NCell = size(DA, 1);
E = zeros(NCell, 1); nE = zeros(NCell, 1);
for c = coupures
    sel = SA(:, c) & isfinite(DT(:, c));
    m = sum(sel);
    if m < 5, continue; end
    pa = (tiedrank(DA(sel, c)) - 1) / (m - 1);
    pt = (tiedrank(DT(sel, c)) - 1) / (m - 1);
    E(sel) = E(sel) + abs(pa - pt);
    nE(sel) = nE(sel) + 1;
end
end
