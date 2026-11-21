function [res, Kcol, Nref] = stats_clonales_champ3(D, opts_clonal, reference, mode_permutation, graine)
%% =========================================================================
%  BIAIS CLONAL DES ASSEMBLÉES — analyse d'UN champ d'imagerie
%  -------------------------------------------------------------------------
%  Version fonction de stats_clonales_champ2 : MÊMES calculs (min-p de
%  Westfall & Young, p_min_atteignable, testable, chi2 par permutation,
%  test global du champ, descripteurs). Changements :
%   - reference ('assemblees' ou 'champ') et mode_permutation ('globale'
%     ou 'blocs') sont des arguments ; etape3_clonal lance toutes les
%     variantes de config.clonal.variantes (la 1re = test principal)
%   - 'blocs' : les couleurs ne sont mélangées qu'entre cellules d'un même
%     carré de D.taille_bloc_px pixels -> contrôle de la proximité spatiale
%   - 'territoire' : pour CHAQUE assemblée, l'urne est faite des cellules
%     colorées de sa propre surface (enveloppe convexe de ses cellules,
%     élargie de D.marge_territoire_px) ; la composition de référence ET le
%     mélange sont ceux de cette urne -> une assemblée qui n'est qu'une
%     sélection de cellules locales ne sort pas. Une assemblée qui remplit
%     toute sa surface devient non testable (rien à échanger).
%   - 'activite' : les couleurs ne sont mélangées qu'entre cellules de même
%     niveau d'activité (5 niveaux : nombre de SCE de chaque cellule, D.activite)
%     -> contrôle du biais de détection (une couleur mieux détectée a plus
%     d'événements et peut dominer une assemblée de cellules très actives)
%     (la statistique ne change pas, seule l'hypothèse nulle change)
%   - taux_nul : probabilité, sous H0, que p_minp < seuilstat pour CETTE
%     assemblée (<= 5 %, souvent moins car p_minp est discret). L'agrégation
%     compare le nombre d'assemblées biaisées à la somme de ces taux.
%   - contrôle d'alignement exact (numel(colorcell) == nb de cellules)
%   - compacité comparée aux cellules ACTIVES (au moins un SCE), pas à
%     toutes les ROI, et calculée sur gcamp_centroid (masques suite2p)
%   - les figures sont dans figures_clonales.m
%
%  QUESTION POSÉE POUR CHAQUE ASSEMBLÉE
%    « Cette assemblée contient-elle trop de cellules d'UNE couleur, quelle
%      qu'elle soit, par rapport à un tirage au hasard dans l'urne ? »
%  L'urne : 'champ' = toutes les cellules colorées du champ ;
%           'assemblees' = les cellules colorées membres d'une assemblée
%           (insensible au biais de détection : teste si les couleurs se
%           répartissent au hasard ENTRE les assemblées)
%  =========================================================================

%% ---------- 0. DONNÉES ET PARAMÈTRES ------------------------------------
colorcell      = D.colorcell(:);
assembly       = D.assembly;
gcamp_centroid = D.gcamp_centroid;
Race           = D.Race;
PCl_final      = D.PCl_final;
score_sequence = D.score_sequence;
id_animal      = D.id_animal;
id_champ       = D.id_champ;
id_session     = D.id_session;
id_condition   = D.id_condition;

nPerm          = opts_clonal.nPerm;
seuilstat      = opts_clonal.seuilstat;
nmin_couleur   = opts_clonal.nmin_couleur;
nmin_assemblee = opts_clonal.nmin_assemblee;
% 21 nov. (demandé par l'utilisateur) : champs absents = ancien calcul
%  midp   : p_minp « mid-p » (Lancaster) : permutations ex aequo avec
%           l'observé comptées à moitié : (n_< + 0,5 (n_= + 1)) / (nPerm + 1) ;
%           taux_nul et p_min_atteignable calculés de la même façon
%  tarone : pour chaque assemblée, seules les couleurs qui PEUVENT donner
%           p < seuilstat (toutes ses cellules colorées de cette couleur)
%           entrent dans le min-p ; choisi sans regarder les comptes observés
midp   = isfield(opts_clonal, 'midp') && opts_clonal.midp;
tarone = isfield(opts_clonal, 'tarone') && opts_clonal.tarone;

nom_couleur = {'R', 'G', 'B', 'RG', 'RB', 'GB'};
nC = 6;
nA = numel(assembly);

if numel(colorcell) ~= size(Race, 1)
    error('colorcell (%d) et Race (%d cellules) ne décrivent pas les mêmes cellules.', ...
        numel(colorcell), size(Race, 1));
end

%% ---------- 1. Composition du champ (toutes les cellules) ----------------
nbR  = sum(colorcell == 1);  nbG  = sum(colorcell == 2);  nbB  = sum(colorcell == 3); %#ok<NASGU>
nbRG = sum(colorcell == 4);  nbRB = sum(colorcell == 5);  nbGB = sum(colorcell == 6); %#ok<NASGU>

%% ---------- 2. Appartenance des cellules aux assemblées ------------------
paires_asm  = [];
paires_cell = [];
for a = 1:nA
    cells_a = unique(assembly{a}(:));
    paires_asm  = [paires_asm;  repmat(a, numel(cells_a), 1)]; %#ok<AGROW>
    paires_cell = [paires_cell; cells_a];                      %#ok<AGROW>
end
nb_appart = accumarray(paires_cell, 1, [numel(colorcell) 1]);
if any(nb_appart > 1)
    warning('%d cellule(s) appartiennent à plusieurs assemblées.', sum(nb_appart > 1));
end

%% ---------- 3. Population de référence pour le tirage --------------------
est_coloree = colorcell >= 1 & colorcell <= 6;   % 7 (blanc) et 8 (aucune) exclues
switch reference
    case 'champ'
        idx_ref = find(est_coloree);
    case 'assemblees'
        idx_ref = find(est_coloree & nb_appart > 0);
    otherwise
        error('reference doit valoir ''champ'' ou ''assemblees''.');
end
Nref    = numel(idx_ref);                       % nb de cellules dans l'urne
col_ref = colorcell(idx_ref);                   % leurs couleurs (1..6)
Kcol    = accumarray(col_ref, 1, [nC 1])';      % nb de cellules de chaque couleur
couleur_valide = Kcol >= nmin_couleur;
if Nref == 0
    error('Aucune cellule colorée dans la population de référence.');
end

pos_dans_ref = zeros(numel(colorcell), 1);
pos_dans_ref(idx_ref) = 1:Nref;
garde = pos_dans_ref(paires_cell) > 0;
pA    = paires_asm(garde);
pPos  = pos_dans_ref(paires_cell(garde));
if isempty(pA)
    error('Aucune cellule colorée dans les assemblées.');
end

n_col = accumarray(pA, 1, [nA 1]);                         % cellules colorées par assemblée
n_tot = cellfun(@(x) numel(unique(x)), assembly(:));       % taille totale de l'assemblée
Obs = accumarray([pA, col_ref(pPos)], 1, [nA nC]);         % Obs(a, c)

if ~strcmp(mode_permutation, 'territoire')
    %% ---------- 4. Table des p-values hypergéométriques (calculée une fois) ---
    % Ptab(a, c, k+1) = P(au moins k cellules de couleur c dans l'assemblée a)
    nmax = max([n_col; 1]);
    Ptab = nan(nA, nC, nmax + 1);
    for a = 1:nA
        for c = 1:nC
            if couleur_valide(c)
                kv = 0:n_col(a);
                Ptab(a, c, 1:n_col(a)+1) = hygecdf(kv - 1, Nref, Kcol(c), n_col(a), 'upper');
                Ptab(a, c, 1) = 1;               % P(au moins 0 cellule) = 1
            else
                Ptab(a, c, 1:n_col(a)+1) = 1;    % couleur exclue : ne gagne jamais
            end
        end
    end
    base = reshape(1:nA*nC, nA, nC);             % indice linéaire de (a, c, k) = base + k*nA*nC
    p_brut = Ptab(base + Obs*nA*nC);             % p hypergéo de chaque couleur (non corrigée, avant Tarone)
    if tarone
        % couleur exclue pour cette assemblée si même le cas le plus favorable
        % ne donne pas p < seuilstat (ne dépend que de n_col et de Kcol)
        m_max0 = min(repmat(n_col, 1, nC), repmat(Kcol, nA, 1));
        exclue_tarone = Ptab(base + m_max0*nA*nC) >= seuilstat;
        [ia, ic] = find(exclue_tarone);
        for j = 1:numel(ia)
            Ptab(ia(j), ic(j), :) = 1;           % ne gagne jamais le min
        end
    end

    T_obs  = min(Ptab(base + Obs*nA*nC), [], 2);    % la couleur la plus surprenante (couleurs gardées)

    m_max = min(repmat(n_col, 1, nC), repmat(Kcol, nA, 1));
    T_min = min(Ptab(base + m_max*nA*nC), [], 2);    % cas le plus favorable

    Eatt = n_col * (Kcol / Nref);                % effectifs attendus
    Eatt(:, ~couleur_valide) = NaN;
    chi2_obs = sum((Obs - Eatt).^2 ./ Eatt, 2, 'omitnan');

    %% ---------- 5. Permutations : on mélange les couleurs dans l'urne --------
    % 'globale' : un seul groupe (toute l'urne) ; 'blocs' : un groupe par carré ;
    % 'activite' : un groupe par niveau d'activité (quintiles dans l'urne)
    switch mode_permutation
        case 'globale'
            groupes_bloc = {(1:Nref)'};
        case 'blocs'
            xy = gcamp_centroid(idx_ref, :);
            [~, ~, num_bloc] = unique(floor(xy / D.taille_bloc_px), 'rows');
            groupes_bloc = accumarray(num_bloc, (1:Nref)', [], @(v) {v});
        case 'activite'
            act_ref = D.activite(idx_ref);
            bornes_act = unique(prctile(act_ref, [20 40 60 80]));
            [~, ~, num_niv] = histcounts(act_ref, [-Inf, bornes_act(:)', Inf]);
            groupes_bloc = accumarray(num_niv(:), (1:Nref)', [], @(v) {v});
            groupes_bloc = groupes_bloc(~cellfun(@isempty, groupes_bloc));
        otherwise
            error('mode_permutation doit valoir ''globale'', ''blocs'', ''activite'' ou ''territoire''.');
    end

    rng(graine);
    T_perm    = zeros(nA, nPerm);
    chi2_perm = zeros(nA, nPerm);
    for b = 1:nPerm
        col_melange = col_ref;                   % même mélange pour toutes les assemblées du champ
        for g = 1:numel(groupes_bloc)
            j = groupes_bloc{g};
            col_melange(j) = col_ref(j(randperm(numel(j))));
        end
        cnt = accumarray([pA, col_melange(pPos)], 1, [nA nC]);
        T_perm(:, b)    = min(Ptab(base + cnt*nA*nC), [], 2);
        chi2_perm(:, b) = sum((cnt - Eatt).^2 ./ Eatt, 2, 'omitnan');
    end

    %% ---------- 6. p-values ---------------------------------------------------
    tol = 1e-10;
    if midp
        n_lt = sum(T_perm < T_obs * (1 - tol), 2); n_le = sum(T_perm <= T_obs * (1 + tol), 2);
        p_minp = (n_lt + 0.5 * (n_le - n_lt + 1)) / (nPerm + 1);
        n_lt = sum(T_perm < T_min * (1 - tol), 2); n_le = sum(T_perm <= T_min * (1 + tol), 2);
        p_min_atteignable = (n_lt + 0.5 * (n_le - n_lt + 1)) / (nPerm + 1);
    else
        p_minp            = (1 + sum(T_perm    <= T_obs * (1 + tol), 2)) / (nPerm + 1);
        p_min_atteignable = (1 + sum(T_perm    <= T_min * (1 + tol), 2)) / (nPerm + 1);
    end
    p_chi2            = (1 + sum(chi2_perm >= chi2_obs * (1 - tol), 2)) / (nPerm + 1);
    testable = n_col >= nmin_assemblee & p_min_atteignable < seuilstat;

    % taux_nul : fraction des permutations qui auraient elles-mêmes donné p_minp < seuilstat
    % score_minp : p_minp en score continu (z), centré-réduit sur les permutations
    %   de CETTE assemblée : 0 en moyenne et écart-type 1 sans effet, quelle que
    %   soit la taille (ex aequo comptés à moitié, « mid-p »). Sert au test de
    %   l'étape 4 qui utilise toutes les p (A5).
    % p0_sup : probabilité sans effet que p_minp > 0,5 (proportion d'assemblées
    %   clonales, étape 4 A6)
    taux_nul   = nan(nA, 1);
    taux_nul_initial = mean(T_perm < seuilstat, 2);     % P(p initial < seuilstat) sous H0 (8 nov.)
    score_minp = nan(nA, 1);
    p0_sup     = nan(nA, 1);
    bornes_p   = [0.5, nPerm - 0.5] / nPerm;
    for a = 1:nA
        [~, ~, rang] = unique(T_perm(a, :));
        rang = rang(:);
        cnt = accumarray(rang, 1);
        n_inf_egal = cumsum(cnt);                        % nb de permutations <= chaque valeur
        p_perm = n_inf_egal(rang) / nPerm;               % p_minp qu'aurait eue chaque permutation
        p_mid_perm = (n_inf_egal(rang) - 0.5 * cnt(rang)) / nPerm;
        if midp, p_perm = p_mid_perm; end                % même règle que p_minp
        taux_nul(a) = mean(p_perm < seuilstat);
        p0_sup(a)   = mean(p_perm > 0.5);
        n_lt = sum(T_perm(a, :) < T_obs(a) * (1 - tol));
        n_le = sum(T_perm(a, :) <= T_obs(a) * (1 + tol));
        p_mid_obs = (n_lt + 0.5 * (n_le - n_lt)) / nPerm;
        z0 = -norminv(min(max(p_mid_perm, bornes_p(1)), bornes_p(2)));
        zo = -norminv(min(max(p_mid_obs, bornes_p(1)), bornes_p(2)));
        if std(z0) > 0
            score_minp(a) = (zo - mean(z0)) / std(z0);
        end
    end

    % Test global du champ (somme des chi2)
    S_obs   = sum(chi2_obs);
    S_perm  = sum(chi2_perm, 1);
    p_champ = (1 + sum(S_perm >= S_obs * (1 - tol))) / (nPerm + 1);
    z_champ = (S_obs - mean(S_perm)) / std(S_perm);
    Kcol_ref = repmat(Kcol, nA, 1);          % composition de l'urne de référence, par assemblée
    n_urne_territoire = nan(nA, 1);
else
    % 'territoire' : urne propre à chaque assemblée = cellules colorées de sa
    % surface (enveloppe de ses cellules + marge) ; référence ET mélange dans
    % cette urne seulement. Question : l'assemblée a-t-elle plus d'une couleur
    % que les cellules qui occupent la même surface ? (voir minp_territoire)
    TT = minp_territoire(colorcell, assembly, gcamp_centroid, est_coloree, D.marge_territoire_px, ...
        nPerm, seuilstat, nmin_couleur, graine, midp, tarone);
    tol = 1e-10; %#ok<NASGU>
    p_minp = TT.p_minp; p_min_atteignable = TT.p_min_atteignable; p_chi2 = TT.p_chi2;
    taux_nul = TT.taux_nul; score_minp = TT.score_minp; p0_sup = TT.p0_sup;
    p_brut = TT.p_brut; Obs = TT.Obs; n_col = TT.n_col;
    Kcol_ref = TT.Kcol_loc; n_urne_territoire = TT.n_urne;
    taux_nul_initial = nan(nA, 1);
    testable = n_col >= nmin_assemblee & p_min_atteignable < seuilstat;
    couleur_valide = true(1, nC);            % couleurs exclues déjà à NaN, assemblée par assemblée
    p_champ = NaN; z_champ = NaN;            % pas de test global du champ pour cette variante
end

% Couleur la plus enrichie (annotation seulement)
p_brut(:, ~couleur_valide) = NaN;
[~, idx_dom] = min(p_brut, [], 2);
couleur_dominante = nom_couleur(idx_dom)';
couleur_dominante(n_col == 0) = {'-'};
n_couleur_dominante    = Obs(sub2ind(size(Obs), (1:nA)', idx_dom(:)));
part_couleur_dominante = n_couleur_dominante ./ n_col;
part_couleur_dominante(n_col == 0) = NaN;
part_colorees          = n_col ./ n_tot;

if numel(score_sequence) == nA
    sc = score_sequence(:);
else
    sc = nan(nA, 1);
end

%% ---------- 6b. Descripteurs des assemblées (exploration) ---------------
n_cell_champ = numel(colorcell);
pct_reseau   = 100 * n_tot / n_cell_champ;

% Étalement : distance moyenne entre cellules de l'assemblée, comparée à des
% groupes de même taille tirés parmi les cellules ACTIVES (au moins un SCE)
etalement_px    = nan(nA, 1);
compacite_pct   = nan(nA, 1);
ratio_etalement = nan(nA, 1);
actives = find(sum(Race, 2) > 0 & all(~isnan(gcamp_centroid), 2));
rng(graine + 1);
n_desc = nA;                         % opts_clonal.descripteurs = false : sautés (interface, 8 nov.)
if isfield(opts_clonal, 'descripteurs') && ~opts_clonal.descripteurs, n_desc = 0; end
for a = 1:n_desc
    cells_a = intersect(unique(assembly{a}(:)), actives);
    m_a = numel(cells_a);
    if m_a < 3, continue; end
    etalement_px(a) = mean(pdist(gcamp_centroid(cells_a, :)));
    d_alea = zeros(1000, 1);
    for b = 1:1000
        d_alea(b) = mean(pdist(gcamp_centroid(actives(randperm(numel(actives), m_a)), :)));
    end
    compacite_pct(a)   = 100 * mean(d_alea <= etalement_px(a));
    ratio_etalement(a) = etalement_px(a) / mean(d_alea);
end

% SCE : activations significatives (RACEORTHO), SCE où l'assemblée est la
% plus active, fraction moyenne de ses cellules actives dans ces SCE
n_SCE_total     = size(Race, 2);
n_SCE_signif    = nan(nA, 1);
n_SCE_dominante = nan(nA, 1);
frac_active_moy = nan(nA, 1);
if size(PCl_final, 1) == nA
    n_SCE_signif = sum(PCl_final > 0, 2);
end
S_sce = Race > 0;
frac = zeros(nA, n_SCE_total);
for a = 1:nA
    cells_a = unique(assembly{a}(:));
    if ~isempty(cells_a)
        frac(a, :) = mean(S_sce(cells_a, :), 1);
    end
end
[vmax, lab] = max(frac, [], 1);
lab(vmax == 0) = 0;
for a = 1:nA
    n_SCE_dominante(a) = sum(lab == a);
    if any(lab == a)
        frac_active_moy(a) = mean(frac(a, lab == a));
    end
end

%% ---------- 7. Table de résultats ----------------------------------------
res = table();
res.animal            = repmat({id_animal},    nA, 1);
res.champ             = repmat({id_champ},     nA, 1);
res.session           = repmat({id_session},   nA, 1);
res.condition         = repmat({id_condition}, nA, 1);
res.assemblee         = (1:nA)';
res.n_total           = n_tot;
res.n_colorees        = n_col;
res.testable          = testable;
res.p_min_atteignable = p_min_atteignable;
res.p_minp            = p_minp;      % TEST PRINCIPAL (corrigé pour les 6 couleurs)
res.p_initial         = min(p_brut, [], 2);   % p hypergéo de la couleur la plus enrichie, NON corrigé (8 nov.)
res.taux_nul_initial  = taux_nul_initial;     % P(p_initial < seuilstat) sous H0 : jusqu'à ~6 x seuilstat
res.Kcol_ref          = Kcol_ref;    % nA x 6 : composition de l'urne de référence de chaque assemblée
res.n_urne_territoire = n_urne_territoire;   % 'territoire' : cellules colorées de sa surface
res.taux_nul          = taux_nul;    % P(p_minp < seuilstat) sous H0, pour cette assemblée
res.score_minp        = score_minp;  % p_minp en score continu (0 sans effet, écart-type 1)
res.p0_sup            = p0_sup;      % P(p_minp > 0,5) sous H0
res.p_chi2            = p_chi2;      % composition globale (secondaire)
res.couleur_dominante      = couleur_dominante;
res.n_couleur_dominante    = n_couleur_dominante;
res.part_couleur_dominante = part_couleur_dominante;
res.part_colorees          = part_colorees;
res.compte            = Obs;         % nA x 6 : R G B RG RB GB
res.p_brut            = p_brut;      % nA x 6 : p hypergéo NON corrigées
res.score_sequence    = sc;
res.p_champ           = repmat(p_champ, nA, 1);
res.z_champ           = repmat(z_champ, nA, 1);
res.p_detection       = repmat(D.p_detection, nA, 1);
res.reference         = repmat({reference}, nA, 1);
res.mode_permutation  = repmat({mode_permutation}, nA, 1);
res.source_couleur    = repmat({D.source_couleur}, nA, 1);
res.nPerm             = repmat(nPerm, nA, 1);
res.K                 = repmat(D.id_K, nA, 1);
res.n_cellules_champ  = repmat(n_cell_champ, nA, 1);
res.pct_reseau        = pct_reseau;
res.etalement_px      = etalement_px;
res.compacite_pct     = compacite_pct;
res.ratio_etalement   = ratio_etalement;
res.n_SCE_total       = repmat(n_SCE_total, nA, 1);
res.n_SCE_signif      = n_SCE_signif;
res.n_SCE_dominante   = n_SCE_dominante;
res.frac_active_moy   = frac_active_moy;

fprintf('   [%s / %s] %s_%s : %d assemblées, %d testables, %d à p_minp < %.2f (attendu %.2f) | champ p = %.4f\n', ...
    id_animal, id_session, reference, mode_permutation, nA, sum(testable), ...
    sum(p_minp(testable) < seuilstat), seuilstat, sum(taux_nul(testable)), p_champ);
end


%% =========================================================================
function T = minp_territoire(colorcell, assembly, xy, est_coloree, marge_px, nPerm, seuilstat, nmin_couleur, graine, midp, tarone)
% p_minp de chaque assemblée contre les cellules colorées de SA surface.
% Même statistique que le test principal (couleur la plus surprenante,
% hypergéométrique), mais la référence est la composition de la surface et
% le hasard tire les couleurs des membres parmi les cellules de la surface.
nA = numel(assembly); nC = 6;
T.p_minp = ones(nA, 1); T.p_min_atteignable = ones(nA, 1); T.p_chi2 = ones(nA, 1);
T.taux_nul = zeros(nA, 1); T.score_minp = nan(nA, 1); T.p0_sup = nan(nA, 1);
T.p_brut = nan(nA, nC); T.Obs = zeros(nA, nC); T.n_col = zeros(nA, 1);
T.Kcol_loc = zeros(nA, nC); T.n_urne = zeros(nA, 1);
colorees = find(est_coloree(:) & all(~isnan(xy), 2));
tol = 1e-10;
bornes_p = [0.5, nPerm - 0.5] / nPerm;
rng(graine);
for a = 1:nA
    membres = unique(assembly{a}(:));
    mc = membres(est_coloree(membres));
    n = numel(mc);
    T.n_col(a) = n;
    if n == 0, continue; end
    dans = territoire(xy, membres, colorees, marge_px);
    urne = union(colorees(dans), mc);
    col_u = colorcell(urne);
    N = numel(urne);
    K = accumarray(col_u(:), 1, [nC 1])';
    T.Kcol_loc(a, :) = K;
    T.n_urne(a) = N;
    obs = accumarray(colorcell(mc), 1, [nC 1])';
    T.Obs(a, :) = obs;
    valide = K >= nmin_couleur;
    if ~any(valide), continue; end
    Pt = ones(nC, n + 1);                    % Pt(c, k+1) = P(au moins k cellules de couleur c)
    for c = find(valide)
        Pt(c, :) = hygecdf((0:n) - 1, N, K(c), n, 'upper');
        Pt(c, 1) = 1;
    end
    pb = Pt(sub2ind(size(Pt), 1:nC, obs + 1));
    pb(~valide) = NaN;
    T.p_brut(a, :) = pb;
    mmax = min(n, K);
    if tarone                                % couleurs qui ne peuvent pas donner p < seuilstat
        for c = find(valide)
            if Pt(c, mmax(c) + 1) >= seuilstat, Pt(c, :) = 1; end
        end
    end
    Tobs = min(Pt(sub2ind(size(Pt), find(valide), obs(valide) + 1)));
    Tmin = min(Pt(sub2ind(size(Pt), find(valide), mmax(valide) + 1)));
    E = n * K / N; E(~valide) = NaN;
    chi2obs = sum((obs - E).^2 ./ E, 'omitnan');
    % hasard : les n membres colorés reçoivent les couleurs de n cellules
    % tirées au hasard dans la surface
    cnt = zeros(nPerm, nC);
    for b = 1:nPerm
        cnt(b, :) = accumarray(col_u(randperm(N, n)), 1, [nC 1])';
    end
    V = Pt(sub2ind(size(Pt), repmat(1:nC, nPerm, 1), cnt + 1));
    V(:, ~valide) = Inf;
    Tp = min(V, [], 2)';
    chi2p = sum((cnt - E).^2 ./ E, 2, 'omitnan')';
    if midp
        T.p_minp(a) = (sum(Tp < Tobs * (1 - tol)) + 0.5 * (sum(abs(Tp - Tobs) <= Tobs * tol) + 1)) / (nPerm + 1);
        T.p_min_atteignable(a) = (sum(Tp < Tmin * (1 - tol)) + 0.5 * (sum(abs(Tp - Tmin) <= Tmin * tol) + 1)) / (nPerm + 1);
    else
        T.p_minp(a)            = (1 + sum(Tp <= Tobs * (1 + tol))) / (nPerm + 1);
        T.p_min_atteignable(a) = (1 + sum(Tp <= Tmin * (1 + tol))) / (nPerm + 1);
    end
    T.p_chi2(a)            = (1 + sum(chi2p >= chi2obs * (1 - tol))) / (nPerm + 1);
    [~, ~, rang] = unique(Tp);
    rang = rang(:);
    nb = accumarray(rang, 1);
    n_inf_egal = cumsum(nb);
    p_perm = n_inf_egal(rang) / nPerm;
    p_mid_perm = (n_inf_egal(rang) - 0.5 * nb(rang)) / nPerm;
    if midp, p_perm = p_mid_perm; end
    T.taux_nul(a) = mean(p_perm < seuilstat);
    T.p0_sup(a)   = mean(p_perm > 0.5);
    n_lt = sum(Tp < Tobs * (1 - tol));
    n_le = sum(Tp <= Tobs * (1 + tol));
    z0 = -norminv(min(max(p_mid_perm, bornes_p(1)), bornes_p(2)));
    zo = -norminv(min(max((n_lt + 0.5 * (n_le - n_lt)) / nPerm, bornes_p(1)), bornes_p(2)));
    if std(z0) > 0, T.score_minp(a) = (zo - mean(z0)) / std(z0); end
end
end


function dans = territoire(xy, membres, candidats, marge)
% cellules candidates situées dans la surface de l'assemblée : enveloppe
% convexe des centres de ses cellules, élargie de « marge » pixels, plus
% toutes les cellules à moins de « marge » d'un membre (assemblée alignée)
P = xy(membres, :);
P = P(all(~isnan(P), 2), :);
Q = xy(candidats, :);
dans = false(numel(candidats), 1);
if isempty(P), return; end
if size(P, 1) >= 3
    etat = warning('off', 'all');
    try
        k = convhull(P(:, 1), P(:, 2));
        ps = polybuffer(polyshape(P(k, 1), P(k, 2)), marge);
        dans = isinterior(ps, Q(:, 1), Q(:, 2));
    catch
    end
    warning(etat);
end
d = min(pdist2(Q, P), [], 2);
dans = dans(:) | d <= marge;
end

