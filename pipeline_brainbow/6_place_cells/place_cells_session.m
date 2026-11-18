function R = place_cells_session(session, dossier_sortie, opts)
%% =========================================================================
%  PLACE_CELLS_SESSION (18 nov. 2026, demandée par l'utilisateur)
%  Détection des place cells d'une session et test de biais clonal.
%  Reprend place_cells_2024_08_28.mlx, avec ces corrections :
%   - bins : 100 bins de 1 % du tapis (bords 0:0.01:1). L'ancien
%     func_on_bin(linspace(0,1,100)) donnait 99 bins décalés d'un cran, le
%     bin 1 recevant les positions NaN ;
%   - bins sans passage dans un tour : NaN (et non 0 remplacé ensuite) ;
%   - comparaison de deux place cells : l'ancienne « corrélation circulaire »
%     prenait les valeurs de ΔF comme des ANGLES (atan2(mean(sin(activité))))
%     et n'avait pas de sens. Remplacée par (a) la distance circulaire entre
%     les centres des champs (0 et 100 se rejoignent) et (b) la corrélation
%     de Pearson des courbes d'activité moyennes (TrPC) ;
%   - statistique : l'ancien Wilcoxon sur les paires comptait chaque cellule
%     dans des dizaines de paires (paires non indépendantes, p trop petit).
%     Remplacé par des permutations des couleurs entre les place cells de
%     la session (toutes, puis à l'intérieur de blocs spatiaux pour enlever
%     l'effet du voisinage) ;
%   - contrôle de distance : seulement les cellules colorées (classes 1 à 6).
%  Inchangé : F / polynôme de degré 3, Savitzky-Golay (3, 5), z-score ;
%  course = vitesse > 2 cm/s ; tours = retours de la position (> 0,9 tapis) ;
%  premier et dernier tour retirés ; TrPC lissé (gaussienne, 3 bins) ;
%  seuil par cellule = max(99e percentile de 10000 décalages circulaires
%  indépendants de chaque tour, 2 x SD de TrPClap, 1) ; champ = le pic le plus
%  haut de TrPC, largeur 2 à 30 bins.
%
%  Tests de couleur (cellules colorées 1 à 6 seulement, classes 7 et 8 exclues) :
%   A. être place cell dépend-il de la couleur ? (khi2 sur la table couleur x
%      PC, p par permutation des couleurs entre cellules colorées)
%   B. deux place cells de même couleur ont-elles des champs plus proches ?
%      stat = distance circulaire moyenne des paires de même couleur ;
%      p = permutations des couleurs entre place cells colorées (unilatéral :
%      plus proche que le hasard)
%   C. idem avec la corrélation des courbes (plus corrélées que le hasard)
%   B et C aussi « par blocs » : couleurs permutées dans des blocs de
%   opts.taille_bloc_px pixels (contrôle du voisinage clonal).
%  =========================================================================

if nargin < 3, opts = struct(); end
opts = defaut(opts, 'longueur_tapis', 192);
opts = defaut(opts, 'seuil_vitesse', 2);
opts = defaut(opts, 'n_bins', 100);
opts = defaut(opts, 'n_shuffle', 10000);
opts = defaut(opts, 'percentile', 99);
opts = defaut(opts, 'seuil_min', 1);
opts = defaut(opts, 'largeur_min', 2);
opts = defaut(opts, 'largeur_max', 30);
opts = defaut(opts, 'tours_min', 5);
opts = defaut(opts, 'n_perm', 10000);
opts = defaut(opts, 'taille_bloc_px', 128);
opts = defaut(opts, 'figures', true);
if ~isfolder(dossier_sortie), mkdir(dossier_sortie); end

R = struct('nom', session.nom, 'animal', session.id_animal, 'session', session.id_session, ...
    'opts', opts, 'statut', 'ok', 'NCell', NaN, 'n_pc', NaN, 'pct_pc', NaN);

%% activité (comme le script d'origine)
S = load(session.fichier_data, 'data_recording');
F = double(S.data_recording.F);
speed = double(S.data_recording.speed(:)');
position = double(S.data_recording.position(:)');
R.sampling_rate = S.data_recording.sampling_rate;
clear S
[NCell, Nz] = size(F);
speed = ajuster(speed, Nz); position = ajuster(position, Nz);
DF = F;
ws = warning('off', 'all');
for i = 1:NCell
    p0 = polyfit(1:Nz, DF(i, :), 3);
    DF(i, :) = DF(i, :) ./ polyval(p0, 1:Nz);
end
warning(ws);
DF = sgolayfilt(DF', 3, 5)';
DF = normalize(DF, 2, 'zscore');
clear F

%% tours, pendant la course
PosT = position / opts.longueur_tapis;
course = speed > opts.seuil_vitesse;
PosRun = PosT(course);
TrRun = DF(:, course);
[~, locs] = findpeaks(-diff(PosRun), 'MinPeakHeight', 0.9);      % retours au début du tapis
bornes = [[1, locs + 1]', [locs, numel(PosRun)]'];
nb_lap = size(bornes, 1);
bords = linspace(0, 1, opts.n_bins + 1);
n_tours = nb_lap - 2;                                              % premier et dernier retirés
R.NCell = NCell; R.n_tours = n_tours;
if n_tours < opts.tours_min
    R.statut = sprintf('trop peu de tours (%d < %d)', n_tours, opts.tours_min);
    save(fullfile(dossier_sortie, 'place_cells.mat'), '-struct', 'R');
    return
end
TrPClap = nan(n_tours, NCell, opts.n_bins);
for t = 2:nb_lap - 1
    idx = bornes(t, 1):bornes(t, 2);
    b = discretize(PosRun(idx), bords);
    ok = ~isnan(b);
    for n = 1:NCell
        TrPClap(t - 1, n, :) = accumarray(b(ok)', TrRun(n, idx(ok))', [opts.n_bins 1], @mean, NaN);
    end
end
TrPC = squeeze(mean(TrPClap, 1, 'omitnan'));
if NCell == 1, TrPC = TrPC(:)'; end
TrPC = smoothdata(TrPC, 2, 'gaussian', 3);

%% seuil par décalages circulaires indépendants de chaque tour
th = zeros(NCell, 1);
nsh = opts.n_shuffle; pct = opts.percentile; nb = opts.n_bins;
parfor n = 1:NCell
    L = squeeze(TrPClap(:, n, :));                 % tours x bins
    if size(L, 2) == 1, L = L'; end
    pk = zeros(nsh, 1);
    for s = 1:nsh
        k = randi(nb, n_tours, 1);
        Ls = L;
        for t = 1:n_tours, Ls(t, :) = circshift(L(t, :), k(t)); end
        m = smoothdata(mean(Ls, 1, 'omitnan'), 'gaussian', 3);
        p = findpeaks(m);
        pk(s) = max([p, 0]);
    end
    th(n) = prctile(pk, pct);
end

%% détection (comme le script d'origine)
pc = struct('amplitude', nan(NCell, 1), 'position', nan(NCell, 1), 'largeur', nan(NCell, 1), 'prominence', nan(NCell, 1));
seuil = zeros(NCell, 1);
for n = 1:NCell
    tmp = TrPClap(:, n, :);
    seuil(n) = max([th(n), 2 * std(tmp(:), 'omitnan'), opts.seuil_min]);
    [pks, locs, w, p] = findpeaks(TrPC(n, :), 'MinPeakProminence', seuil(n), ...
        'MinPeakWidth', opts.largeur_min, 'MaxPeakWidth', opts.largeur_max);
    if ~isempty(pks)
        [~, j] = max(pks);
        pc.amplitude(n) = pks(j); pc.position(n) = locs(j); pc.largeur(n) = w(j); pc.prominence(n) = p(j);
    end
end
est_pc = ~isnan(pc.position);
R.est_pc = est_pc; R.pc = pc; R.seuil = seuil; R.seuil_shuffle = th;
R.TrPC = TrPC; R.pct_pc = 100 * mean(est_pc); R.n_pc = sum(est_pc);

%% couleurs et positions
try
    [colorcell, ~, R.source_couleur] = charger_couleurs([session.dossier_suite2p filesep], session.fichier_colorcell);
    if numel(colorcell) ~= NCell, error('colorcell (%d) et activité (%d) : nombres différents', numel(colorcell), NCell); end
    [~, ~, cen] = charger_geometrie([session.dossier_suite2p filesep]);
catch err
    R.statut = ['pas de couleur : ' err.message];
    save(fullfile(dossier_sortie, 'place_cells.mat'), '-struct', 'R');
    return
end
R.colorcell = colorcell;
col = colorcell >= 1 & colorcell <= 6;
R.n_colorees = sum(col); R.n_pc_colorees = sum(col & est_pc);

% A. proportion de place cells par couleur
cc = colorcell(col); pp = est_pc(col);
R.pc_par_couleur = arrayfun(@(c) sum(pp(cc == c)), 1:6);
R.n_par_couleur = arrayfun(@(c) sum(cc == c), 1:6);
chi_obs = khi2(cc, pp);
rng(1);
chi_perm = zeros(opts.n_perm, 1);
for b = 1:opts.n_perm, chi_perm(b) = khi2(cc(randperm(numel(cc))), pp); end
R.p_proportion = (1 + sum(chi_perm >= chi_obs)) / (opts.n_perm + 1);

% B, C. paires de place cells colorées
ip = find(col & est_pc);
R.n_paires_meme = 0;
if numel(ip) >= 4
    c = colorcell(ip);
    pos = pc.position(ip);
    dpos = abs(pos - pos'); dpos = min(dpos, opts.n_bins - dpos);      % distance circulaire (bins)
    rc = corr(TrPC(ip, :)', 'Rows', 'pairwise');
    haut = triu(true(numel(ip)), 1);
    meme = (c == c') & haut;
    R.n_paires_meme = sum(meme(:));
    if R.n_paires_meme >= 1
        blocs = floor(cen(ip, 1) / opts.taille_bloc_px) * 100 + floor(cen(ip, 2) / opts.taille_bloc_px);
        [R.dist_meme, R.dist_diff, R.p_distance, R.p_distance_blocs] = test_paires(dpos, c, haut, blocs, opts.n_perm, 'petit');
        [R.corr_meme, R.corr_diff, R.p_correlation, R.p_correlation_blocs] = test_paires(rc, c, haut, blocs, opts.n_perm, 'grand');
        % voisinage : distance anatomique des paires et distance des champs
        danat = squareform(pdist(cen(ip, :)));
        R.rho_anatomie_champ = corr(danat(haut), dpos(haut), 'Type', 'Spearman');
    end
end
save(fullfile(dossier_sortie, 'place_cells.mat'), '-struct', 'R');

%% figure de la session
if opts.figures && R.n_pc > 0
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1100 700]);
    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
    ax = nexttile(tl);
    [~, o] = sort(pc.position(est_pc)); id = find(est_pc); id = id(o);
    M = TrPC(id, :); M = (M - min(M, [], 2)) ./ max(max(M, [], 2) - min(M, [], 2), eps);
    imagesc(ax, M); colormap(ax, hot); xlabel(ax, 'position (% du tapis)'); ylabel(ax, 'place cells (triées)');
    title(ax, sprintf('%d place cells / %d cellules (%.0f %%), %d tours', R.n_pc, NCell, R.pct_pc, n_tours));
    ax = nexttile(tl); hold(ax, 'on');
    cols = [1 0 0; 0 .7 0; 0 0 1; .9 .8 0; 1 0 1; 0 .8 .8];
    for k = 1:numel(id)
        if colorcell(id(k)) >= 1 && colorcell(id(k)) <= 6
            plot(ax, pc.position(id(k)), k, 'o', 'MarkerFaceColor', cols(colorcell(id(k)), :), 'MarkerEdgeColor', 'none');
        end
    end
    set(ax, 'YDir', 'reverse', 'XLim', [0 100], 'YLim', [0 numel(id) + 1]); xlabel(ax, 'centre du champ');
    t = sprintf('couleur des place cells (p proportion %.3g', R.p_proportion);
    if isfield(R, 'p_distance'), t = sprintf('%s ; p distance %.3g, blocs %.3g', t, R.p_distance, R.p_distance_blocs); end
    title(ax, [t ')']);
    title(tl, session.nom, 'Interpreter', 'none');
    exportgraphics(fig, fullfile(dossier_sortie, 'place_cells_session.png'), 'Resolution', 130);
    close(fig);
end
end


function [s_meme, s_diff, p, p_blocs] = test_paires(M, c, haut, blocs, n_perm, sens)
% stat = moyenne de M sur les paires de même couleur ; permutations des
% couleurs (toutes, puis dans chaque bloc spatial)
v = M(haut);
f = @(cc) mean(v(meme_couleur(cc, haut)));
s_meme = f(c); s_diff = mean(v(~meme_couleur(c, haut)));
nul = zeros(n_perm, 1); nulb = zeros(n_perm, 1);
[~, ~, gb] = unique(blocs);
for b = 1:n_perm
    nul(b) = f(c(randperm(numel(c))));
    cb = c;
    for g = 1:max(gb)
        k = find(gb == g);
        cb(k) = c(k(randperm(numel(k))));
    end
    nulb(b) = f(cb);
end
if strcmp(sens, 'petit')
    p = (1 + sum(nul <= s_meme)) / (n_perm + 1); p_blocs = (1 + sum(nulb <= s_meme)) / (n_perm + 1);
else
    p = (1 + sum(nul >= s_meme)) / (n_perm + 1); p_blocs = (1 + sum(nulb >= s_meme)) / (n_perm + 1);
end
end

function m = meme_couleur(cc, haut)
M = cc(:) == cc(:)';
m = M(haut);
end

function x = khi2(cc, pp)
T = zeros(6, 2);
for c = 1:6, T(c, 1) = sum(cc == c & pp); T(c, 2) = sum(cc == c & ~pp); end
T = T(sum(T, 2) > 0, :);
E = sum(T, 2) * sum(T, 1) / sum(T(:));
x = sum((T(:) - E(:)).^2 ./ max(E(:), eps));
end

function v = ajuster(v, Nz)
if numel(v) >= Nz, v = v(1:Nz); else, v(end+1:Nz) = NaN; end
end

function s = defaut(s, champ, valeur)
if ~isfield(s, champ) || isempty(s.(champ)), s.(champ) = valeur; end
end
