function T = test_couleur_paires(E, V, dist_um, colorcell, bornes_um, nPerm, groupes, min_paires, graine, cumulatif)
%TEST_COULEUR_PAIRES  Les paires de même couleur sont-elles plus synchrones
%   que les paires de couleurs différentes, à distance égale ?
%
%   Pour chaque distance et chaque couleur c (R G B RG RB GB) :
%     delta = moyenne de E sur les paires (c, c)
%           - moyenne de E sur les paires (c, autre couleur)
%   Ligne 'toutes' : toutes les paires de même couleur contre toutes les
%   paires de couleurs différentes. Seules les cellules colorées (1 à 6)
%   comptent ; blanc (7) et sans couleur (8) sont ignorés.
%
%   Hasard : couleurs permutées entre les cellules colorées, à l'intérieur
%   de chaque groupe (groupes = même numéro pour toutes les cellules :
%   permutation globale ; blocs spatiaux ; niveaux d'activité...).
%   p unilatéral (delta plus grand que le hasard), z par rapport au hasard,
%   delta_hasard = delta moyen des permutations,
%   enrichissement = paires de même couleur observées / attendues (avec une
%   permutation globale : regroupement spatial des couleurs, donc des clones).
%   Une case est testable si elle a au moins min_paires paires de même
%   couleur ET de couleurs différentes.
%
%   cumulatif = true (5 nov., distances_paires) : deux jeux de lignes,
%     forme 'disque'  : toutes les paires à moins de bornes_um(b+1) µm
%     forme 'anneau'  : paires entre bornes_um(b) et bornes_um(b+1) µm
%                       (disjoints : servent au modèle mixte sur toutes les
%                       distances de agregation_paires)
%     Les anneaux sont obtenus par différence des sommes des disques, dans
%     la même boucle de permutations (aucun calcul en plus).
%   cumulatif = false : tranches [bornes(b), bornes(b+1)), forme 'tranche'.
%   moyenne_meme / moyenne_diff : moyennes de E des paires de même couleur /
%   de couleurs différentes (delta = moyenne_meme - moyenne_diff).

if nargin < 10, cumulatif = false; end
nC = 6;
noms = {'R', 'G', 'B', 'RG', 'RB', 'GB', 'toutes'};
bornes_um = bornes_um(:)';
col = colorcell(:);
colorees = find(col >= 1 & col <= nC);
nb = numel(bornes_um) - 1;
Ba = discretize(dist_um, bornes_um);

W = cell(nb, 1);
M = cell(nb, 1);
for b = 1:nb
    if cumulatif
        m = V & (dist_um <= bornes_um(b + 1));
    else
        m = V & (Ba == b);
    end
    Wb = E;
    Wb(~m) = 0;
    W{b} = Wb;
    M{b} = double(m);
end

[S, N] = sommes(W, M, col, nC);
[S, N] = avec_anneaux(S, N, cumulatif);
[d_obs, n_meme, n_diff, m_meme, m_diff] = stats(S, N, nC);

fixer_flux_aleatoire(graine, 1, 1);
grp = groupes(colorees);
grp = grp(:);
ug = unique(grp);
d_perm = nan([size(d_obs), nPerm]);
n_perm = zeros([size(d_obs), nPerm]);
for k = 1:nPerm
    c2 = col;
    for ig = 1:numel(ug)
        idx = colorees(grp == ug(ig));
        c2(idx) = col(idx(randperm(numel(idx))));
    end
    [Sk, Nk] = sommes(W, M, c2, nC);
    [Sk, Nk] = avec_anneaux(Sk, Nk, cumulatif);
    [d_perm(:, :, k), n_perm(:, :, k)] = stats(Sk, Nk, nC);
end

p = (1 + sum(d_perm >= d_obs, 3)) ./ (1 + sum(~isnan(d_perm), 3));   % permutations valides seulement
delta_hasard = mean(d_perm, 3, 'omitnan');
z = (d_obs - delta_hasard) ./ std(d_perm, 0, 3, 'omitnan');
enrichissement = n_meme ./ mean(n_perm, 3);
testable = n_meme >= min_paires & n_diff >= min_paires & ~isnan(d_obs);
p(~testable) = NaN;
z(~testable) = NaN;

% Distances et forme de chaque ligne de d_obs
if cumulatif
    dmin = [zeros(1, nb), bornes_um(1:nb)];
    dmax = [bornes_um(2:end), bornes_um(2:end)];
    forme = [repmat({'disque'}, 1, nb), repmat({'anneau'}, 1, nb)];
else
    dmin = bornes_um(1:nb);
    dmax = bornes_um(2:end);
    forme = repmat({'tranche'}, 1, nb);
end

% Tableau long : une ligne par (distance, couleur)
nl = size(d_obs, 1);
[ib, ic] = ndgrid(1:nl, 1:nC + 1);
ib = ib(:); ic = ic(:);
T = table();
T.couleur         = noms(ic)';
T.forme           = forme(ib)';
T.dist_min_um     = dmin(ib)';
T.dist_max_um     = dmax(ib)';
T.n_paires_meme   = n_meme(:);
T.n_paires_diff   = n_diff(:);
T.delta           = d_obs(:);
T.delta_hasard    = delta_hasard(:);
T.moyenne_meme    = m_meme(:);
T.moyenne_diff    = m_diff(:);
T.p               = p(:);
T.z               = z(:);
T.enrichissement  = enrichissement(:);
T.testable        = testable(:);
end


function [S, N] = sommes(W, M, col, nC)
% Sommes de E et nombres de paires entre couleurs (6 x 6 x distances)
nb = numel(W);
X = zeros(numel(col), nC);
ok = find(col >= 1 & col <= nC);
X(sub2ind(size(X), ok, col(ok))) = 1;
S = zeros(nC, nC, nb);
N = zeros(nC, nC, nb);
for b = 1:nb
    S(:, :, b) = X' * (W{b} * X);
    N(:, :, b) = X' * (M{b} * X);
end
end


function [S, N] = avec_anneaux(S, N, cumulatif)
% Disques -> disques puis anneaux (différences de disques successifs)
if ~cumulatif, return; end
Sa = S - cat(3, zeros(size(S, 1), size(S, 2)), S(:, :, 1:end - 1));
Na = N - cat(3, zeros(size(N, 1), size(N, 2)), N(:, :, 1:end - 1));
S = cat(3, S, Sa);
N = cat(3, N, Na);
end


function [d, n_meme, n_diff, m_meme, m_diff] = stats(S, N, nC)
% delta (distances x 7), nombres de paires et moyennes (même couleur /
% couleurs différentes)
nb = size(S, 3);
d = nan(nb, nC + 1);
n_meme = zeros(nb, nC + 1);
n_diff = zeros(nb, nC + 1);
m_meme = nan(nb, nC + 1);
m_diff = nan(nb, nC + 1);
for b = 1:nb
    Sb = S(:, :, b);  Nb = N(:, :, b);
    sm = diag(Sb) / 2;   nm = diag(Nb) / 2;                    % paires (c, c)
    sdf = sum(Sb, 2) - diag(Sb);  ndf = sum(Nb, 2) - diag(Nb); % paires (c, autre)
    d(b, 1:nC)      = (sm ./ nm - sdf ./ ndf)';
    m_meme(b, 1:nC) = (sm ./ nm)';
    m_diff(b, 1:nC) = (sdf ./ ndf)';
    n_meme(b, 1:nC) = nm';
    n_diff(b, 1:nC) = ndf';
    tm = sum(sm);  tnm = sum(nm);
    td = (sum(Sb(:)) - trace(Sb)) / 2;  tnd = (sum(Nb(:)) - trace(Nb)) / 2;
    d(b, nC + 1)      = tm / tnm - td / tnd;
    m_meme(b, nC + 1) = tm / tnm;
    m_diff(b, nC + 1) = td / tnd;
    n_meme(b, nC + 1) = tnm;
    n_diff(b, nC + 1) = tnd;
end
end
