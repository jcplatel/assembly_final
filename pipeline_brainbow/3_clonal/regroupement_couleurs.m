function G = regroupement_couleurs(xy_um, couleur, nshuffle, graine, rayon_um)
% 12 nov. (demande de l'utilisateur : une seule mesure dans la figure) :
% aussi un COMPTAGE, la mesure des panneaux D et E de la figure :
%   obs_attendu   = nombre de paires de même couleur à moins de rayon_um
%                   (défaut 100 µm) / nombre moyen après mélange des couleurs
%                   (1 = hasard ; > 1 = cellules de même couleur regroupées)
%   z_rayon, p_rayon : z et p unilatéral (plus de paires que le hasard)
%   obs_attendu_c, z_rayon_c, p_rayon_c : paires c-c, pour chaque couleur
% L'indice de distance ci-dessous (< 1 = regroupées) reste calculé.
%REGROUPEMENT_COULEURS  Dans un champ, les cellules de même couleur sont-elles
%   plus proches entre elles que de cellules d'autres couleurs ?
%   Aucune activité, aucune assemblée : seulement les positions et les couleurs.
%
%   xy_um   : positions des cellules (µm), une ligne par cellule
%   couleur : 1 à 6 (R G B RG RB GB) ; les autres valeurs (non colorées,
%             blanches, NaN) sont ignorées
%   Toutes couleurs ensemble :
%     rapport = distance moyenne des paires de MÊME couleur
%             / distance moyenne des paires de couleurs DIFFÉRENTES
%   Pour chaque couleur c :
%     rapport_c = distance moyenne entre deux cellules c
%               / distance moyenne entre une cellule c et une cellule d'une autre couleur
%   Hasard : couleurs mélangées entre les cellules colorées du champ
%   (nshuffle fois ; positions et nombre de cellules de chaque couleur gardés).
%   indice = rapport observé / rapport moyen du hasard
%            1 = comme le hasard ; 0,8 = paires de même couleur 20 % plus proches
%   p      = (1 + nb de mélanges aussi regroupés) / (nshuffle + 1)  (unilatéral)
%   z      = (rapport observé - moyenne du hasard) / écart-type du hasard
%   Une couleur de moins de 2 cellules : NaN.
%   (Même question que le test par couleur de distance_calculation2, avec
%   en plus « même couleur contre couleurs différentes », sous forme d'indice
%   comparable d'un champ à l'autre pour l'agrégation sur les animaux.)

nC = 6;
if nargin < 5 || isempty(rayon_um), rayon_um = 100; end
couleur = couleur(:);
ok = couleur >= 1 & couleur <= nC & all(isfinite(xy_um), 2);
xy = xy_um(ok, :);
col = couleur(ok);
n = numel(col);
G = struct('n_colorees', n, 'n_couleur', accumarray(col, 1, [nC 1])', 'nshuffle', nshuffle, ...
    'rapport', NaN, 'rapport_nul', NaN, 'indice', NaN, 'z', NaN, 'p', NaN, ...
    'rapport_c', nan(1, nC), 'rapport_nul_c', nan(1, nC), 'indice_c', nan(1, nC), ...
    'z_c', nan(1, nC), 'p_c', nan(1, nC), 'rayon_um', rayon_um, 'n_paires_rayon', NaN, ...
    'obs_attendu', NaN, 'z_rayon', NaN, 'p_rayon', NaN, ...
    'obs_attendu_c', nan(1, nC), 'z_rayon_c', nan(1, nC), 'p_rayon_c', nan(1, nC));
if n < 4 || sum(G.n_couleur >= 2) < 1 || sum(G.n_couleur >= 1) < 2
    return
end

D = squareform(pdist(xy));
dsum = sum(D, 2);
[r_obs, rc_obs] = rapports(D, dsum, col, nC);
[pi_, pj] = find(triu(D <= rayon_um, 1));          % paires à moins de rayon_um
G.n_paires_rayon = numel(pi_);
[k_obs, kc_obs] = comptes(pi_, pj, col, nC);
rng(graine);
r_nul  = nan(nshuffle, 1);
rc_nul = nan(nshuffle, nC);
k_nul  = nan(nshuffle, 1);
kc_nul = nan(nshuffle, nC);
for b = 1:nshuffle
    cb = col(randperm(n));
    [r_nul(b), rc_nul(b, :)] = rapports(D, dsum, cb, nC);
    [k_nul(b), kc_nul(b, :)] = comptes(pi_, pj, cb, nC);
end
if G.n_paires_rayon > 0
    G.obs_attendu = k_obs / mean(k_nul);
    G.z_rayon     = (k_obs - mean(k_nul)) / ecart_type(k_nul);
    G.p_rayon     = (1 + sum(k_nul >= k_obs)) / (nshuffle + 1);
    att_c = mean(kc_nul, 1);
    okc = G.n_couleur >= 2 & att_c > 0;
    G.obs_attendu_c(okc) = kc_obs(okc) ./ att_c(okc);
    sd_c = ecart_type(kc_nul);
    G.z_rayon_c(okc) = (kc_obs(okc) - att_c(okc)) ./ sd_c(okc);
    G.p_rayon_c(okc) = (1 + sum(kc_nul(:, okc) >= kc_obs(okc), 1)) / (nshuffle + 1);
end

G.rapport     = r_obs;
G.rapport_nul = mean(r_nul, 'omitnan');
G.indice      = r_obs / G.rapport_nul;
G.z           = (r_obs - G.rapport_nul) / ecart_type(r_nul);
G.p           = (1 + sum(r_nul <= r_obs * (1 + 1e-12))) / (nshuffle + 1);

G.rapport_c     = rc_obs;
G.rapport_nul_c = mean(rc_nul, 1, 'omitnan');
G.indice_c      = rc_obs ./ G.rapport_nul_c;
G.z_c           = (rc_obs - G.rapport_nul_c) ./ ecart_type(rc_nul);
G.p_c           = (1 + sum(rc_nul <= rc_obs * (1 + 1e-12), 1)) / (nshuffle + 1);
sans = isnan(rc_obs);
G.indice_c(sans) = NaN;
G.z_c(sans)      = NaN;
G.p_c(sans)      = NaN;
end


function [r, rc] = rapports(D, dsum, col, nC)
% Sommes des distances : à l'intérieur de chaque couleur, et de chaque
% couleur vers toutes les cellules (paires ordonnées : chaque paire compte 2 fois,
% ce qui s'annule dans les moyennes)
n = numel(col);
nc = accumarray(col, 1, [nC 1])';
s_ligne = accumarray(col, dsum, [nC 1])';
s_meme = zeros(1, nC);
for a = 1:nC
    if nc(a) >= 2
        idx = col == a;
        s_meme(a) = sum(sum(D(idx, idx)));
    end
end
paires_meme = nc .* (nc - 1);
n_meme = sum(paires_meme);
n_diff = n * (n - 1) - n_meme;
if n_meme > 0 && n_diff > 0
    r = (sum(s_meme) / n_meme) / ((sum(s_ligne) - sum(s_meme)) / n_diff);
else
    r = NaN;
end
rc = nan(1, nC);
for a = 1:nC
    if nc(a) >= 2 && nc(a) < n
        rc(a) = (s_meme(a) / paires_meme(a)) / ((s_ligne(a) - s_meme(a)) / (nc(a) * (n - nc(a))));
    end
end
end


function [k, kc] = comptes(pi_, pj, col, nC)
% paires de même couleur parmi les paires proches : toutes couleurs, et par couleur
meme = col(pi_) == col(pj);
kc = accumarray(col(pi_(meme)), 1, [nC 1])';
k = sum(kc);
end


function s = ecart_type(X)
% écart-type de chaque colonne, sans les NaN
s = nan(1, size(X, 2));
for j = 1:size(X, 2)
    x = X(~isnan(X(:, j)), j);
    if numel(x) > 1, s(j) = std(x); end
end
end
