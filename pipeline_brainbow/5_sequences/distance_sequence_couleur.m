function C = distance_sequence_couleur(ordre, colorcell, xy_um, nPerm, min_paires, graine)
%DISTANCE_SEQUENCE_COULEUR  Les cellules de même couleur sont-elles plus
%   proches dans la séquence ? (test de l'utilisateur, 6 nov. : « la distance
%   dans la séquence, même couleur contre couleur différente »)
%   ordre     : cellules de la séquence, dans leur ordre (ordre_tout)
%   colorcell : couleur de chaque cellule (1 à 6 = colorée)
%   xy_um     : centre de chaque cellule (µm), pour le contrôle spatial
%   Position dans la séquence : rang / (nombre de cellules de la séquence - 1),
%   de 0 (première) à 1 (dernière). Distance d'une paire = écart de position
%   (fraction de la séquence), comme la « distance intra-assemblée » du code
%   d'origine (rang / NCell).
%   Paires : cellules COLORÉES de la séquence ; même couleur contre couleurs
%   différentes, sans aucune correction : delta = moyenne(même) -
%   moyenne(différentes) ; < 0 = même couleur plus proche.
%   p : couleurs mélangées entre les cellules colorées de la séquence (nPerm
%   fois), unilatéral (delta mélangé <= delta observé).
%   rho_espace (contrôle, pas une correction) : Spearman entre la distance
%   dans la séquence et la distance dans l'espace, sur les mêmes paires.
%   Si les voisines dans l'espace étaient voisines dans la séquence, un effet
%   de couleur pourrait venir de l'espace (les couleurs sont regroupées).

C = struct('n_sequence', numel(ordre), 'n_colorees', 0, 'n_meme', 0, 'n_diff', 0, 'd_meme', NaN, ...
    'd_diff', NaN, 'delta', NaN, 'pct', NaN, 'p', NaN, 'rho_espace', NaN, 'p_espace', NaN);
n = numel(ordre);
if n < 3, return; end
position = zeros(max(ordre), 1);
position(ordre) = (0:n - 1)' / (n - 1);
col = colorcell(ordre);
col = col(:);
ok = col >= 1 & col <= 6;
cel = ordre(ok);
col = col(ok);
C.n_colorees = numel(cel);
if numel(cel) < 3, return; end
[i, j] = find(triu(true(numel(cel)), 1));
d_seq = abs(position(cel(i)) - position(cel(j)));
meme = col(i) == col(j);
C.n_meme = sum(meme);
C.n_diff = sum(~meme);
if C.n_meme < min_paires || C.n_diff < min_paires, return; end
C.d_meme = mean(d_seq(meme));
C.d_diff = mean(d_seq(~meme));
C.delta = C.d_meme - C.d_diff;
C.pct = 100 * C.delta / C.d_diff;

fixer_flux_aleatoire(graine, 6, 1);
delta_perm = zeros(nPerm, 1);
for k = 1:nPerm
    c2 = col(randperm(numel(col)));
    m2 = c2(i) == c2(j);
    delta_perm(k) = mean(d_seq(m2)) - mean(d_seq(~m2));
end
C.p = (1 + sum(delta_perm <= C.delta + 1e-12)) / (nPerm + 1);

% contrôle : la séquence suit-elle l'espace ?
xy = xy_um(cel, :);
d_esp = sqrt(sum((xy(i, :) - xy(j, :)) .^ 2, 2));
[C.rho_espace, C.p_espace] = corr(d_esp, d_seq, 'Type', 'Spearman');
end
