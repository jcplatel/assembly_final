function L = table_minp(Kcol, nmax, nTirages, seuil, nmin_couleur, graine, midp, tarone)
%TABLE_MINP  p_minp d'une assemblée en fonction de son nombre de cellules
%   colorées, pour un champ donné (13 nov., screening_reglages).
%
%   Même test que stats_clonales_champ3, variante 'champ' / 'globale' :
%   statistique T = plus petite p hypergéométrique des 6 couleurs ; hasard =
%   couleurs mélangées entre les cellules colorées du champ. Sous ce hasard,
%   les n cellules colorées d'une assemblée reçoivent n couleurs tirées sans
%   remise dans le champ : la loi de T ne dépend que de n. Elle est donc
%   calculée UNE fois par champ pour chaque n (nTirages tirages), au lieu de
%   refaire 10 000 mélanges pour chaque assemblée de chaque réglage.
%
%   Kcol : nombre de cellules de chaque couleur dans le champ (1 x 6)
%   L.PT(n, c, k+1) : P(au moins k cellules de couleur c parmi n)
%   L.Ts(:, n)      : T du hasard, triés
%   L.borne(n)      : p_minp < seuil  <=>  T * (1 + tol) < borne(n)
%   L.taux_nul(n)   : P(p_minp < seuil) sous le hasard (<= seuil)
%   L.p_min_att(n)  : plus petite p_minp possible (toutes les cellules de la
%                     couleur la plus favorable)
%   L.z_moy, L.z_sd : pour le score continu (comme score_minp)
%
%   23 nov. (demandé par l'utilisateur : « à partir de maintenant tout calculer
%   en mid-p Tarone ») : mêmes règles que stats_clonales_champ3
%     midp   : ex aequo comptés à moitié : p = (n_< + 0,5 (n_= + 1)) / (nTirages + 1) ;
%              taux_nul et p_min_att de la même façon (L.midp)
%     tarone : pour chaque n, une couleur qui ne PEUT pas donner p < seuil
%              (toutes les n cellules de cette couleur, ou toutes celles du
%              champ) n'entre pas dans le min-p (L.PT mis à 1 pour elle)
%   Absents ou false = ancien calcul.

if nargin < 7 || isempty(midp), midp = false; end
if nargin < 8 || isempty(tarone), tarone = false; end
Kcol = Kcol(:)';
nC = numel(Kcol);
N = sum(Kcol);
nmax = min(nmax, N);
L = struct('Kcol', Kcol, 'N', N, 'nmax', nmax, 'nTirages', nTirages, 'seuil', seuil, 'tol', 1e-10, ...
    'midp', midp, 'tarone', tarone);
L.PT = ones(max(nmax, 1), nC, max(nmax, 1) + 1);
L.Ts = ones(nTirages, max(nmax, 1));
L.borne = -inf(max(nmax, 1), 1);
L.taux_nul = zeros(max(nmax, 1), 1);
L.p_min_att = ones(max(nmax, 1), 1);
L.z_moy = zeros(max(nmax, 1), 1);
L.z_sd = zeros(max(nmax, 1), 1);
valide = Kcol >= nmin_couleur;
if nmax < 1 || ~any(valide), return; end

for n = 1:nmax
    for c = find(valide)
        L.PT(n, c, 1:n + 1) = hygecdf((0:n) - 1, N, Kcol(c), n, 'upper');
        L.PT(n, c, 1) = 1;
    end
    if tarone        % couleurs qui ne peuvent pas atteindre le seuil pour cette taille
        for c = find(valide)
            if L.PT(n, c, min(n, Kcol(c)) + 1) >= seuil, L.PT(n, c, :) = 1; end
        end
    end
end

% Hasard : les n premières cellules d'un tirage sans remise de nmax cellules
% sont un tirage sans remise de n cellules -> un tirage sert à tous les n
urne = repelem(1:nC, Kcol)';
rng(graine);
lin = repmat((1:nmax)', 1, nC) + repmat(((1:nC) - 1) * nmax, nmax, 1);
Tnul = ones(nTirages, nmax);
for d = 1:nTirages
    tir = urne(randperm(N, nmax));
    cnt = cumsum(repmat(tir, 1, nC) == repmat(1:nC, nmax, 1), 1);    % nmax x nC
    Tnul(d, :) = min(L.PT(lin + cnt * (nmax * nC)), [], 2)';
end

tol = L.tol;
bornes_p = [0.5, nTirages - 0.5] / nTirages;
for n = 1:nmax
    t = sort(Tnul(:, n));
    L.Ts(:, n) = t;
    [u, ~, rang] = unique(t);
    cnt = accumarray(rang(:), 1);
    cum = cumsum(cnt);
    % p_minp si le T observé vaut u(k) : (1 + nb de tirages <= u(k)) / (nTirages + 1)
    % (mid-p : (nb < u(k) + 0,5 (nb = u(k) + 1)) / (nTirages + 1))
    if midp
        p_u = (cum - cnt + 0.5 * (cnt + 1)) / (nTirages + 1);
    else
        p_u = (1 + cum) / (nTirages + 1);
    end
    k = find(p_u < seuil, 1, 'last');
    if isempty(k), k = 0; end
    if k < numel(u), L.borne(n) = u(k + 1); else, L.borne(n) = Inf; end
    % taux_nul et score : mêmes définitions que stats_clonales_champ3
    if midp
        L.taux_nul(n) = sum(cnt((cum - 0.5 * cnt) / nTirages < seuil)) / nTirages;
    else
        L.taux_nul(n) = sum(cnt(cum / nTirages < seuil)) / nTirages;
    end
    m_max = min(n, Kcol);
    Tmin = Inf;
    for c = find(valide)
        Tmin = min(Tmin, L.PT(n, c, m_max(c) + 1));
    end
    n_le = sum(t <= Tmin * (1 + tol));
    if midp
        n_lt = sum(t < Tmin * (1 - tol));
        L.p_min_att(n) = (n_lt + 0.5 * (n_le - n_lt + 1)) / (nTirages + 1);
    else
        L.p_min_att(n) = (1 + n_le) / (nTirages + 1);
    end
    p_mid = (cum - 0.5 * cnt) / nTirages;
    z = sqrt(2) * erfcinv(2 * min(max(p_mid, bornes_p(1)), bornes_p(2)));   % = -norminv(p)
    L.z_moy(n) = sum(cnt .* z) / nTirages;
    L.z_sd(n) = sqrt(sum(cnt .* (z - L.z_moy(n)).^2) / (nTirages - 1));
end
end
