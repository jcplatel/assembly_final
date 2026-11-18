function F = fca_arbre(ev, Tn, sigma_tr, nSurr, jusqu_racine, graine)
%FCA_ARBRE  Functional Clustering Algorithm (Feldt et al. 2009, Phys Rev E
%   80:056104 ; utilisé par Feldt Muldoon et al. 2013, PNAS 110:3567).
%   ev       : 1 x N cellules, trames des événements de chaque cellule sur
%              la ligne de temps du REPOS (1..Tn)
%   sigma_tr : écart-type du décalage des surrogates, en trames
%   nSurr    : nombre de surrogates par comparaison
%   jusqu_racine : après la coupe, continuer à fusionner (par significativité
%              décroissante) jusqu'à une seule unité, pour avoir l'arbre entier
%
%   Méthode (comme l'article) :
%   1. Similarité de deux trains a et b : AMD = moyenne de D(a->b) et
%      D(b->a), D(a->b) = moyenne, sur les événements de a, de la distance au
%      plus proche événement de b, divisée par Tn / (n_b + 1) (distance
%      attendue pour n_b événements au hasard). Petit = synchrone.
%   2. Surrogates : chaque événement décalé au hasard (loi normale, sigma),
%      ce qui garde le nombre d'événements et à peu près les intervalles.
%      Significativité « mise à l'échelle » = (médiane - observé) /
%      (médiane - 5e percentile) des AMD des surrogates : 1 = seuil de 95 %.
%   3. On fusionne la paire la plus significative : les deux trains sont
%      réunis en UN train (le groupe devient une unité), puis on recalcule.
%   4. Arrêt (coupe) quand plus aucune paire n'est significative (> 1).
%   Écarts à l'article, pour le temps de calcul (300 à 500 cellules ici,
%   environ 100 dans l'article) :
%   - nSurr surrogates (500 par défaut) au lieu de 5000 à 10 000 ;
%   - après une fusion, seule la ligne du nouveau groupe est recalculée (les
%     autres paires ne changent pas) ;
%   - les surrogates d'un groupe réutilisent ceux de ses membres (même loi) ;
%   - chaque sens de l'AMD décale le train de départ (l'autre reste fixe).
%   Temps en trames entières ; décalage réfléchi aux bords.
%
%   F.Z       : arbre au format linkage (N - 1 lignes : [noeud, noeud, étape])
%   F.sig     : significativité de chaque fusion
%   F.coupe   : nombre de fusions significatives (avant la coupe)
%   F.amas    : clusters FCA = groupes formés avant la coupe, non réunis
%               ensuite, d'au moins 3 cellules (indices de cellules)
%   F.noeuds_amas : numéro de ligne de Z de chaque cluster

N = numel(ev);
fixer_flux_aleatoire(graine, 8, 1);
DIST = zeros(N, Tn, 'single');
J = cell(1, N);
for u = 1:N
    e = sort(double(ev{u}(:)'));
    ev{u} = e;
    DIST(u, :) = distance_evenements(e, Tn);
    J{u} = decaler(e, nSurr, sigma_tr, Tn);
end
nb = cellfun(@numel, ev);
actif = true(1, N);
id = 1:N;                        % noeud de l'arbre porté par chaque place
SIG = -Inf(N, N);
for a = 1:N - 1
    autres = a + 1:N;
    s = ligne_sig(a, autres, ev, nb, DIST, J, Tn);
    SIG(a, autres) = s;
    SIG(autres, a) = s;
end
Z = zeros(N - 1, 3);
sig_fusion = nan(N - 1, 1);
coupe = N - 1;
coupe_trouvee = false;
for etape = 1:N - 1
    S2 = SIG;
    S2(~actif, :) = -Inf;
    S2(:, ~actif) = -Inf;
    S2(1:N + 1:end) = -Inf;
    [best, k] = max(S2(:));
    [a, b] = ind2sub([N N], k);
    if ~coupe_trouvee && ~(best > 1)
        coupe = etape - 1;
        coupe_trouvee = true;
        if ~jusqu_racine
            Z = Z(1:etape - 1, :);
            sig_fusion = sig_fusion(1:etape - 1);
            break
        end
    end
    if a > b, [a, b] = deal(b, a); end
    Z(etape, :) = [id(a), id(b), etape];
    sig_fusion(etape) = best;
    % le groupe prend la place de a ; b disparaît
    ev{a} = sort([ev{a}, ev{b}]);
    nb(a) = numel(ev{a});
    DIST(a, :) = min(DIST(a, :), DIST(b, :));
    J{a} = [J{a}, J{b}];
    J{b} = [];
    ev{b} = [];
    actif(b) = false;
    id(a) = N + etape;
    autres = find(actif);
    autres(autres == a) = [];
    if ~isempty(autres)
        s = ligne_sig(a, autres, ev, nb, DIST, J, Tn);
        SIG(a, autres) = s;
        SIG(autres, a) = s;
    end
end
F.Z = Z;
F.sig = sig_fusion;
F.coupe = coupe;
% clusters : noeuds formés avant la coupe et non réunis avant la coupe
[F.amas, F.noeuds_amas] = amas_avant_coupe(Z, N, coupe);
end


function s = ligne_sig(a, autres, ev, nb, DIST, J, Tn)
% Significativité de l'unité a contre chaque unité de autres
nS = size(J{a}, 1);
no = numel(autres);
na = nb(a);
w_ao = (nb(autres)' + 1) / Tn;            % normalisation D(a -> o) : / (Tn / (n_o + 1))
w_oa = (na + 1) / Tn;                     % normalisation D(o -> a)
% observé
d_ao = mean(double(DIST(autres, ev{a})), 2);                       % no x 1
tous = [ev{autres}];
prop = repelem(1:no, nb(autres));
d_oa = accumarray(prop(:), double(DIST(a, tous))', [no 1]) ./ nb(autres)';
amd = (d_ao .* w_ao + d_oa * w_oa) / 2;
% surrogates, sens a -> o : a décalé, o fixe (par paquets de surrogates)
f_ao = zeros(no, nS);
pas = max(1, floor(5e6 / max(1, no * na)));   % mémoire : <= 5e6 valeurs à la fois
for s0 = 1:pas:nS
    ss = s0:min(nS, s0 + pas - 1);
    X = double(DIST(autres, reshape(J{a}(ss, :), 1, [])));          % no x (|ss| * na)
    f_ao(:, ss) = squeeze(mean(reshape(X, no, numel(ss), na), 3));
end
% sens o -> a : chaque o décalé, a fixe
Jt = [J{autres}];                                                   % nS x somme des n_o
P = sparse(1:numel(prop), prop, 1 ./ nb(autres(prop)), numel(prop), no);
f_oa = reshape(double(DIST(a, Jt(:)')), nS, []) * P;                % nS x no
amd_s = (f_ao .* w_ao + f_oa' * w_oa) / 2;                          % no x nS
med = median(amd_s, 2);
q05 = prctile(amd_s, 5, 2);
s = (med - amd) ./ (med - q05);
s(~(med - q05 > 0)) = 0;
s = s';
end


function J = decaler(e, nS, sigma, Tn)
% nS trains décalés : chaque événement + N(0, sigma), réfléchi aux bords
x = round(e + sigma * randn(nS, numel(e)));
x(x < 1) = 2 - x(x < 1);
x(x > Tn) = 2 * Tn - x(x > Tn);
x = min(max(x, 1), Tn);
J = uint32(x);
end


function d = distance_evenements(e, Tn)
% distance (trames) de chaque trame 1..Tn à l'événement le plus proche
e = unique(e);
t = 1:Tn;
k = discretize(t, [-Inf, e, Inf]);
avant = [-Inf, e];
apres = [e, Inf];
d = single(min(t - avant(k), apres(k) - t));
end


function [amas, noeuds] = amas_avant_coupe(Z, N, coupe)
amas = {};
noeuds = [];
if coupe < 1, return; end
membres = cell(1, coupe);
parent_avant = false(1, coupe);
for i = 1:coupe
    m = [];
    for e = 1:2
        c = Z(i, e);
        if c <= N
            m = [m, c]; %#ok<AGROW>
        else
            m = [m, membres{c - N}]; %#ok<AGROW>
            parent_avant(c - N) = true;     % ce noeud a été réuni avant la coupe
        end
    end
    membres{i} = m;
end
for i = 1:coupe
    if ~parent_avant(i) && numel(membres{i}) >= 3
        amas{end+1} = sort(membres{i}); %#ok<AGROW>
        noeuds(end+1) = i; %#ok<AGROW>
    end
end
end
