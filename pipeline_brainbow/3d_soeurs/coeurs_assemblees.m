function C = coeurs_assemblees(session, config, K)
%COEURS_ASSEMBLEES  Étape 3d, UNE session : chaque assemblée a-t-elle un
%   cœur clonal ? (7 nov., demandé par l'utilisateur : pouvoir dire oui ou
%   non pour chaque assemblée, en plus du test général)
%   Score de l'assemblée = somme, sur les paires de ses membres colorés, de
%   pi (probabilité d'être sœurs, carte_soeurs : même couleur ET proches ;
%   0 sinon) = nombre attendu de paires de sœurs dans l'assemblée. Un petit
%   cœur de sœurs proches compte même dans une grosse assemblée ; les
%   couleurs rares pèsent plus ; les paires de même couleur lointaines ne
%   comptent pas.
%   Deux p par assemblée (unilatéraux, opts.nPerm_coeur tirages) :
%     p1 « contient des sœurs » : couleurs mélangées entre les cellules
%        colorées du champ (mêmes membres, mêmes positions). Simulation : une
%        assemblée compacte SANS aucun choix clonal passe ce test une fois
%        sur deux (elle contient les sœurs de son voisinage) ; descriptif ;
%     p2 « cœur clonal » (le verdict) : fausses assemblées où chaque
%        membre coloré est remplacé par une de ses opts.k_voisins cellules
%        colorées les plus proches (couleurs réelles) : l'assemblée
%        recrute-t-elle des sœurs plus que sa forme dans l'espace ne
%        l'explique ? (les clones étant compacts, toute assemblée compacte
%        contient quelques sœurs)
%   Testable : au moins opts.min_colorees membres colorés et au moins une
%   paire de membres colorés à une distance où pi peut être > 0.
%   Le verdict (FDR sur toutes les assemblées) est dans agregation_coeurs.
%   Écrit <session>\soeurs\coeurs.mat (C : une ligne par assemblée).

opts = options_soeurs(config);
fichier = fullfile(session.dossier, 'soeurs', 'coeurs.mat');
R = load(fullfile(session.dossier, 'soeurs', 'soeurs.mat'), 'col', 'xy_um');
col = R.col(:);
xy = R.xy_um;
NCell = numel(col);
[asm, best_K] = assemblees_retenues(session, NCell);
C = table();
if isempty(asm)
    save(fichier, 'C');
    return
end
path = avec_separateur(session.dossier_suite2p);
[~, ~, ~, ~, ~, Coomasqx, Coomasqy] = charger_geometrie(path);
chev = masques_chevauchent(Coomasqx, Coomasqy);
colorees = find(col >= 1 & col <= 6);
ib = find(any(K.pi_cb > 0, 2), 1, 'last');       % dernière tranche où pi peut être > 0
if isempty(ib), ib = 0; end
nA = numel(asm);

% paires de membres colorés de chaque assemblée
I = []; J = []; Bn = []; A = [];
membres_col = cell(1, nA);
n_proches = zeros(1, nA);
for a = 1:nA
    m = intersect(asm{a}, colorees);
    membres_col{a} = m;
    [i, j, b] = paires(m, xy, chev, opts.bornes_um);
    I = [I; i]; J = [J; j]; Bn = [Bn; b]; A = [A; a * ones(numel(i), 1)]; %#ok<AGROW>
    n_proches(a) = sum(b <= ib);
end
S_obs = score(col, I, J, Bn, A, nA, K.pi_cb);

% p1 : couleurs mélangées dans le champ
fixer_flux_aleatoire(session.graine, 11, 1);
nul1 = zeros(opts.nPerm_coeur, nA);
for p = 1:opts.nPerm_coeur
    cp = col;
    cp(colorees) = col(colorees(randperm(numel(colorees))));
    nul1(p, :) = score(cp, I, J, Bn, A, nA, K.pi_cb);
end

% p2 : fausses assemblées de voisines (couleurs réelles)
xc = xy(colorees, :);
Dc = sqrt((xc(:, 1) - xc(:, 1)') .^ 2 + (xc(:, 2) - xc(:, 2)') .^ 2);
Dc(1:numel(colorees) + 1:end) = Inf;
[~, ordre] = sort(Dc, 2);
kv = min(opts.k_voisins, numel(colorees) - 1);
voisins = ordre(:, 1:kv);                          % indices dans colorees
[~, rang] = ismember(1:NCell, colorees);           % rang de chaque cellule dans colorees (0 sinon)
fixer_flux_aleatoire(session.graine, 12, 1);
nul2 = zeros(opts.nPerm_coeur, nA);
for a = 1:nA
    m = membres_col{a};
    if numel(m) < 2, continue; end
    r = rang(m);
    for p = 1:opts.nPerm_coeur
        faux = zeros(1, numel(m));
        pris = false(1, numel(colorees));
        for u = randperm(numel(m))
            cand = voisins(r(u), randperm(kv));
            libre = cand(~pris(cand));
            if isempty(libre), libre = cand; end
            faux(u) = libre(1);
            pris(libre(1)) = true;
        end
        cel = colorees(faux);
        [i, j, b] = paires(cel, xy, chev, opts.bornes_um);
        nul2(p, a) = score(col, i, j, b, ones(numel(i), 1), 1, K.pi_cb);
    end
end

p1 = (1 + sum(nul1 >= S_obs - 1e-12, 1)) / (opts.nPerm_coeur + 1);
p2 = (1 + sum(nul2 >= S_obs - 1e-12, 1)) / (opts.nPerm_coeur + 1);
z1 = (S_obs - mean(nul1, 1)) ./ std(nul1, 0, 1);
z2 = (S_obs - mean(nul2, 1)) ./ std(nul2, 0, 1);
n_col = cellfun(@numel, membres_col);
testable = n_col >= opts.min_colorees & n_proches > 0;

% cellules du cœur : membres colorés dans une paire de même couleur avec pi > 0
coeur = cell(1, nA);
coul_coeur = cell(1, nA);
ci = col(I);
mm = ci == col(J);
v = zeros(size(I));
v(mm) = K.pi_cb(sub2ind(size(K.pi_cb), Bn(mm), ci(mm)));
noms_c = {'R', 'G', 'B', 'RG', 'RB', 'GB'};
for a = 1:nA
    k = A == a & v > 0;
    cel = unique([I(k); J(k)])';
    coeur{a} = strjoin(arrayfun(@num2str, cel, 'UniformOutput', false), ' ');
    cc = unique(col(cel));
    coul_coeur{a} = strjoin(noms_c(cc(cc >= 1 & cc <= 6)), ' ');
end
C = table(repmat({session.id_animal}, nA, 1), repmat({session.nom}, nA, 1), repmat(best_K, nA, 1), ...
    (1:nA)', cellfun(@numel, asm)', n_col', n_proches', S_obs', mean(nul1, 1)', z1', p1', ...
    mean(nul2, 1)', z2', p2', testable', coeur', coul_coeur', 'VariableNames', {'animal', 'session', 'K', ...
    'assemblee', 'n_membres', 'n_colorees', 'n_paires_proches', 'score', 'attendu_champ', 'z_champ', ...
    'p_champ', 'attendu_voisines', 'z_voisines', 'p_voisines', 'testable', 'cellules_coeur', ...
    'couleurs_coeur'});
save(fichier, 'C');
end


function [i, j, b] = paires(cel, xy, chev, bornes)
% paires de cellules (ROI qui ne se chevauchent pas) et leur tranche de distance
cel = cel(:);
n = numel(cel);
if n < 2
    i = zeros(0, 1); j = zeros(0, 1); b = zeros(0, 1);
    return
end
[p, q] = find(triu(true(n), 1));
i = cel(p);
j = cel(q);
ok = ~chev(sub2ind(size(chev), i, j));
i = i(ok);
j = j(ok);
d = sqrt(sum((xy(i, :) - xy(j, :)) .^ 2, 2));
b = discretize(d, bornes);
end


function S = score(col, I, J, B, A, nA, pi_cb)
% somme des pi des paires de même couleur, par assemblée
S = zeros(1, nA);
if isempty(I), return; end
ci = col(I);
m = ci == col(J) & ci >= 1 & ci <= 6;
v = zeros(size(I));
v(m) = pi_cb(sub2ind(size(pi_cb), B(m), ci(m)));
S = accumarray(A, v, [nA 1])';
end
