function T = test_soeurs_session(session, config, K)
%TEST_SOEURS_SESSION  Étape 3d, 2e passage, UNE session : la similarité
%   d'activité (résidus) augmente-t-elle avec la probabilité d'être sœurs ?
%   Paires de cellules colorées (ROI qui ne se chevauchent pas, similarité
%   mesurée). Pour chaque paire : pi = probabilité d'être sœurs (carte_soeurs :
%   sa couleur et sa distance ; 0 si les couleurs diffèrent).
%   Modèle (moindres carrés, toutes les paires) :
%       résidu ~ 1 + pi + [même couleur R] + ... + [même couleur GB]
%   Les 6 indicateurs « même couleur c » absorbent tout effet de la couleur
%   qui ne dépend pas de la distance (par exemple optique) : la pente de pi
%   ne vient que de la comparaison des paires de même couleur proches
%   (probables sœurs) aux paires de même couleur plus lointaines (pas
%   sœurs). beta = différence de similarité entre une paire sûrement sœur
%   (pi = 1) et une paire non sœur. Aussi la pente sans les indicateurs
%   (beta_simple, pour comparaison).
%   p de la session : couleurs mélangées entre les cellules colorées du
%   champ (opts.nPerm fois), pi et indicateurs recalculés, unilatéral.
%   Dose-réponse (figure) : résidu moyen des paires de couleurs différentes,
%   de même couleur avec pi = 0, puis par catégorie de pi.
%   Par couleur (7 nov., demande de l'utilisateur : le STTC peut dépendre de
%   la couleur ; descriptif) : pour la couleur c, paires c / autre couleur,
%   puis paires c-c par catégorie de pi (dose.<mesure>.moyenne_c, n_c) ;
%   pente beta_c du même modèle limité aux paires qui contiennent une
%   cellule c : résidu ~ 1 + pi + [même couleur c] (table Tc).
%   Écrit les résultats dans <session>\soeurs\soeurs.mat (T, Tc, dose).

opts = options_soeurs(config);
f = fullfile(session.dossier, 'soeurs', 'soeurs.mat');
R = load(f);
col = R.col(:);
NCell = numel(col);
colorees = find(col >= 1 & col <= 6);
haut = triu(true(NCell), 1);
xy = R.xy_um;
lignes = cell(0, 13);
lignes_c = cell(0, 8);
noms_c = {'R', 'G', 'B', 'RG', 'RB', 'GB'};
dose = struct();
for k = 1:numel(opts.mesures)
    nom = opts.mesures{k};
    if ~isfield(R.mesures, nom) || isempty(R.mesures.(nom).E)
        continue
    end
    Mk = R.mesures.(nom);
    P = Mk.V & haut;
    P(:, setdiff(1:NCell, colorees)) = false;
    P(setdiff(1:NCell, colorees), :) = false;
    [i, j] = find(P);
    y = double(Mk.E(sub2ind([NCell NCell], i, j)));
    d = sqrt(sum((xy(i, :) - xy(j, :)) .^ 2, 2));
    b = discretize(d, opts.bornes_um);
    [beta, beta_simple, pi_p, meme] = regression(y, col(i), col(j), b, K.pi_cb);
    % hasard : couleurs mélangées entre cellules colorées
    fixer_flux_aleatoire(session.graine + 10 * k, 10, 1);
    nul = nan(opts.nPerm, 1);
    for p = 1:opts.nPerm
        cp = col;
        cp(colorees) = col(colorees(randperm(numel(colorees))));
        nul(p) = regression(y, cp(i), cp(j), b, K.pi_cb);
    end
    p_session = (1 + sum(nul >= beta)) / (opts.nPerm + 1);
    z = (beta - mean(nul, 'omitnan')) / std(nul, 'omitnan');
    effet_pct = 100 * beta / Mk.brut_100;
    if strcmp(nom, 'sce'), effet_pct = NaN; end      % phi : moyenne proche de 0, effet donné par beta
    n_proches = sum(meme & pi_p > 0);
    lignes(end+1, :) = {session.id_animal, session.nom, nom, beta, beta_simple, z, p_session, effet_pct, ...
        Mk.brut_100, numel(y), sum(meme), n_proches, sum(pi_p)}; %#ok<AGROW>
    % dose-réponse
    cat = zeros(size(y));
    cat(~meme) = 1;
    cat(meme & pi_p == 0) = 2;
    e = opts.categories_pi;
    for c = 1:numel(e) - 1
        cat(meme & pi_p > e(c) & pi_p <= e(c + 1)) = 2 + c;
    end
    nC = numel(e) + 1;
    dose.(nom) = struct('moyenne', accumarray(cat(cat > 0), y(cat > 0), [nC 1], @mean, NaN)', ...
        'n', accumarray(cat(cat > 0), 1, [nC 1])', 'pi_moyen', accumarray(cat(cat > 0), pi_p(cat > 0), ...
        [nC 1], @mean, NaN)');
    % par couleur (descriptif)
    ci = col(i);
    cj = col(j);
    moy_c = nan(6, nC);
    n_c = zeros(6, nC);
    for c = 1:6
        autre = xor(ci == c, cj == c);              % une cellule c, l'autre d'une autre couleur
        cc = meme & ci == c;
        catc = zeros(size(y));
        catc(autre) = 1;
        catc(cc) = cat(cc);
        u = catc > 0;
        if ~any(u), continue; end
        moy_c(c, :) = accumarray(catc(u), y(u), [nC 1], @mean, NaN)';
        n_c(c, :) = accumarray(catc(u), 1, [nC 1])';
        sel = autre | cc;
        beta_c = pente_couleur(y(sel), pi_p(sel), cc(sel));
        effet_c = 100 * beta_c / Mk.brut_100;
        if strcmp(nom, 'sce'), effet_c = NaN; end
        lignes_c(end+1, :) = {session.id_animal, session.nom, nom, noms_c{c}, beta_c, effet_c, sum(cc), ...
            sum(cc & pi_p > 0)}; %#ok<AGROW>
    end
    dose.(nom).moyenne_c = moy_c;
    dose.(nom).n_c = n_c;
end
T = cell2table(lignes, 'VariableNames', {'animal', 'session', 'mesure', 'beta', 'beta_simple', 'z', 'p', ...
    'effet_pct', 'brut_100', 'n_paires', 'n_meme', 'n_meme_pi', 'soeurs_attendues'});
Tc = cell2table(lignes_c, 'VariableNames', {'animal', 'session', 'mesure', 'couleur', 'beta', 'effet_pct', ...
    'n_meme', 'n_meme_pi'});
save(f, 'T', 'Tc', 'dose', '-append');
end


function beta = pente_couleur(y, pi_p, cc)
% résidu ~ 1 + pi + [même couleur c], paires qui contiennent une cellule c
beta = NaN;
if sum(pi_p > 0) < 3 || sum(cc & pi_p == 0) < 3, return; end
X = [ones(numel(y), 1), pi_p, double(cc)];
if rank(X) < 3, return; end
coef = X \ y;
beta = coef(2);
end


function [beta, beta_simple, pi_p, meme] = regression(y, ci, cj, b, pi_cb)
meme = ci == cj;
pi_p = zeros(size(y));
pi_p(meme) = pi_cb(sub2ind(size(pi_cb), b(meme), ci(meme)));
beta = NaN;
beta_simple = NaN;
if ~any(pi_p > 0), return; end
n = numel(y);
D = zeros(n, 6);
D(sub2ind([n 6], find(meme), ci(meme))) = 1;
D = D(:, any(D, 1));
X = [ones(n, 1), pi_p, D];
coef = X \ y;
beta = coef(2);
if nargout > 1
    c0 = [ones(n, 1), pi_p] \ y;
    beta_simple = c0(2);
end
end
