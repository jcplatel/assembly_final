function Tses = ressemblance_SCE(cfg, dossier_analyse, dossier_sortie, nNull, filtre_agregation)
%   filtre_agregation (défaut true) : true = assemblées >= 15 cellules et
%   >= 5 SCE (test de couleur) ; false = toutes les assemblées du clustering
%RESSEMBLANCE_SCE  (18 nov., demandée par l'utilisateur) Les SCE « sans
%   assemblée » se ressemblent-ils entre eux plus que par hasard ? Si oui :
%   assemblées manquées par le clustering (hypothèse B). Sinon : plutôt du
%   fond (hypothèse A). AUCUNE couleur.
%
%   Assemblées : K retenu à l'étape 2 (dossier_analyse\best_k_par_session.mat),
%   assemblees_v2.mat, mêmes critères que l'agrégation (>= 15 cellules, >= 5 SCE).
%   Ressemblance de deux SCE : cosinus des vecteurs de cellules actives
%   = cellules communes / sqrt(n1 x n2).
%   Hasard : nNull mélanges curveball des SCE du groupe (chaque cellule garde
%   son nombre de SCE, chaque SCE son nombre de cellules).
%   « Jumeau » d'un SCE : un autre SCE du groupe plus ressemblant qu'un seuil
%   choisi pour qu'au hasard un SCE ait un jumeau dans ~5 % des cas
%   (percentile 1 - 0,05 / (n - 1) des paires du hasard).
%   pct_paires : % de paires au-dessus du 99e percentile du hasard (1 % attendu).
%   Fait pour les SCE sans assemblée et, en contrôle positif, avec assemblée.
%   Écrit ressemblance_SCE_par_session.csv et ressemblance_SCE.png.

if nargin < 4, nNull = 20; end
if nargin < 5, filtre_agregation = true; end
ecart_min_s = 2;   % paires de SCE à moins de 2 s ignorées (un même transient peut compter dans les deux)
if ~isfolder(dossier_sortie), mkdir(dossier_sortie); end
cfg = dossiers_analyse(cfg);
S = lire_table_sessions(cfg);
B = load(fullfile(dossier_analyse, 'best_k_par_session.mat'), 'best_k'); B = B.best_k;
ses = {};
for s = 1:height(S)
    nom = S.nom{s}; fprintf('%d/%d %s\n', s, height(S), nom);
    try
        P = load(fullfile(S.dossier_source{s}, 'pretraitement.mat'), 'Race', 'TRace', 'sampling_rate');
        Race = P.Race > 0; nS = size(Race, 2); tps = double(P.TRace(:)') / P.sampling_rate;
        nAsm = zeros(1, nS); K = NaN;
        i = find(strcmp(B.session, nom), 1);
        if ~isempty(i) && ~isnan(B.K(i))
            K = B.K(i);
            A = load(fullfile(S.dossier_source{s}, sprintf('k%d', K), 'assemblees_v2.mat'), 'CellCl_final', 'PCl_final');
            C = A.CellCl_final > 0; Sa = A.PCl_final > 0;
            if filtre_agregation, garde = sum(C, 2) >= 15 & sum(Sa, 2) >= 5; else, garde = true(size(C, 1), 1); end
            nAsm = sum(Sa(garde, :), 1);
        end
        r = struct('session', nom, 'animal', S.id_animal{s}, 'K', K, 'n_SCE', nS);
        groupes = {'sans', nAsm == 0; 'avec', nAsm > 0};
        rng(s);
        for g = 1:2
            R = Race(:, groupes{g, 2});
            [m, t, f1, f1n, f2, f2n, n] = mesurer(R, nNull, tps(groupes{g, 2}), ecart_min_s);
            r.(['n_' groupes{g, 1}]) = n;
            r.(['sim_' groupes{g, 1}]) = m;               % ressemblance moyenne / hasard
            r.(['pct_paires_' groupes{g, 1}]) = t;        % % de paires > 99e pct du hasard (1 % attendu)
            r.(['pct_jumeau_' groupes{g, 1}]) = f1;       % % de SCE avec >= 1 jumeau
            r.(['pct_jumeau_hasard_' groupes{g, 1}]) = f1n;
            r.(['pct_2jumeaux_' groupes{g, 1}]) = f2;     % % de SCE avec >= 2 jumeaux (petit groupe)
            r.(['pct_2jumeaux_hasard_' groupes{g, 1}]) = f2n;
        end
        ses{end+1} = r; %#ok<AGROW>
    catch err
        fprintf('  ERREUR : %s\n', err.message);
    end
end
Tses = struct2table([ses{:}]);
writetable(Tses, fullfile(dossier_sortie, 'ressemblance_SCE_par_session.csv'));

fig = figure('Color', 'w', 'Position', [50 50 1300 430]);
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl); hold(ax, 'on');
x = [Tses.sim_sans, Tses.sim_avec];
plot(ax, [1 2], x', '-', 'Color', [0.75 0.75 0.75]);
scatter(ax, ones(height(Tses), 1), x(:, 1), 30, [0.3 0.45 0.8], 'filled');
scatter(ax, 2 * ones(height(Tses), 1), x(:, 2), 30, [0.85 0.4 0.2], 'filled');
yline(ax, 1, 'k--', 'hasard');
set(ax, 'XTick', [1 2], 'XTickLabel', {'SCE sans assemblée', 'SCE avec assemblée'}, 'XLim', [0.5 2.5]);
ylabel(ax, 'ressemblance moyenne / hasard'); title(ax, 'A. Ressemblance moyenne entre SCE du groupe');
ax = nexttile(tl); hold(ax, 'on');
scatter(ax, Tses.pct_jumeau_hasard_sans, Tses.pct_jumeau_sans, 35, [0.3 0.45 0.8], 'filled', 'DisplayName', 'sans assemblée');
scatter(ax, Tses.pct_jumeau_hasard_avec, Tses.pct_jumeau_avec, 35, [0.85 0.4 0.2], 'filled', 'DisplayName', 'avec assemblée');
lim = [0 100]; plot(ax, lim, lim, 'k--', 'HandleVisibility', 'off'); xlim(ax, lim); ylim(ax, lim);
xlabel(ax, '% de SCE avec un jumeau : hasard'); ylabel(ax, '% observé'); legend(ax, 'Location', 'southeast', 'Box', 'off');
title(ax, 'B. SCE qui ont au moins un jumeau');
ax = nexttile(tl); hold(ax, 'on');
scatter(ax, Tses.pct_2jumeaux_hasard_sans, Tses.pct_2jumeaux_sans, 35, [0.3 0.45 0.8], 'filled', 'DisplayName', 'sans assemblée');
scatter(ax, Tses.pct_2jumeaux_hasard_avec, Tses.pct_2jumeaux_avec, 35, [0.85 0.4 0.2], 'filled', 'DisplayName', 'avec assemblée');
plot(ax, lim, lim, 'k--', 'HandleVisibility', 'off'); xlim(ax, lim); ylim(ax, lim);
xlabel(ax, '% de SCE avec >= 2 jumeaux : hasard'); ylabel(ax, '% observé'); legend(ax, 'Location', 'southeast', 'Box', 'off');
title(ax, 'C. SCE dans un petit groupe (>= 2 jumeaux)');
title(tl, 'Les SCE sans assemblée se ressemblent-ils plus que par hasard ? (une session = un point ; au-dessus de la diagonale = plus que le hasard)', 'FontWeight', 'bold');
exportgraphics(fig, fullfile(dossier_sortie, 'ressemblance_SCE.png'), 'Resolution', 170);
end


function [m, t, f1, f1n, f2, f2n, n] = mesurer(R, nNull, tps, ecart_min_s)
n = size(R, 2);
m = NaN; t = NaN; f1 = NaN; f1n = NaN; f2 = NaN; f2n = NaN;
if n < 5, return; end
[Mo, haut] = cosinus(R);
loin = abs(tps(:) - tps(:)') >= ecart_min_s;     % paires assez éloignées dans le temps
haut = haut & loin; Mo(~loin) = 0;
nulls = cell(nNull, 1); Mn = cell(nNull, 1);
for k = 1:nNull
    Rn = curveball(R, 10 * sum(any(R, 2)));
    Mn{k} = cosinus(Rn); Mn{k}(~loin) = 0; nulls{k} = Mn{k}(haut);
end
tout = vertcat(nulls{:});
seuil = prctile(tout, 99);
m = mean(Mo(haut)) / mean(tout);
t = 100 * mean(Mo(haut) > seuil);
% jumeau : seuil tel qu'au hasard un SCE a un jumeau dans ~5 % des cas
seuil_j = prctile(tout, 100 * (1 - 0.05 / (n - 1)));
[f1, f2] = jumeaux(Mo, seuil_j);
a = zeros(nNull, 2);
for k = 1:nNull, [a(k, 1), a(k, 2)] = jumeaux(Mn{k}, seuil_j); end
f1n = mean(a(:, 1)); f2n = mean(a(:, 2));
end

function [M, haut] = cosinus(R)
R = double(R); n = sum(R, 1);
M = (R' * R) ./ sqrt(n' * n);
M(1:size(M, 1) + 1:end) = 0;
haut = triu(true(size(M)), 1);
end

function [f1, f2] = jumeaux(M, seuil)
nj = sum(M > seuil, 2);
f1 = 100 * mean(nj >= 1); f2 = 100 * mean(nj >= 2);
end
