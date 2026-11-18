function T = diagnostic_seuil_SCE(dossier_session, synchronous_frames, nNull)
%DIAGNOSTIC_SEUIL_SCE  Aide au choix du seuil des SCE, sur UNE session déjà prétraitée.
%
%   diagnostic_seuil_SCE(dossier_session)
%   diagnostic_seuil_SCE(dossier_session, synchronous_frames)
%
%   Lit pretraitement.mat (aucune couleur, aucun recalcul du prétraitement)
%   et compare, pour chaque seuil s, le nombre de SCE détectés (pics de MAct
%   > s) au nombre attendu par hasard sous deux modèles nuls, qui décalent
%   le train de chaque cellule d'un bloc, à l'intérieur du repos :
%     - +/- 10 trames : les shuffles de la règle 'percentile' actuelle ; ils
%                       gardent les fluctuations lentes communes, mais aussi
%                       une partie des vrais SCE (étalés sur 21 trames)
%     - circulaire    : cellules indépendantes ; sous-estime le hasard s'il
%                       y a des fluctuations lentes communes pendant le repos
%   synchronous_frames (facultatif) : fenêtre de co-activité à tester à la
%   place de celle du prétraitement (MAct = cellules actives entre les trames
%   i et i + synchronous_frames). Une fenêtre plus courte réduit les
%   coïncidences dues au hasard, si les vrais SCE sont assez synchrones.
%
%   Fichier autonome : n'a besoin d'aucune autre fonction du pipeline.
%   Peut tourner dans une deuxième fenêtre MATLAB pendant l'analyse.
%   Écrit diagnostic_seuil_SCE_fenetre<n>.png dans le dossier de la session.

if nargin < 3, nNull = 20; end
etat_rng = rng;                                   % générateur remis en l'état à la fin
remettre_rng = onCleanup(@() rng(etat_rng)); %#ok<NASGU>
fichier_pre = fullfile(dossier_session, 'pretraitement.mat');
P = load(fichier_pre, 'Raster', 'opts', 'WinRest', 'sampling_rate');
opts = P.opts;
fs = P.sampling_rate;
fenetre_pretraitement = opts.synchronous_frames;
if nargin >= 2 && ~isempty(synchronous_frames)
    opts.synchronous_frames = synchronous_frames;
end
Nz = size(P.Raster, 2);
repos = false(1, Nz);
repos(P.WinRest) = true;
opts_courbe = opts;
opts_courbe.sce_n_cells_threshold = min(opts.sce_n_cells_threshold, 5);   % courbes depuis 5 cellules

modeles = {'+/- 10 trames (règle actuelle)', 10; ...
           'circulaire (cellules indépendantes)', Inf};
fprintf('Calcul des modèles nuls (%d surrogates chacun) ...\n', nNull);
N_hasard = cell(1, size(modeles, 1));
for m = 1:size(modeles, 1)
    [seuils, N_obs, N_hasard{m}] = courbe_SCE_hasard(P.Raster, repos, opts_courbe, modeles{m, 2}, nNull, 1000 + m);
end

seuil_actuel = opts.sce_n_cells_threshold_calcule;
fprintf('\nSession : %s\n', dossier_session);
fprintf('%d cellules, %.0f s de repos, %.1f Hz : +/- 10 trames = +/- %.2f s\n', ...
    size(P.Raster, 1), sum(repos) / fs, fs, 10 / fs);
fprintf('Fenêtre de co-activité testée : %d trames = %.0f ms (prétraitement : %d trames)\n', ...
    opts.synchronous_frames + 1, 1000 * (opts.synchronous_frames + 1) / fs, fenetre_pretraitement + 1);
fprintf('Seuil actuel (prétraitement) : %d cellules (règle %s)\n', seuil_actuel, regle(opts));

% Seuil que donnerait la règle 'fdr' pour chaque modèle
fprintf('\nPlus petit seuil où hasard / détectés reste sous la limite :\n');
for m = 1:size(modeles, 1)
    s5 = seuil_fdr(seuils, N_obs, N_hasard{m}, 0.05);
    s10 = seuil_fdr(seuils, N_obs, N_hasard{m}, 0.10);
    fprintf('   %-36s  5 %% : %3d cellules (%4d SCE)   10 %% : %3d cellules (%4d SCE)\n', modeles{m, 1}, ...
        s5, n_au_dessus(seuils, N_obs, s5), s10, n_au_dessus(seuils, N_obs, s10));
end

% Tableau
s_affiches = unique([seuils(1), 10:5:60, seuil_actuel]);
s_affiches = s_affiches(s_affiches >= seuils(1) & s_affiches <= seuils(end));
i_aff = arrayfun(@(s) find(seuils == s, 1), s_affiches);
T = table(seuils(i_aff)', N_obs(i_aff)', N_hasard{1}(i_aff)', N_hasard{2}(i_aff)', ...
    'VariableNames', {'seuil', 'SCE_detectes', 'hasard_10trames', 'hasard_circulaire'});
fprintf('\nSCE détectés au-dessus de chaque seuil, et attendus par hasard :\n');
disp(T)

% Figure
fig = figure('Visible', 'off', 'Position', [100 100 1500 420]);
couleurs = [0.85 0.33 0.1; 0.3 0.6 0.3];
subplot(1, 3, 1); hold on
semilogy(seuils, max(N_obs, 0.1), 'k-', 'LineWidth', 2);
for m = 1:size(modeles, 1)
    semilogy(seuils, max(N_hasard{m}, 0.1), '-', 'Color', couleurs(m, :), 'LineWidth', 1.5);
end
set(gca, 'YScale', 'log');
xline(seuil_actuel, '--k');
xlabel('Seuil (cellules co-actives)'); ylabel('Nombre de SCE au-dessus du seuil');
legend([{'données'}, modeles(:, 1)', {'seuil actuel'}], 'Location', 'northeast');
title(sprintf('Fenêtre de %d trames (%.0f ms)', opts.synchronous_frames + 1, ...
    1000 * (opts.synchronous_frames + 1) / fs));

subplot(1, 3, 2); hold on
for m = 1:size(modeles, 1)
    part = N_hasard{m} ./ max(N_obs, 1);
    plot(seuils, 100 * min(part, 1), '-', 'Color', couleurs(m, :), 'LineWidth', 1.5);
end
yline(5, ':k'); yline(10, ':k'); xline(seuil_actuel, '--k');
ylim([0 100]);
xlabel('Seuil (cellules co-actives)'); ylabel('Part due au hasard (%)');
legend(modeles(:, 1)', 'Location', 'northeast');
title('Hasard / détectés');

subplot(1, 3, 3);
taille_SCE = -diff([N_obs 0]);              % nombre de SCE de taille exactement s+1
bar(seuils + 1, taille_SCE, 1, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
xline(seuil_actuel + 0.5, '--k');
xlabel('Taille du SCE (cellules)'); ylabel('Nombre de SCE');
title('Taille des pics de co-activité');
nom_png = sprintf('diagnostic_seuil_SCE_fenetre%d.png', opts.synchronous_frames + 1);
exportgraphics(fig, fullfile(dossier_session, nom_png), 'Resolution', 150);
close(fig)
fprintf('Figure : %s\n', fullfile(dossier_session, nom_png));
end


function s = seuil_fdr(seuils, N_obs, N_hasard, fdr)
trop = N_hasard > fdr * N_obs;
dernier = find(trop, 1, 'last');
if isempty(dernier), s = seuils(1); else, s = seuils(dernier) + 1; end
end

function n = n_au_dessus(seuils, N_obs, s)
i = find(seuils == s, 1);
if isempty(i), n = 0; else, n = N_obs(i); end
end

function r = regle(opts)
if isfield(opts, 'sce_regle'), r = opts.sce_regle; else, r = 'percentile'; end
end


%% =========================================================================
function [seuils, N_obs, N_hasard] = courbe_SCE_hasard(Raster, repos, opts, decalage_max, nNull, graine)
%COURBE_SCE_HASARD  SCE détectés au-dessus de chaque seuil : données et hasard.
%
%   seuils   : seuil minimum (opts.sce_n_cells_threshold) jusqu'au plus grand pic
%   N_obs    : nombre de SCE détectés dans les données au-dessus de chaque seuil
%   N_hasard : nombre moyen de SCE détectés dans nNull surrogates où le train
%              de chaque cellule est décalé d'un bloc de k trames, k tiré entre
%              -decalage_max et +decalage_max (Inf = décalage circulaire
%              complet). L'activité de chaque cellule est gardée, la synchronie
%              plus fine que le décalage est détruite.
%
%   « Au-dessus du seuil s » = pic de MAct STRICTEMENT supérieur à s, comme
%   findpeaks(..., 'MinPeakHeight', s) dans SCE_detection5.
%
%   Les décalages se font À L'INTÉRIEUR des trames de repos mises bout à bout
%   (repos = masque logique) : les transients ne sont détectés qu'au repos.
%   Décaler sur toute la session étalerait l'activité sur les périodes de
%   course, où le raster est vide, et sous-estimerait le hasard.

synchronous_frames = opts.synchronous_frames;
MinPeakDistancesce = opts.MinPeakDistancesce;
seuil_min          = opts.sce_n_cells_threshold;
if isfield(opts, 'max_cells_SCE'), max_cells = opts.max_cells_SCE; else, max_cells = 150; end

R = Raster(:, logical(repos));          % trames de repos mises bout à bout
[NCell, Nr] = size(R);
pics_obs = pics_SCE(activite_cosynchrone(R, synchronous_frames), seuil_min, MinPeakDistancesce, max_cells);
seuils = seuil_min:max([seuil_min; pics_obs(:)]);
N_obs = arrayfun(@(s) sum(pics_obs > s), seuils);

N_hasard = zeros(size(seuils));
if Nr <= 2 * synchronous_frames
    return
end
for n = 1:nNull
    fixer_flux(graine, n, nNull);
    Rn = zeros(NCell, Nr);
    for c = 1:NCell
        if isinf(decalage_max)
            k = randi(Nr);
        else
            k = randi([-decalage_max decalage_max]);
        end
        Rn(c, :) = circshift(R(c, :), k);
    end
    p = pics_SCE(activite_cosynchrone(Rn, synchronous_frames), seuil_min, MinPeakDistancesce, max_cells);
    N_hasard = N_hasard + arrayfun(@(s) sum(p > s), seuils) / nNull;
end
end


function pks = pics_SCE(MAct, seuil_min, MinPeakDistancesce, max_cells)
% Pics de MAct au-dessus de seuil_min, avec la même règle que SCE_detection5.
% Avec MinPeakDistance, findpeaks garde les pics du plus haut au plus bas :
% les pics au-dessus d'un seuil s plus élevé sont donc exactement ceux qu'on
% obtiendrait en relançant findpeaks avec s.
if numel(MAct) <= MinPeakDistancesce + 1 || max(MAct) <= seuil_min
    pks = zeros(1, 0);
    return
end
pks = findpeaks(double(MAct), 'MinPeakHeight', seuil_min, 'MinPeakDistance', MinPeakDistancesce);
pks = pks(pks <= max_cells);
end


function MAct = activite_cosynchrone(Raster, synchronous_frames)
%ACTIVITE_COSYNCHRONE  MAct(i) = nombre de cellules actives au moins une fois
%   entre les trames i et i + synchronous_frames (même calcul que la boucle
%   de SCE_detection4, en beaucoup plus rapide).
Nz = size(Raster, 2);
if Nz <= synchronous_frames
    MAct = zeros(1, 0);
    return
end
fenetre_max = movmax(Raster, [0 synchronous_frames], 2);
MAct = sum(fenetre_max(:, 1:Nz - synchronous_frames), 1);
end


function fixer_flux(graine, indice, nombre)
% Flux aléatoire propre à (graine, indice) : résultats reproductibles.
flux = RandStream.create('mrg32k3a', 'NumStreams', max(nombre, indice), ...
    'StreamIndices', indice, 'Seed', graine);
RandStream.setGlobalStream(flux);
end
