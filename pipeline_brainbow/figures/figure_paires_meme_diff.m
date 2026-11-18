function fig = figure_paires_meme_diff(dossier, rayon_um, mesures)
%FIGURE_PAIRES_MEME_DIFF  Même couleur contre couleurs différentes, une
%   colonne par mesure (STTC et Pearson du dF/F par défaut) ; en haut avec
%   les résidus, en bas la mesure brute (sans aucune correction),
%   sans rien recalculer : lit paires_agregation.mat (agregation_paires).
%
%   figure_paires_meme_diff('E:\...\paires_agregation')             % rayon principal
%   figure_paires_meme_diff('E:\...\paires_agregation', 100)        % paires <= 100 µm
%   figure_paires_meme_diff('E:\...\paires_agregation', 100, {'sttc', 'dff', 'sttc_500ms'})
%
%   Permutation principale de config.paires (en général 'couplage').
%   Écrit paires_meme_diff_<rayon>um.png.

dossier = regexprep(char(dossier), '[\\/]+$', '');
S = load(fullfile(dossier, 'paires_agregation.mat'));
if nargin < 2 || isempty(rayon_um), rayon_um = S.principal_dist(2); end
if nargin < 3 || isempty(mesures), mesures = {'sttc', 'dff'}; end
mesures = cellstr(mesures);
if S.cumulatif
    pd = [0, rayon_um];
else
    pd = S.principal_dist;
    warning('Distances en tranches : rayon ignoré, tranche %g-%g µm.', pd(1), pd(2));
end
perm0 = S.opts.permutations{1};
nM = numel(mesures);
fig = figure('Color', 'w', 'Position', [60 40 560 * nM 980]);
tl = tiledlayout(fig, 2, nM, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Same vs different colour at rest (%s permutation; %d animals, %d sessions)', ...
    perm0, S.n_animaux, S.n_sessions), 'FontWeight', 'bold', 'FontSize', 11);
% 1re ligne : résidus (distance et niveau de chaque cellule retirés) ; 2e : brut
lettres = 'ABCDEFGH';
for k = 1:nM
    ax = nexttile(tl, k);
    panneau_paires_meme_diff(ax, S, mesures{k}, perm0, pd, [lettres(k) '1. '], 'residu');
    ax = nexttile(tl, nM + k);
    panneau_paires_meme_diff(ax, S, mesures{k}, perm0, pd, [lettres(k) '2. '], 'brut');
end
exportgraphics(fig, fullfile(dossier, sprintf('paires_meme_diff_%gum.png', rayon_um)), 'Resolution', 200);
end
