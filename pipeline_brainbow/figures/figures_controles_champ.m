function figures_controles_champ(dossier)
%FIGURES_CONTROLES_CHAMP  Figures des contrôles du champ (étape 4), sans
%   rien recalculer : lit controles_champ.mat (écrit par controles_champ).
%
%   figures_controles_champ('E:\...\agregation\controles_champ')
%
%   colour_spatial_clustering.png
%     A. chaque session (son propre champ de vue) : indice de regroupement
%        (paires de même couleur / paires de couleurs différentes, divisé par
%        le même rapport après mélange des couleurs) ; 1 = hasard, < 1 = même
%        couleur plus proche. Rond plein : p < 0,05. Trait noir : moyenne de l'animal.
%     B. tous les animaux : un point par animal, moyenne +/- sem ; p = Wilcoxon
%        sur les animaux.
%     C. chaque couleur : un point par animal ; p = Wilcoxon sur les animaux.
%     D. selon la distance : paires de même couleur observées / attendues si
%        les couleurs étaient mélangées dans le champ (calcul de l'étape 3b,
%        lu dans <analyse>\paires_agregation\paires_agregation.mat), moyenne
%        des animaux, toutes couleurs et chaque couleur.
%     E. par anneau (0-50, 50-100, 100-200, 200-300, 300-500 µm, sans
%        cumuler) : paires de même couleur observées / attendues ; trait
%        gris : chaque animal ; noir : moyenne des animaux +/- sem ; trait
%        rouge pointillé et titre : le modèle mixte de agregation_paires
%        (log(observé / attendu) ~ log2(distance), animal et session
%        aléatoires) : excès à 100 µm et changement quand la distance double.
%   detection_by_colour.png
%     pour chaque couleur, un point par animal : moyenne de la couleur /
%     moyenne des cellules colorées du champ (1 = pareil), pour le bruit,
%     F0, les transients au repos, la participation aux SCE et la part de
%     membres d'assemblées ; et les cellules hyperactives (observées /
%     attendues d'après la part de chaque couleur).

dossier = regexprep(char(dossier), '[\\/]+$', '');
S = load(fullfile(dossier, 'controles_champ.mat'));
esp = S.esp; det = S.det; noms_coul = S.noms_coul;
nC = 6;
cols = [0.86 0.16 0.16; 0.16 0.66 0.22; 0.20 0.38 0.92; ...
        0.93 0.78 0.10; 0.82 0.22 0.80; 0.10 0.74 0.80];
GRIS = [0.45 0.45 0.45];
% 12 nov. : A-C en comptage (paires de même couleur à moins de esp.rayon_um,
% observées / attendues ; > 1 = regroupées), comme D et E. Anciens fichiers
% (sans esp.mesure) : indice de distance (< 1 = regroupées).
comptage = isfield(esp, 'mesure') && strcmp(esp.mesure, 'comptage');
if comptage
    lab_A = sprintf('same-colour pairs < %g \\mum: observed / expected (> 1: grouped)', esp.rayon_um);
    lab_C = sprintf('pairs of this colour < %g \\mum: observed / expected', esp.rayon_um);
else
    lab_A = 'clustering index (< 1: same colour closer)';
    lab_C = 'clustering index (colour vs other colours)';
end

%% Figure 1 : couleurs regroupées dans l'espace
% A : chaque session rangée par animal ; B : tous les animaux (moyenne +/- sem,
% demandé le 5 nov.) ; C : chaque couleur
% (5 nov.) 2 lignes : A, B, C en haut ; D (disques) et E (anneaux) en bas
fig1 = figure('Color', 'w', 'Position', [40 40 1800 1000]);
tl = tiledlayout(fig1, 2, 10, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Spatial clustering of Brainbow colours (%d animals, %d sessions)', ...
    numel(esp.animaux), numel(esp.champs)), 'FontWeight', 'bold', 'FontSize', 12);

ax = nexttile(tl, [1 4]); hold(ax, 'on');
% Un point = une session = un champ de vue (un même plan revu un autre jour
% n'est en général pas à la même profondeur) ; trait = l'animal.
animaux = esp.animaux;
nA = numel(animaux);
[~, ia] = ismember(esp.animal_du_champ, animaux);
rng(1);
for a = 1:nA
    m = find(ia == a);
    if isempty(m), continue; end
    y = exp(esp.champ(m));
    x = a + 0.16 * (rand(numel(m), 1) - 0.5) * 2;
    if comptage, sig = esp.z_champ(m) > 1.6449; else, sig = esp.z_champ(m) < -1.6449; end
    scatter(ax, x(~sig), y(~sig), 40, GRIS, 'LineWidth', 1);
    scatter(ax, x(sig), y(sig), 40, GRIS, 'filled');
    plot(ax, a + [-0.25 0.25], exp(esp.animal(a)) * [1 1], 'k-', 'LineWidth', 2);
end
yline(ax, 1, '--', 'Color', [0.6 0.6 0.6]);
set(ax, 'XTick', 1:nA, 'XTickLabel', animaux, 'XLim', [0.4 nA + 0.6], 'TickLabelInterpreter', 'none', ...
    'FontSize', 9, 'TickDir', 'out');
xtickangle(ax, 30);
ylabel(ax, lab_A);
title(ax, {'A. Each session, by animal', sprintf('%d/%d sessions p < 0.05', esp.n_champs_sig, ...
    numel(esp.champs))}, 'FontSize', 10);
box(ax, 'off');
yl_A = ylim(ax);

% B. tous les animaux : un point par animal, moyenne +/- sem
ax = nexttile(tl, [1 2]); hold(ax, 'on');
va = exp(esp.animal(:));
x = 1 + 0.10 * (rand(numel(va), 1) - 0.5) * 2;
scatter(ax, x, va, 40, GRIS, 'filled');
errorbar(ax, 1, mean(va), std(va) / sqrt(numel(va)), 'o', 'Color', 'k', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 7, 'LineWidth', 2, 'CapSize', 8);
yline(ax, 1, '--', 'Color', [0.6 0.6 0.6]);
set(ax, 'XTick', 1, 'XTickLabel', {'all animals'}, 'XLim', [0.5 1.5], 'YLim', yl_A, 'FontSize', 9, 'TickDir', 'out');
title(ax, {'B. All animals', sprintf('p = %s', texte_p(esp.p_wilcoxon))}, 'FontSize', 10);
box(ax, 'off');

% C. chaque couleur
ax = nexttile(tl, [1 4]); hold(ax, 'on');
for k = 1:nC
    y = exp(esp.c_animal{k});
    if isempty(y), continue; end
    x = k + 0.12 * (rand(numel(y), 1) - 0.5) * 2;
    scatter(ax, x, y, 45, cols(k, :), 'filled', 'MarkerFaceAlpha', 0.8);
    plot(ax, k + [-0.25 0.25], esp.c_indice_median(k) * [1 1], 'k-', 'LineWidth', 2);
end
yline(ax, 1, '--', 'Color', [0.6 0.6 0.6]);
yl = ylim(ax);
h_y = yl(2) - yl(1);
ylim(ax, [yl(1), yl(2) + 0.10 * h_y]);      % place pour les p, sous le titre
for k = 1:nC
    if isnan(esp.c_p(k)), continue; end
    text(ax, k, yl(2) + 0.09 * h_y, sprintf('p = %s', texte_p(esp.c_p(k))), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'FontSize', 8, 'Color', couleur_p(esp.c_p(k)));
end
set(ax, 'XTick', 1:nC, 'XTickLabel', noms_coul, 'XLim', [0.4 nC + 0.6], 'FontSize', 9, 'TickDir', 'out');
ylabel(ax, lab_C);
title(ax, 'C. Each colour (one dot per animal)', 'FontSize', 10);
box(ax, 'off');

% D. selon la distance : paires de même couleur observées / attendues si les
% couleurs étaient mélangées dans le champ (calcul de l'étape 3b, permutation
% 'globale', moyenne des animaux ; lu dans paires_agregation.mat s'il existe)
ax = nexttile(tl, [1 5]); hold(ax, 'on');
f_paires = fullfile(fileparts(fileparts(dossier)), 'paires_agregation', 'paires_agregation.mat');
Ag = [];
P = struct();
if isfile(f_paires)
    vars = intersect({'A', 'opts', 'MR', 'AniR', 'AR'}, who('-file', f_paires));
    P = load(f_paires, vars{:});
end
if isfield(P, 'A') && isfield(P, 'opts')
    Ag = P.A(strcmp(P.A.mesure, P.opts.mesures{1}) & strcmp(P.A.permutation, 'globale'), :);
end
if ~isempty(Ag)
    Agt = sortrows(Ag(strcmp(Ag.couleur, 'toutes'), :), {'dist_min_um', 'dist_max_um'});
    dist = [Agt.dist_min_um, Agt.dist_max_um];
    nd = size(dist, 1);
    if all(dist(:, 1) == 0)
        etiq = arrayfun(@(r) sprintf('%g', r), dist(:, 2), 'UniformOutput', false);
        etiq = strrep(etiq, 'Inf', 'all');
        xlab = 'pairs closer than ... \mum';
    else
        etiq = arrayfun(@(i) sprintf('%g-%g', dist(i, 1), dist(i, 2)), (1:nd)', 'UniformOutput', false);
        etiq = strrep(etiq, '-Inf', '+');
        xlab = 'distance between cells (\mum)';
    end
    h = gobjects(0); t = {};
    for k = 1:nC
        Agc = sortrows(Ag(strcmp(Ag.couleur, noms_coul{k}), :), {'dist_min_um', 'dist_max_um'});
        if isempty(Agc), continue; end
        [~, ix] = ismember([Agc.dist_min_um, Agc.dist_max_um], dist, 'rows');
        h(end+1) = plot(ax, ix, Agc.enrichissement, '-o', 'Color', cols(k, :), 'LineWidth', 1, ...
            'MarkerFaceColor', cols(k, :), 'MarkerSize', 3); %#ok<AGROW>
        t{end+1} = noms_coul{k}; %#ok<AGROW>
    end
    h(end+1) = plot(ax, 1:nd, Agt.enrichissement, '-o', 'Color', 'k', 'LineWidth', 2.5, ...
        'MarkerFaceColor', 'k', 'MarkerSize', 5);
    t{end+1} = 'all colours';
    yline(ax, 1, '--', 'Color', [0.6 0.6 0.6]);
    legend(ax, h, t, 'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
    set(ax, 'XTick', 1:nd, 'XTickLabel', etiq, 'XLim', [0.5 nd + 0.5], 'FontSize', 9, 'TickDir', 'out');
    xlabel(ax, xlab);
    ylabel(ax, 'same-colour pairs: observed / expected');
    title(ax, {'D. Same-colour pairs vs distance', '(mean of animals)'}, 'FontSize', 10);
else
    text(ax, 0.5, 0.5, {'Run step 3b (main\_paires) and its', 'aggregation to get this panel'}, ...
        'Units', 'normalized', 'HorizontalAlignment', 'center');
    axis(ax, 'off');
end
box(ax, 'off');

% E. par anneau (sans cumuler), 0-500 µm comme le test ; demandé le 5 nov.
% Valeurs et test : regroupement_anneaux (agregation_paires).
ax = nexttile(tl, [1 5]); hold(ax, 'on');
if isfield(P, 'AR') && ~isempty(P.AR)
    dans = @(T) T(isfinite(T.dist_max_um), :);           % l'anneau ouvert (> 500 µm) n'est pas testé
    ARt = sortrows(dans(P.AR(strcmp(P.AR.couleur, 'toutes'), :)), 'dist_min_um');
    centres = (ARt.dist_min_um + ARt.dist_max_um) / 2;
    etiq_r = arrayfun(@(i) sprintf('%g-%g', ARt.dist_min_um(i), ARt.dist_max_um(i)), (1:height(ARt))', ...
        'UniformOutput', false);
    % chaque animal
    Ar = dans(P.AniR(strcmp(P.AniR.couleur, 'toutes'), :));
    ua = unique(Ar.animal);
    for a = 1:numel(ua)
        Ra = sortrows(Ar(strcmp(Ar.animal, ua{a}), :), 'dist_min_um');
        plot(ax, (Ra.dist_min_um + Ra.dist_max_um) / 2, exp(Ra.lr), '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.8);
    end
    % chaque couleur (moyenne des animaux)
    h = gobjects(0); t = {};
    for k = 1:nC
        Rk = sortrows(dans(P.AR(strcmp(P.AR.couleur, noms_coul{k}), :)), 'dist_min_um');
        if isempty(Rk), continue; end
        h(end+1) = plot(ax, (Rk.dist_min_um + Rk.dist_max_um) / 2, exp(Rk.lr_moyen), '-o', 'Color', cols(k, :), ...
            'LineWidth', 1, 'MarkerFaceColor', cols(k, :), 'MarkerSize', 3); %#ok<AGROW>
        t{end+1} = noms_coul{k}; %#ok<AGROW>
    end
    % toutes couleurs : moyenne des animaux +/- sem (calculées en log)
    yc = exp(ARt.lr_moyen);
    h(end+1) = errorbar(ax, centres, yc, yc - exp(ARt.lr_moyen - ARt.lr_sem), exp(ARt.lr_moyen + ARt.lr_sem) - yc, ...
        '-o', 'Color', 'k', 'LineWidth', 2.5, 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'CapSize', 5);
    t{end+1} = 'all colours';
    % modèle mixte (toutes couleurs)
    titre_E = {'E. Same-colour pairs per ring'};
    if isfield(P, 'MR') && ~isempty(P.MR) && any(strcmp(P.MR.couleur, 'toutes'))
        r = P.MR(strcmp(P.MR.couleur, 'toutes'), :);
        dd = logspace(log10(min(centres)), log10(max(centres)), 50);
        h(end+1) = plot(ax, dd, exp(r.est_100um + r.pente * log2(dd / 100)), '--', 'Color', [0.75 0.1 0.1], ...
            'LineWidth', 1.5);
        t{end+1} = 'mixed model';
        titre_E = [titre_E, {sprintf('100 \\mum: %+.0f%%, p = %s', 100 * (exp(r.est_100um) - 1), texte_p(r.p_100um)), ...
            sprintf('per doubling: \\times%.2f, p = %s', exp(r.pente), texte_p(r.p_pente))}];
    end
    yline(ax, 1, '--', 'Color', [0.6 0.6 0.6]);
    legend(ax, h, t, 'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
    set(ax, 'XScale', 'log', 'XTick', centres, 'XTickLabel', etiq_r, ...
        'XLim', [min(centres) / 1.4, max(centres) * 1.4], 'FontSize', 9, 'TickDir', 'out');
    xlabel(ax, 'distance between cells (\mum)');
    ylabel(ax, 'same-colour pairs: observed / expected');
    title(ax, titre_E, 'FontSize', 10);
else
    text(ax, 0.5, 0.5, {'Run the pairs aggregation again', '(agregation\_paires) to get this panel'}, ...
        'Units', 'normalized', 'HorizontalAlignment', 'center');
    axis(ax, 'off');
end
box(ax, 'off');
exportgraphics(fig1, fullfile(dossier, 'colour_spatial_clustering.png'), 'Resolution', 200);

%% Figure 2 : détection et activité selon la couleur
nM = numel(det.mesures);
fig2 = figure('Color', 'w', 'Position', [80 40 1400 760]);
tl = tiledlayout(fig2, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, {'Does the colour change detection or activity? (one dot per animal; 1 = same as the coloured cells of the field)', ...
    'p = Wilcoxon on animals, not corrected (for a control, the cautious choice)'}, ...
    'FontWeight', 'bold', 'FontSize', 12);
for m = 1:nM
    ax = nexttile(tl); hold(ax, 'on');
    for k = 1:nC
        y = 2 .^ det.animal{m, k};
        if isempty(y), continue; end
        x = k + 0.12 * (rand(numel(y), 1) - 0.5) * 2;
        scatter(ax, x, y, 40, cols(k, :), 'filled', 'MarkerFaceAlpha', 0.8);
        plot(ax, k + [-0.25 0.25], det.rapport_median(m, k) * [1 1], 'k-', 'LineWidth', 2);
    end
    yline(ax, 1, '--', 'Color', [0.6 0.6 0.6]);
    set(ax, 'YScale', 'log');
    yl = ylim(ax);
    yl = [min(yl(1), 0.5), max(yl(2), 2)];
    ylim(ax, yl);
    for k = 1:nC
        if isnan(det.p(m, k)), continue; end
        text(ax, k, yl(2), texte_p(det.p(m, k)), 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top', 'FontSize', 8, 'Color', couleur_p(det.p(m, k)));
    end
    set(ax, 'XTick', 1:nC, 'XTickLabel', noms_coul, 'XLim', [0.4 nC + 0.6], 'FontSize', 9, 'TickDir', 'out');
    ylabel(ax, 'colour / all coloured cells');
    title(ax, sprintf('%s. %s', char('A' + m - 1), det.noms{m}), 'FontSize', 10);
    if strcmp(det.mesures{m}, 'membre')
        subtitle(ax, '(detection bias OR clonal effect)', 'FontSize', 8);
    end
    box(ax, 'off');
end
ax = nexttile(tl); hold(ax, 'on');
hb = bar(ax, 1:nC, [det.hyper_obs(:), det.hyper_attendu(:)], 'grouped', 'EdgeColor', 'none');
hb(1).FaceColor = 'flat'; hb(1).CData = cols;
hb(2).FaceColor = [0.82 0.82 0.82];
set(ax, 'XTick', 1:nC, 'XTickLabel', noms_coul, 'FontSize', 9, 'TickDir', 'out');
ylabel(ax, 'hyperactive cells (all sessions)');
legend(ax, {'observed', 'expected from colour share'}, 'Box', 'off', 'Location', 'northeast', 'FontSize', 8);
title(ax, sprintf('%s. Hyperactive cells: %d (%.1f%% of all cells)', char('A' + nM), det.n_hyper_tot, ...
    100 * det.n_hyper_tot / max(det.n_cell_tot, 1)), 'FontSize', 10);
box(ax, 'off');
exportgraphics(fig2, fullfile(dossier, 'detection_by_colour.png'), 'Resolution', 200);
end


function t = texte_p(p)
if isnan(p)
    t = 'NA';
elseif p < 0.001
    t = '< 0.001';
else
    t = sprintf('%.3f', p);
end
end


function c = couleur_p(p)
if p < 0.05
    c = [0.75 0.1 0.1];
else
    c = [0.45 0.45 0.45];
end
end
