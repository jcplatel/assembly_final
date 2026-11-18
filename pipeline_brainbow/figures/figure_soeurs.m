function fig = figure_soeurs(dossier)
%FIGURE_SOEURS  Étape 3d, sans rien recalculer : lit carte_soeurs.mat et
%   soeurs_agregation.mat dans <PathSave>\soeurs_agregation\.
%   En haut (couleurs et positions seulement) :
%     A. part de vraies sœurs parmi les paires de même couleur, selon la
%        distance (pi = 1 - 1 / (observé / attendu)) : chaque couleur, et
%        toutes couleurs en noir ; zone grise = au-delà de opts.pi_max_um,
%        pi mis à 0 dans le test (valeurs montrées avant cette limite) ;
%     B. paires de même couleur observées / attendues selon l'écart le long
%        des bandes de couleur et en travers (chaque session alignée sur ses
%        bandes, chaque couleur sur l'angle trouvé avec les autres couleurs ;
%        toutes sessions) : forme des clones.
%   En bas (une mesure d'activité par panneau) : dose-réponse. Résidu moyen
%   des paires de couleurs différentes, de même couleur avec pi = 0 (même
%   couleur mais trop loin pour être sœurs), puis par catégorie de pi. Gris :
%   chaque animal ; noir : moyenne des animaux +/- sem, toutes couleurs.
%   Couleurs (7 nov.) : la même chose couleur par couleur (moyenne des
%   animaux ; « diff. colour » = paires entre une cellule de cette couleur
%   et une cellule d'une autre couleur ; points avec moins de 20 paires en
%   tout non montrés) ; légende : effet chez les sœurs de cette couleur
%   (descriptif). Titre : effet chez les vraies sœurs et Wilcoxon sur les
%   animaux (toutes couleurs, le test).
%   Écrit soeurs.png.

dossier = regexprep(char(dossier), '[\\/]+$', '');
K = load(fullfile(dossier, 'carte_soeurs.mat'));
S = load(fullfile(dossier, 'soeurs_agregation.mat'));
opts = S.opts;
cols = [0.86 0.16 0.16; 0.16 0.66 0.22; 0.20 0.38 0.92; ...
        0.93 0.78 0.10; 0.82 0.22 0.80; 0.10 0.74 0.80];
noms_c = {'R', 'G', 'B', 'RG', 'RB', 'GB'};
mesures = S.B.mesure;
nM = numel(mesures);
fig = figure('Color', 'w', 'Position', [40 40 1500 950]);
tl = tiledlayout(fig, 2, max(2, nM), 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Are probable sister pairs more alike in activity? (%d animals, %d sessions)', ...
    S.B.n_animaux(1), K.n_sessions), 'FontWeight', 'bold', 'FontSize', 12);

% A. pi selon la distance
ax = nexttile(tl); hold(ax, 'on');
nb = numel(opts.bornes_um) - 1;
etiq = arrayfun(@(i) texte_tranche(opts.bornes_um(i), opts.bornes_um(i + 1)), 1:nb, 'UniformOutput', false);
pi_cb = K.pi_cb; pi_b = K.pi_b;
if isfield(K, 'pi_brut_cb'), pi_cb = K.pi_brut_cb; pi_b = K.pi_brut_b; end
debut = find(opts.bornes_um(1:nb) >= opts.pi_max_um, 1);
if ~isempty(debut)
    haut_y = 100 * max([pi_cb(:); pi_b(:); 0.1]) * 1.1;
    patch(ax, [debut - 0.5, nb + 0.5, nb + 0.5, debut - 0.5], [0 0 haut_y haut_y], [0.92 0.92 0.92], ...
        'EdgeColor', 'none');
    text(ax, debut - 0.4, haut_y * 0.5, {sprintf('> %g \\mum:', opts.pi_max_um), '\pi set to 0'}, ...
        'VerticalAlignment', 'middle', 'FontSize', 8, 'Color', [0.4 0.4 0.4]);
end
h = gobjects(0); t = {};
for c = 1:6
    h(end+1) = plot(ax, 1:nb, 100 * pi_cb(:, c), '-o', 'Color', cols(c, :), 'MarkerFaceColor', cols(c, :), ...
        'MarkerSize', 4, 'LineWidth', 1); %#ok<AGROW>
    t{end+1} = noms_c{c}; %#ok<AGROW>
end
h(end+1) = plot(ax, 1:nb, 100 * pi_b, '-o', 'Color', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 2.5);
t{end+1} = 'all colours';
legend(ax, h, t, 'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
set(ax, 'XTick', 1:nb, 'XTickLabel', etiq, 'XLim', [0.5 nb + 0.5], 'FontSize', 9, 'TickDir', 'out');
xtickangle(ax, 30);
xlabel(ax, 'distance between the two cells (\mum)');
ylabel(ax, 'probable sisters among same-colour pairs (%)');
title(ax, 'A. Probability of being sisters (colours and positions only)', 'FontSize', 10);
box(ax, 'off');

% B. carte 2D alignée sur les bandes
ax = nexttile(tl);
nk = size(K.r_uv, 1);
centres = ((1:nk) - 0.5) * opts.carte_pas_um;
imagesc(ax, centres, centres, K.r_uv);
set(ax, 'YDir', 'normal', 'FontSize', 9, 'TickDir', 'out');
axis(ax, 'image');
colormap(ax, parula);
v = K.r_uv(isfinite(K.r_uv));
if ~isempty(v), caxis(ax, [min(0.9, min(v)), max(1.1, prctile(v, 98))]); end
cb = colorbar(ax);
cb.Label.String = 'same-colour pairs: observed / expected';
xlabel(ax, 'distance along the colour bands (\mum)');
ylabel(ax, 'distance across the bands (\mum)');
title(ax, 'B. Shape of clones (sessions aligned on their bands)', 'FontSize', 10);
for k = 3:max(2, nM), nexttile(tl); axis off; end

% C... dose-réponse par mesure
noms_cat = [{'diff. colour', 'same, \pi = 0'}, arrayfun(@(i) texte_pi(opts.categories_pi(i), ...
    opts.categories_pi(i + 1)), 1:numel(opts.categories_pi) - 1, 'UniformOutput', false)];
nC = numel(noms_cat);
lettres = 'CDEF';
for k = 1:nM
    ax = nexttile(tl); hold(ax, 'on');
    a = S.Ad(strcmp(S.Ad.mesure, mesures{k}), :);
    Y = nan(height(a), nC);
    for r = 1:height(a)
        y = a.moyenne{r};
        n = a.n{r};
        y(n == 0) = NaN;
        Y(r, :) = y;
        plot(ax, 1:nC, y, '-', 'Color', [0.75 0.75 0.75]);
    end
    hc = gobjects(0); tc = {};
    if ismember('moyenne_c', a.Properties.VariableNames)
        Mc = cat(3, a.moyenne_c{:});                 % couleur x catégorie x animal
        Nc = cat(3, a.n_c{:});
        Mc(Nc == 0) = NaN;
        mc = mean(Mc, 3, 'omitnan');
        mc(sum(Nc, 3) < 20) = NaN;
        for c = 1:6
            if all(isnan(mc(c, :))), continue; end
            hc(end+1) = plot(ax, 1:nC, mc(c, :), '-o', 'Color', cols(c, :), 'MarkerFaceColor', cols(c, :), ...
                'MarkerSize', 3, 'LineWidth', 1); %#ok<AGROW>
            eff = NaN; bet = NaN;
            if isfield(S, 'Bc') && ~isempty(S.Bc)
                r = strcmp(S.Bc.mesure, mesures{k}) & strcmp(S.Bc.couleur, noms_c{c});
                if any(r)
                    eff = S.Bc.effet_pct(find(r, 1));
                    if ismember('beta_moyen', S.Bc.Properties.VariableNames), bet = S.Bc.beta_moyen(find(r, 1)); end
                end
            end
            tc{end+1} = sprintf('%s: sisters %s', noms_c{c}, texte_effet(eff, bet)); %#ok<AGROW>
        end
    end
    nA = sum(isfinite(Y), 1);
    hk = errorbar(ax, 1:nC, mean(Y, 1, 'omitnan'), std(Y, 0, 1, 'omitnan') ./ sqrt(max(nA, 1)), '-o', ...
        'Color', 'k', 'LineWidth', 2, 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'CapSize', 5);
    if ~isempty(hc)
        legend(ax, [hc, hk], [tc, {'all colours'}], 'Location', 'best', 'Box', 'off', 'FontSize', 7, ...
            'AutoUpdate', 'off');
    end
    yline(ax, 0, ':', 'Color', [0.6 0.6 0.6]);
    set(ax, 'XTick', 1:nC, 'XTickLabel', noms_cat, 'XLim', [0.5 nC + 0.5], 'FontSize', 9, 'TickDir', 'out');
    xtickangle(ax, 30);
    ax.YAxis.Exponent = 0;
    xlabel(ax, 'pair category (\pi = probability of being sisters)');
    ylabel(ax, 'activity similarity (residual)');
    b = S.B(strcmp(S.B.mesure, mesures{k}), :);
    bet = NaN;
    if ismember('beta_moyen', b.Properties.VariableNames), bet = b.beta_moyen; end
    title(ax, {sprintf('%s. %s', lettres(k), nom_mesure(mesures{k})), sprintf('sisters: %s, %s', ...
        texte_effet(b.effet_pct, bet), texte_p(b.p_wilcoxon))}, 'FontSize', 10);
    box(ax, 'off');
end
exportgraphics(fig, fullfile(dossier, 'soeurs.png'), 'Resolution', 200);
end


function t = texte_effet(pct, beta)
% en % si la mesure a une échelle positive, sinon la différence brute (phi des SCE)
if isfinite(pct), t = sprintf('%+.0f%%', pct); else, t = sprintf('\\Delta\\phi = %+.4f', beta); end
end


function t = nom_mesure(m)
switch m
    case 'sttc',        t = 'STTC at rest';
    case 'coassemblee', t = 'Same assembly';
    case 'sce',         t = 'SCE co-participation';
    otherwise,          t = m;
end
end


function t = texte_tranche(a, b)
if isinf(b), t = sprintf('>%g', a); else, t = sprintf('%g-%g', a, b); end
end


function t = texte_pi(a, b)
if isinf(b), t = sprintf('\\pi > %g%%', 100 * a); else, t = sprintf('\\pi %g-%g%%', 100 * a, 100 * b); end
end


function t = texte_p(p)
if isnan(p)
    t = 'p = NA';
elseif p < 0.001
    t = 'p < 0.001';
else
    t = sprintf('p = %.3f', p);
end
end
