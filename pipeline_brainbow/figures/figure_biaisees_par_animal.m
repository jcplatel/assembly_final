function fig = figure_biaisees_par_animal(dossier_variante, fdr, alpha, nom_variante)
%FIGURE_BIAISEES_PAR_ANIMAL  Assemblées biaisées de chaque animal, une par rond,
%   coloré par la couleur Brainbow en excès.
%
%   figure_biaisees_par_animal('E:\...\agregation\champ_globale')
%   (lancée aussi par l'étape 4 pour chaque variante)
%
%   Lit resultats_biais_clonal.mat (table R de l'étape 4 : assemblées
%   testables après les critères de taille et d'activité). Écrit
%   biased_assemblies_per_animal.png.
%     A. une colonne par animal, un rond par assemblée :
%          rond plein    biaisée (q FDR < fdr), couleur en excès
%          rond vide     p_minp < alpha seulement (pas après FDR)
%          trait noir    nombre attendu par hasard (somme des taux_nul)
%          sous l'axe    nombre d'assemblées testables
%     B. proportion à p_minp < alpha, observée contre attendue, par animal
%        (ce que compare le Wilcoxon sur les animaux, test A1)
%   Aucun test ici : seulement les résultats de l'étape 4.

if nargin < 2 || isempty(fdr), fdr = 0.10; end
if nargin < 3 || isempty(alpha), alpha = 0.05; end
if nargin < 4, nom_variante = ''; end
dossier_variante = regexprep(char(dossier_variante), '[\\/]+$', '');
S = load(fullfile(dossier_variante, 'resultats_biais_clonal.mat'), 'R');
R = S.R;
if isempty(nom_variante)
    [~, nom_variante] = fileparts(dossier_variante);
end
noms_coul = {'R', 'G', 'B', 'RG', 'RB', 'GB'};
cols = [0.86 0.16 0.16; 0.16 0.66 0.22; 0.20 0.38 0.92; ...
        0.93 0.78 0.10; 0.82 0.22 0.80; 0.10 0.74 0.80];
GRIS = [0.6 0.6 0.6];

Rt = R(R.testable, :);
animaux = unique(R.animal);
nA = numel(animaux);
n_t = zeros(nA, 1); n_att = zeros(nA, 1);
n_b = zeros(nA, 1); n_p = zeros(nA, 1); n_sig = zeros(nA, 1);
coul_b = cell(nA, 1); coul_p = cell(nA, 1);
for a = 1:nA
    m = strcmp(Rt.animal, animaux{a});
    n_t(a) = sum(m);
    n_att(a) = sum(Rt.taux_nul(m));
    b = m & Rt.q < fdr;
    p = m & Rt.p_minp < alpha & ~b;
    n_b(a) = sum(b); n_p(a) = sum(p);
    n_sig(a) = sum(m & Rt.p_minp < alpha);
    [~, coul_b{a}] = ismember(Rt.couleur_dominante(b), noms_coul);
    [~, coul_p{a}] = ismember(Rt.couleur_dominante(p), noms_coul);
end

fig = figure('Color', 'w', 'Position', [60 60 1350 560]);
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, {sprintf('Biased assemblies per animal: %d biased (FDR %g%%) + %d at p < %g only, out of %d testable (%.1f expected at p < %g by chance)', ...
    sum(n_b), 100 * fdr, sum(n_p), alpha, sum(n_t), sum(n_att), alpha), texte_variante(nom_variante)}, ...
    'FontWeight', 'bold', 'FontSize', 12, 'Interpreter', 'none');

% A. un rond par assemblée
ax = nexttile(tl, [1 2]); hold(ax, 'on');
ymax = max([n_b + n_p; n_att; 1]);
for a = 1:nA
    y = 0;
    for i = 1:n_b(a)
        y = y + 1;
        c = couleur_de(coul_b{a}(i), cols, GRIS);
        plot(ax, a, y, 'o', 'MarkerSize', 13, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
    end
    for i = 1:n_p(a)
        y = y + 1;
        c = couleur_de(coul_p{a}(i), cols, GRIS);
        plot(ax, a, y, 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', c, 'LineWidth', 2.2);
    end
    plot(ax, a + [-0.32 0.32], n_att(a) * [1 1], 'k-', 'LineWidth', 2);
    text(ax, a, -0.35, sprintf('n = %d', n_t(a)), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'FontSize', 8, 'Color', [0.35 0.35 0.35]);
end
set(ax, 'XTick', 1:nA, 'XTickLabel', animaux, 'XLim', [0.4 nA + 0.6], 'YLim', [-1.1 ymax + 1], ...
    'YTick', 0:max(1, ceil(ymax / 8)):ymax + 1, 'TickLabelInterpreter', 'none', 'FontSize', 10, 'TickDir', 'out');
ax.YGrid = 'on';
ylabel(ax, 'number of assemblies');
title(ax, 'A. One circle per assembly, coloured by the Brainbow colour in excess', 'FontSize', 10);
% légende
h = gobjects(0); t = {};
h(end+1) = plot(ax, nan, nan, 'o', 'MarkerSize', 11, 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'k');
t{end+1} = sprintf('biased (FDR q < %g)', fdr);
h(end+1) = plot(ax, nan, nan, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', [0.5 0.5 0.5], 'LineWidth', 2);
t{end+1} = sprintf('p < %g only', alpha);
h(end+1) = plot(ax, nan, nan, 'k-', 'LineWidth', 2);
t{end+1} = sprintf('expected at p < %g by chance', alpha);
for k = 1:6
    h(end+1) = plot(ax, nan, nan, 'o', 'MarkerSize', 9, 'MarkerFaceColor', cols(k, :), 'MarkerEdgeColor', 'none'); %#ok<AGROW>
    t{end+1} = noms_coul{k}; %#ok<AGROW>
end
legend(ax, h, t, 'Location', 'northoutside', 'Orientation', 'horizontal', 'NumColumns', 9, ...
    'Box', 'off', 'FontSize', 8);
box(ax, 'off');

% B. proportion observée contre attendue (test A1)
ax = nexttile(tl); hold(ax, 'on');
obs = 100 * n_sig ./ max(n_t, 1);
att = 100 * n_att ./ max(n_t, 1);
ok = n_t > 0;
for a = find(ok)'
    plot(ax, [1 2], [att(a) obs(a)], '-', 'Color', [0.7 0.7 0.7]);
end
scatter(ax, ones(sum(ok), 1), att(ok), 45, [0.45 0.45 0.45], 'filled');
scatter(ax, 2 * ones(sum(ok), 1), obs(ok), 45, [0.75 0.15 0.15], 'filled');
plot(ax, [0.8 1.2], mean(att(ok)) * [1 1], 'k-', 'LineWidth', 2);
plot(ax, [1.8 2.2], mean(obs(ok)) * [1 1], 'k-', 'LineWidth', 2);
set(ax, 'XTick', [1 2], 'XTickLabel', {'expected', 'observed'}, 'XLim', [0.5 2.5], ...
    'FontSize', 10, 'TickDir', 'out');
ylabel(ax, sprintf('%% of testable assemblies at p < %g', alpha));
title(ax, {'B. Each animal (one line)', 'what the Wilcoxon on animals compares'}, 'FontSize', 10);
box(ax, 'off');

exportgraphics(fig, fullfile(dossier_variante, 'biased_assemblies_per_animal.png'), 'Resolution', 200);
end


function c = couleur_de(i, cols, defaut)
if i >= 1 && i <= size(cols, 1)
    c = cols(i, :);
else
    c = defaut;
end
end


function t = texte_variante(nom)
principal = contains(nom, 'PRINCIPAL');
base = strtrim(regexprep(nom, '\(.*\)', ''));
base = regexprep(base, '_PROVISOIRE$', '');
switch base
    case 'champ_globale',       t = 'field reference, colors shuffled across the field';
    case 'champ_blocs',         t = 'spatial control: colors shuffled only between nearby cells';
    case 'champ_territoire',    t = 'territory control: reference = colored cells within the assembly''s own area';
    case 'champ_activite',      t = 'activity control: colors shuffled between cells of similar activity';
    otherwise,                  t = base;
end
if principal, t = ['MAIN TEST: ' t]; end
end
