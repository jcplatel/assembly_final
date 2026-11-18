function fig = figure_hierarchique(dossier)
%FIGURE_HIERARCHIQUE  Étape 3c, toutes les sessions, sans rien recalculer :
%   lit <PathSave>\hierarchique_agregation\hierarchique_agregation.mat.
%   Une ligne par arbre, une colonne par hasard (couleurs mélangées dans le
%   champ ; par blocs). x : coupe de l'arbre (groupes de ~5, ~10, ~20, ~40
%   cellules en moyenne ; clusters FCA) ; y : % des groupes testés avec une
%   couleur en excès. Gris : chaque animal (observé) ; noir : moyenne des
%   animaux +/- sem ; pointillé : attendu (couleurs mélangées).
%   En haut : p (Wilcoxon unilatéral sur les animaux, sans correction) et
%   nombre de groupes significatifs après FDR / groupes testés.
%   Écrit hierarchique.png.

dossier = regexprep(char(dossier), '[\\/]+$', '');
S = load(fullfile(dossier, 'hierarchique_agregation.mat'));
A = S.A; B = S.B; opts = S.opts;
arbres = unique(B.arbre, 'stable');
perms = unique(B.permutation, 'stable');
GRIS = [0.75 0.75 0.75];
fig = figure('Color', 'w', 'Position', [40 40 1400 120 + 320 * numel(arbres)]);
tl = tiledlayout(fig, numel(arbres), numel(perms), 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Hierarchical clustering of cells: groups with a colour in excess (%d animals)', ...
    numel(unique(A.animal))), 'FontWeight', 'bold', 'FontSize', 12);
for a = 1:numel(arbres)
    for q = 1:numel(perms)
        ax = nexttile(tl); hold(ax, 'on');
        Aa = A(strcmp(A.arbre, arbres{a}) & strcmp(A.permutation, perms{q}), :);
        Bb = sortrows(B(strcmp(B.arbre, arbres{a}) & strcmp(B.permutation, perms{q}), :), 'niveau');
        niv = Bb.niveau(:)';
        nN = numel(niv);
        [~, ia] = ismember(Aa.niveau, niv);
        for an = unique(Aa.animal)'
            m = strcmp(Aa.animal, an{1});
            [xs, o] = sort(ia(m));
            y = Aa.pct_observe(m);
            plot(ax, xs, y(o), '-o', 'Color', GRIS, 'MarkerSize', 3, 'MarkerFaceColor', GRIS);
        end
        plot(ax, 1:nN, Bb.pct_attendu, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.5);
        errorbar(ax, 1:nN, Bb.pct_observe, Bb.pct_sem, '-o', 'Color', 'k', 'LineWidth', 2, ...
            'MarkerFaceColor', 'k', 'MarkerSize', 5, 'CapSize', 5);
        ymax = max([Aa.pct_observe; Bb.pct_observe + Bb.pct_sem; Bb.pct_attendu; 1]) * 1.2;
        for i = 1:nN
            couleur_p = [0.45 0.45 0.45];
            if Bb.p_wilcoxon(i) < 0.05, couleur_p = [0.80 0.10 0.10]; end
            text(ax, i, ymax, {texte_p(Bb.p_wilcoxon(i)), sprintf('FDR %d/%d', Bb.n_fdr(i), Bb.n_groupes(i))}, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 8, 'Color', couleur_p);
        end
        etiq = arrayfun(@texte_niveau, niv, 'UniformOutput', false);
        set(ax, 'XTick', 1:nN, 'XTickLabel', etiq, 'XLim', [0.5 nN + 0.5], 'YLim', [0 ymax * 1.05], ...
            'FontSize', 9, 'TickDir', 'out');
        xlabel(ax, 'cut of the tree (mean group size, cells)');
        ylabel(ax, '% of groups with a colour in excess');
        title(ax, sprintf('%s, %s', nom_arbre(arbres{a}), nom_perm(perms{q}, opts)), 'FontSize', 10);
        box(ax, 'off');
    end
end
exportgraphics(fig, fullfile(dossier, 'hierarchique.png'), 'Resolution', 200);
end


function t = nom_arbre(a)
switch a
    case 'sce',  t = 'SCE co-participation';
    case 'sttc', t = 'STTC residual';
    case 'dff',  t = 'Pearson dF/F residual';
    case 'fca',  t = 'FCA tree';
    case 'fca_amas', t = 'FCA clusters only';
    otherwise,   t = strrep(strrep(a, 'sttc_', 'STTC residual \pm'), 'ms', ' ms');
end
end


function t = nom_perm(p, opts)
switch p
    case 'globale', t = 'colours shuffled in the field';
    case 'blocs',   t = sprintf('shuffled in %g \\mum blocks', opts.taille_bloc_um);
    otherwise,      t = p;
end
end


function t = texte_niveau(n)
if n == 0, t = 'FCA clusters'; else, t = sprintf('~%g', n); end
end


function t = texte_p(p)
if isnan(p)
    t = '';
elseif p < 0.001
    t = 'p<0.001';
else
    t = sprintf('p=%.3f', p);
end
end
