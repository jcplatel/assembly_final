function fig = figure_arbre_session(dossier)
%FIGURE_ARBRE_SESSION  Étape 3c, UNE session, sans rien recalculer : lit
%   <session>\hierarchique\hierarchique.mat.
%   Une ligne par arbre :
%     à gauche, l'arbre ; sous les feuilles, la couleur Brainbow de chaque
%     cellule (gris clair = non colorée) ; pointillés gris = les 4 coupes
%     (groupes de ~5, ~10, ~20, ~40 cellules) ; rond = groupe d'une coupe
%     SIGNIFICATIF après FDR (couleurs mélangées dans le champ), rempli de
%     sa couleur en excès (plus gros = coupe plus haute) ; losange = cluster
%     FCA significatif ; pointillé rouge = coupe de la FCA ;
%     à droite, par coupe : groupes avec une couleur en excès (p < alpha),
%     observés contre attendus (couleurs mélangées dans le champ ; par
%     blocs) ; * = p < 0,05 (champ) ; chiffre = groupes significatifs après FDR.
%   Écrit arbres.png dans le même dossier.

dossier = regexprep(char(dossier), '[\\/]+$', '');
H = load(fullfile(dossier, 'hierarchique.mat'));
cols = [0.86 0.16 0.16; 0.16 0.66 0.22; 0.20 0.38 0.92; ...
        0.93 0.78 0.10; 0.82 0.22 0.80; 0.10 0.74 0.80];
fig = figure('Color', 'w', 'Position', [40 40 1800 900]);
if ~isfield(H, 'arbres') || isempty(H.arbres)
    annotation(fig, 'textbox', [0.1 0.4 0.8 0.2], 'String', {H.nom, H.statut}, 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'Interpreter', 'none', 'FontSize', 12);
    exportgraphics(fig, fullfile(dossier, 'arbres.png'), 'Resolution', 150);
    return
end
opts = H.opts;
nA = numel(H.arbres);
tl = tiledlayout(fig, nA, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s: %d cells in >= %d SCE (%d coloured)', H.nom, numel(H.cellules), opts.min_sce, ...
    H.n_colorees), 'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 11);
for a = 1:nA
    A = H.arbres(a);
    Z = A.Z;
    N = size(Z, 1) + 1;
    xf = zeros(1, N);
    xf(A.ordre) = 1:N;
    h = A.hauteur;
    % arbre
    ax = nexttile(tl, [1 3]); hold(ax, 'on');
    X = nan(5 * (N - 1), 1);
    Y = nan(5 * (N - 1), 1);
    for i = 1:N - 1
        [x1, y1] = position(Z(i, 1), N, xf, A.x, h);
        [x2, y2] = position(Z(i, 2), N, xf, A.x, h);
        X(5 * i - 4:5 * i - 1) = [x1; x1; x2; x2];
        Y(5 * i - 4:5 * i - 1) = [y1; h(i); h(i); y2];
    end
    plot(ax, X, Y, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.5);
    hmin = min(h);
    hmax = max(h);
    plage = hmax - min(0, hmin);
    y0 = min(0, hmin) - 0.03 * plage;
    y1b = y0 - 0.07 * plage;
    % bande de couleurs sous les feuilles
    c = H.couleur(A.ordre);
    C = repmat([0.88 0.88 0.88], N, 1);
    ok = c >= 1 & c <= 6;
    C(ok, :) = cols(c(ok), :);
    V = [reshape([(1:N) - 0.5; (1:N) + 0.5; (1:N) + 0.5; (1:N) - 0.5], [], 1), ...
        repmat([y0; y0; y1b; y1b], N, 1)];
    F = reshape(1:4 * N, 4, N)';
    patch(ax, 'Faces', F, 'Vertices', V, 'FaceVertexCData', C, 'FaceColor', 'flat', 'EdgeColor', 'none');
    % coupes et groupes significatifs après FDR (hasard principal)
    n_sig = 0;
    coupes = A.coupes;
    if isfield(A, 'amas') && ~isempty(A.amas), coupes = [coupes, A.amas]; end
    for iv = 1:numel(coupes)
        Cv = coupes(iv);
        if Cv.niveau > 0
            yline(ax, Cv.hauteur, ':', sprintf('~%g', Cv.niveau), 'Color', [0.55 0.55 0.55], ...
                'LabelHorizontalAlignment', 'right', 'FontSize', 8);
        end
        R = Cv.tests(1);
        g = find(R.fdr & Cv.noeud > 0);
        n_sig = n_sig + numel(g);
        if isempty(g), continue; end
        nd = Cv.noeud(g);
        if Cv.niveau == 0
            scatter(ax, A.x(nd), h(nd), 70, cols(R.couleur(g), :), 'd', 'filled', 'MarkerEdgeColor', 'k');
        else
            scatter(ax, A.x(nd), h(nd), 15 + 15 * iv, cols(R.couleur(g), :), 'filled', ...
                'MarkerEdgeColor', [0.15 0.15 0.15]);
        end
    end
    set(ax, 'XLim', [0.5, N + 0.5], 'YLim', [y1b - 0.01 * plage, hmax + 0.02 * plage], 'XTick', [], ...
        'FontSize', 9, 'TickDir', 'out');
    switch A.nom
        case 'sce',  titre = 'SCE co-participation'; yl = 'distance (1 - phi)';
        case 'sttc', titre = 'STTC residual at rest'; yl = 'distance (max - residual)';
        case 'dff',  titre = 'Pearson dF/F residual at rest'; yl = 'distance (max - residual)';
        case 'fca',  titre = 'FCA (Feldt et al.)'; yl = 'merging step';
        otherwise
            titre = strrep(strrep(A.nom, 'sttc_', 'STTC residual, \pm'), 'ms', ' ms');
            yl = 'distance (max - residual)';
    end
    ylabel(ax, yl);
    if strcmp(A.nom, 'fca') && ~isempty(A.fca)
        yline(ax, A.fca.coupe + 0.5, '--', 'Color', [0.8 0.1 0.1], 'LineWidth', 1.5);
        titre = sprintf('%s, cut after %d merges, %d clusters', titre, A.fca.coupe, numel(A.fca.amas));
    end
    title(ax, sprintf('%s: %d groups significant after FDR %g%%', titre, n_sig, 100 * opts.fdr_q), ...
        'FontSize', 10);
    box(ax, 'off');
    % observés contre attendus, par coupe
    ax = nexttile(tl); hold(ax, 'on');
    nC = numel(coupes);
    nT = numel(coupes(1).tests);
    Y = zeros(nC, 1 + nT);
    noms = {'observed'};
    for t = 1:nT
        noms{end+1} = sprintf('expected (%s)', texte_perm(coupes(1).tests(t).permutation)); %#ok<AGROW>
    end
    etiq = cell(1, nC);
    for iv = 1:nC
        Y(iv, 1) = coupes(iv).tests(1).obs;
        for t = 1:nT, Y(iv, 1 + t) = coupes(iv).tests(t).attendu; end
        if coupes(iv).niveau == 0, etiq{iv} = 'FCA'; else, etiq{iv} = sprintf('~%g', coupes(iv).niveau); end
    end
    b = bar(ax, 1:nC, Y, 'grouped', 'EdgeColor', 'none');
    gris = [0.15 0.15 0.15; 0.65 0.65 0.65; 0.85 0.85 0.85];
    for k = 1:numel(b), b(k).FaceColor = gris(min(k, 3), :); end
    for iv = 1:nC
        R = coupes(iv).tests(1);
        etoile = '';
        if R.p < 0.05, etoile = '*'; end
        text(ax, iv, max(Y(iv, :)), {etoile, sprintf('%d', R.n_fdr)}, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'FontSize', 10);
    end
    ylim(ax, [0, max(Y(:)) * 1.35 + 1]);
    set(ax, 'XTick', 1:nC, 'XTickLabel', etiq, 'FontSize', 9, 'TickDir', 'out');
    xlabel(ax, 'cut (mean group size, cells)');
    ylabel(ax, 'groups with a colour in excess');
    legend(ax, b, noms, 'Location', 'northwest', 'Box', 'off', 'FontSize', 8);
    title(ax, {'* p < 0.05 (field); number = significant after FDR'}, 'FontSize', 9, 'FontWeight', 'normal');
    box(ax, 'off');
end
exportgraphics(fig, fullfile(dossier, 'arbres.png'), 'Resolution', 150);
end


function [x, y] = position(c, N, xf, x_noeud, h)
if c <= N
    x = xf(c);
    y = 0;
else
    x = x_noeud(c - N);
    y = h(c - N);
end
end


function t = texte_perm(p)
switch p
    case 'globale', t = 'field';
    case 'blocs',   t = 'blocks';
    otherwise,      t = p;
end
end
