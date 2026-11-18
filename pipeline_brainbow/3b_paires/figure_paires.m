function figure_paires(res, opts, fichier_png, titre)
%FIGURE_PAIRES  delta (même couleur - couleurs différentes) selon la distance,
%   pour chaque couleur et toutes couleurs, une colonne par mesure
%   (permutation principale). Point plein : p < 0.05.

noms = {'R', 'G', 'B', 'RG', 'RB', 'GB', 'toutes'};
rgb = [1 0.15 0.15; 0.1 0.75 0.1; 0.2 0.4 1; 0.85 0.75 0; 0.9 0.2 0.9; 0 0.75 0.75; 0 0 0];
[bornes, ~, etiquettes] = distances_paires(opts);   % disques ou tranches
nb = numel(bornes) - 1;
mesures = liste_mesures_paires(opts);
nM = numel(mesures);
fig = figure('Visible', 'off', 'Position', [100 100 520 * nM 420]);
for im = 1:nM
    subplot(1, nM, im); hold on
    sel = strcmp(res.mesure, mesures{im}) & strcmp(res.permutation, opts.permutations{1});
    R = res(sel, :);
    h = gobjects(1, numel(noms));
    for c = 1:numel(noms)
        Rc = R(strcmp(R.couleur, noms{c}), :);
        [~, ordre] = sortrows([Rc.dist_min_um, Rc.dist_max_um]);
        Rc = Rc(ordre, :);
        epaisseur = 1; if c == numel(noms), epaisseur = 2.5; end
        h(c) = plot(1:height(Rc), Rc.delta, '-', 'Color', rgb(c, :), 'LineWidth', epaisseur);
        sig = Rc.p < 0.05;
        plot(find(sig), Rc.delta(sig), 'o', 'MarkerFaceColor', rgb(c, :), 'MarkerEdgeColor', rgb(c, :));
        plot(find(~sig), Rc.delta(~sig), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', rgb(c, :));
    end
    yline(0, ':k');
    set(gca, 'XTick', 1:nb, 'XTickLabel', etiquettes);
    xlabel('Distance entre cellules (µm)');
    if isfield(opts, 'distances') && strcmp(opts.distances, 'disques')
        xlabel('Paires à moins de ... µm');
    end
    ylabel('\Delta similarité (même couleur - différentes)');
    title(sprintf('%s (permutation %s)', mesures{im}, opts.permutations{1}), 'Interpreter', 'none');
    legend(h, noms, 'Location', 'northeast', 'Box', 'off');
end
sgtitle(sprintf('%s : point plein = p < 0.05', titre), 'Interpreter', 'none');
exportgraphics(fig, fichier_png, 'Resolution', 150);
close(fig)
end
