function figure_angles_bandes(sessions, dossier)
%FIGURE_ANGLES_BANDES  Étape 3d : vérifier à l'œil la direction des bandes de
%   couleur trouvée dans chaque session (etape_soeurs_session). Une figure
%   par animal, un panneau par session : cellules colorées à leur place
%   (comme l'image GCaMP de l'interface, y vers le bas), grises = non
%   colorées ; trait noir = direction trouvée (toutes les couleurs).
%   Titre : angle (0° = horizontale de l'image), paires de même couleur
%   observées / attendues le long des bandes et en travers (chaque couleur
%   mesurée sur l'angle trouvé avec les autres couleurs : contrôle croisé ;
%   le long > en travers = direction confirmée par des clones différents).
%   Écrit angles_bandes_<animal>.png dans dossier.

dossier = regexprep(char(dossier), '[\\/]+$', '');
cols = [0.86 0.16 0.16; 0.16 0.66 0.22; 0.20 0.38 0.92; ...
        0.93 0.78 0.10; 0.82 0.22 0.80; 0.10 0.74 0.80];
animaux = unique(sessions.id_animal, 'stable');
for a = 1:numel(animaux)
    lignes = find(strcmp(sessions.id_animal, animaux{a}));
    R = {};
    for s = lignes(:)'
        f = fullfile(sessions.dossier{s}, 'soeurs', 'soeurs.mat');
        if ~isfile(f), continue; end
        r = load(f, 'col', 'xy_um', 'comptes', 'nom');
        if ~isfield(r.comptes, 'angle'), continue; end
        R{end+1} = r; %#ok<AGROW>
    end
    if isempty(R), continue; end
    n = numel(R);
    nc = ceil(sqrt(n));
    nl = ceil(n / nc);
    fig = figure('Color', 'w', 'Position', [40 40 min(1800, 420 * nc) min(1000, 400 * nl + 60)]);
    tl = tiledlayout(fig, nl, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('%s: direction of the colour bands in each session', animaux{a}), ...
        'FontWeight', 'bold', 'FontSize', 12, 'Interpreter', 'none');
    for k = 1:n
        r = R{k};
        ax = nexttile(tl); hold(ax, 'on');
        col = r.col(:);
        xy = r.xy_um;
        gris = ~(col >= 1 & col <= 6);
        scatter(ax, xy(gris, 1), xy(gris, 2), 4, [0.82 0.82 0.82], 'filled');
        for c = 1:6
            m = col == c;
            scatter(ax, xy(m, 1), xy(m, 2), 12, cols(c, :), 'filled');
        end
        A = r.comptes.angle;
        centre = mean(xy, 1, 'omitnan');
        L = 0.45 * max(max(xy, [], 1) - min(xy, [], 1));
        t = A.theta_deg * pi / 180;
        plot(ax, centre(1) + L * [-cos(t), cos(t)], centre(2) + L * [-sin(t), sin(t)], '-', 'Color', 'k', ...
            'LineWidth', 2);
        set(ax, 'YDir', 'reverse', 'FontSize', 8, 'TickDir', 'out');
        axis(ax, 'equal');
        axis(ax, 'tight');
        title(ax, {strrep(r.nom, '_', '\_'), sprintf('%.0f\\circ; same colour along %.2f / across %.2f', ...
            A.theta_deg, A.long, A.travers)}, 'FontSize', 8, 'Interpreter', 'tex');
        xlabel(ax, '\mum');
        box(ax, 'on');
    end
    exportgraphics(fig, fullfile(dossier, sprintf('angles_bandes_%s.png', animaux{a})), 'Resolution', 150);
    close(fig);
end
end

