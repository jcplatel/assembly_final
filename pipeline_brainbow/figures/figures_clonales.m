function figures_clonales(colorcell, assembly, calcium, outline_gcampx, outline_gcampy, res, namefull)
% Figures couleur d'un champ (reprises de brainbowassemblies2026_05_07 et
% de stats_clonales_champ2) :
%   color_all_cells.png         camembert des couleurs du champ
%   allcells.png                toutes les cellules sur l'image calcique
%   assemblyN.png               chaque assemblée sur l'image calcique
%   piechart_noblackN.png       camembert des couleurs de chaque assemblée
% res : table de la variante principale (p_minp, p_chi2, testable, ...)

couleurs_rgb = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1];
couleurs8    = [couleurs_rgb; 1 1 1; 0 0 0];   % + blanc (7) et noir (8)
nC = 6;
colorcell = colorcell(:);
nA = numel(assembly);

%% 1. Composition du champ
color_all = arrayfun(@(c) sum(colorcell == c), 1:8);
fig = figure('visible', 'off');
ax = axes(fig);
dessiner_camembert(ax, color_all, [], couleurs8);
titre_obj = sgtitle(fig, 'color All cells', 'FontSize', 12, 'FontWeight', 'bold');
titre_obj.HorizontalAlignment = 'left';
exportgraphics(fig, strcat(namefull, 'color_all_cells.png'), 'Resolution', 150, 'ContentType', 'image');
close(fig)

%% 2. Toutes les cellules sur l'image calcique
if ~isempty(calcium)
    fig = figure('visible', 'off');
    colormap(gray); imagesc(calcium); hold on
    for n = 1:numel(colorcell)
        cc = colorcell(n);
        if cc >= 1 && cc <= 8 && cc == round(cc)
            col = couleurs8(cc, :);
        else
            col = [0.5 0.5 0.5];        % gris = code couleur inconnu
        end
        dessiner_contour(outline_gcampx{n}, outline_gcampy{n}, col, 'k');
    end
    hold off; axis off; axis image;
    title('all cells', 'FontSize', 16);
    exportgraphics(fig, strcat(namefull, 'allcells.png'), 'Resolution', 150, ...
        'BackgroundColor', 'black', 'ContentType', 'image');
    close(fig)
end

%% 3. Une carte et un camembert par assemblée
for a = 1:nA
    if res.testable(a), suffixe = ''; else, suffixe = ' | non testable'; end

    % --- Carte de l'assemblée sur l'image calcique ---
    if ~isempty(calcium)
        fig = figure('visible', 'off');
        colormap(gray); imagesc(calcium); hold on
        for n = unique(assembly{a}(:))'
            cc = colorcell(n);
            if cc >= 1 && cc <= 6, col = couleurs_rgb(cc, :); else, col = [0 0 0]; end
            dessiner_contour(outline_gcampx{n}, outline_gcampy{n}, col, 'w');
        end
        hold off; axis off; axis image;
        title(sprintf('Assemblée %d : %d cellules (%d colorées) | p_{min-p} = %.4f%s', ...
            a, res.n_total(a), res.n_colorees(a), res.p_minp(a), suffixe), 'FontSize', 12);
        exportgraphics(fig, strcat(namefull, ['assembly', num2str(a), '.png']), ...
            'Resolution', 150, 'Padding', 30, 'ContentType', 'image');
        close(fig)
    end

    % --- Camembert des couleurs de l'assemblée ---
    if res.n_colorees(a) > 0
        Obs_a = res.compte(a, :);
        labels = cell(1, nC);
        for c = 1:nC
            if ~isnan(res.p_brut(a, c))
                labels{c} = sprintf('%.1f%% (p brut=%.3f)', 100 * Obs_a(c) / res.n_colorees(a), res.p_brut(a, c));
            else
                labels{c} = sprintf('%.1f%% (exclue)', 100 * Obs_a(c) / res.n_colorees(a));
            end
        end
        fig = figure('visible', 'off');
        ax = axes(fig);
        dessiner_camembert(ax, Obs_a, labels, couleurs_rgb);
        titre_obj = sgtitle(fig, sprintf(['Assemblée %d : %d cellules colorées | ' ...
            'p_{min-p} = %.4f | p_{chi2} = %.4f%s'], a, res.n_colorees(a), res.p_minp(a), ...
            res.p_chi2(a), suffixe), 'FontSize', 11, 'FontWeight', 'bold');
        titre_obj.HorizontalAlignment = 'left';
        exportgraphics(fig, strcat(namefull, ['piechart_noblack', num2str(a), '.png']), ...
            'Resolution', 150, 'ContentType', 'image');
        close(fig)
    end
end
end


function dessiner_contour(x, y, col, couleur_bord)
% Contour d'une cellule (enveloppe convexe des pixels du bord suite2p)
try
    kk = convhull(x, y);
catch
    kk = 1:numel(x);      % trop peu de points pour une enveloppe
end
patch(x(kk), y(kk), col, 'FaceAlpha', 1, 'EdgeColor', couleur_bord, 'LineWidth', 1)
end
