function plot_bicluster_raster(Race, results_bc, namefull)

cellLabels = results_bc.cellLabels;   % [NCell × 1]
sceLabels  = results_bc.sceLabels;    % [NSCE × 1]

[NCell, NSCE] = size(Race);
nCellClusters = max(cellLabels);

% 1. Tri des cellules et des SCE par cluster
[~, idxCells_sorted] = sort(cellLabels);   % ordre des cellules (y)
[~, idxSCE_sorted]   = sort(sceLabels);    % ordre des SCE (x)

Race_sorted       = Race(idxCells_sorted, idxSCE_sorted);
cellLabels_sorted = cellLabels(idxCells_sorted);

% 2. Préparation de la figure type "rasterall"
figure('visible','on');
set(gcf, 'Color', 'k');   % fond de figure noir
axes('Color', 'k');       % fond des axes noir
hold on;

% palette de couleurs pour les assemblées
colors = lines(nCellClusters);   % n couleurs différentes

% 3. Points blancs (toutes les activations, comme "fond")
[row_all, col_all] = find(Race_sorted > 0);
plot(col_all, row_all, '.', 'Color', [0.9 0.9 0.9], 'MarkerSize', 4);

% 4. Par-dessus, points colorés par assemblée (cluster de cellules)
for c = 1:nCellClusters
    mask_cells = (cellLabels_sorted == c);      % lignes appartenant à l'assemblée c
    
    % indices (ligne, colonne) où Race_sorted == 1 ET cellule dans cluster c
    [row_c, col_c] = find(Race_sorted(mask_cells, :) > 0);
    
    % Attention : row_c est relatif au sous-ensemble mask_cells, on le remet
    % dans les coordonnées globales (ligne réelle)
    cell_indices_c = find(mask_cells);
    row_c_global = cell_indices_c(row_c);
    
    plot(col_c, row_c_global, '.', 'Color', colors(c,:), 'MarkerSize', 6);
end

% 5. Mise en forme
set(gca, 'YDir', 'normal');   % 1 en bas, NCell en haut (comme ton raster)
xlabel('sorted SCE #', 'Color', 'w');
ylabel('sorted Cell #', 'Color', 'w');
title('sorted Cell × SCE (biclustering)', 'Color', 'w');

xlim([1 NSCE]);
ylim([1 NCell]);

% ticks en blanc
set(gca, 'XColor', 'w', 'YColor', 'w');

if exist('namefull','var') && ~isempty(namefull)
    exportgraphics(gcf, fullfile(namefull, 'bicluster_rasterall.png'), ...
                   'Resolution', 300);
end

hold off;
% close gcf;

end