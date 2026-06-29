function plot_bicluster_race_colored(Race, results_bc, namefull)

cellLabels = results_bc.cellLabels;   % [NCell × 1]
sceLabels  = results_bc.sceLabels;    % [NSCE × 1]

[NCell, NSCE] = size(Race);
nCellClusters = max(cellLabels);

% 1. Tri des cellules et des SCE par cluster
[~, idxCells_sorted] = sort(cellLabels);   % ordre des cellules (y)
[~, idxSCE_sorted]   = sort(sceLabels);    % ordre des SCE (x)

Race_sorted        = Race(idxCells_sorted, idxSCE_sorted);
cellLabels_sorted  = cellLabels(idxCells_sorted);
sceLabels_sorted   = sceLabels(idxSCE_sorted); %#ok<NASGU> % si tu veux t'en servir plus tard

% 2. Construire une matrice "couleur" :
%    - 0 = fond (pas d'activité)
%    - c = cluster de cellule c pour les 1 de Race
Race_color = zeros(size(Race_sorted));
for c = 1:nCellClusters
    mask_cells = (cellLabels_sorted == c);           % lignes de ce cluster
    % mettre la valeur c là où Race_sorted == 1 pour ces cellules
    Race_color(mask_cells, :) = Race_color(mask_cells, :) + ...
        (Race_sorted(mask_cells, :) == 1) * c;
end

% 3. Plot unique : cellules triées (y), SCE triés (x), couleur = assemblée
figure('visible','on');

imagesc(Race_color);
axis image;
xlabel('SCE triés');
ylabel('Cellules triées');
title('Assemblées (code couleur)');

% Colormap : fond gris pour 0, puis couleurs pour les assemblées
cmap = [0.8 0.8 0.8; lines(nCellClusters)];  % 1ère ligne = gris
colormap(cmap);
caxis([0 nCellClusters]);    % 0 = fond, 1..nCellClusters = assemblées

colorbar('Ticks', 1:nCellClusters, 'TickLabels', ...
    arrayfun(@(c) sprintf('Assembly %d', c), 1:nCellClusters, 'UniformOutput', false));

if exist('namefull','var') && ~isempty(namefull)
    exportgraphics(gcf, fullfile(namefull, 'bicluster_race_color.png'), ...
                   'Resolution', 300);
end

% close gcf;

end