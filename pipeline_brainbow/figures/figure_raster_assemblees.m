function figure_raster_assemblees(A, Race, namefull)
% = rastercolor : raster des SCE coloré par ASSEMBLÉE (pas par couleur
% brainbow : cette figure appartient donc à l'étape 1).
% Mêmes fichiers qu'avant : rasterall, rasterall_blanc, rasterall_blanc_core,
% rasterall_blanc_thick, raster_assemblies_only, rasterallblack.
% Changements :
%  - autant de couleurs que d'assemblées, jamais une couleur brainbow
%    (couleurs_assemblees, comme l'interface)
%  - chaque figure est fermée (elles s'accumulaient d'un K à l'autre)
%  - plus de « return » au milieu du script
%  - images construites sans double boucle

NCl = A.NCl;
if NCl < 1, return; end
assembly  = A.assemblyortho;
x1        = A.x1;
RList     = A.RList;
PCl_final = A.PCl_final;
NCell     = size(Race, 1);

% Couleurs des assemblées : les mêmes que dans l'interface, entre les
% couleurs brainbow (utilitaires/couleurs_assemblees.m) ; plus foncées sur
% fond blanc
colors_noir = couleurs_assemblees(NCl, 'noir');
colors      = couleurs_assemblees(NCl, 'blanc');

% Assemblée de chaque cellule (NCl + 1 = hors assemblée), dans l'ordre x1
cell_id = (NCl + 1) * ones(1, NCell);
for n = 1:NCl
    cell_id(assembly{n}) = n;
end
cell_id_sorted = cell_id(x1);
Race_ordered = Race(x1, RList);
actif = Race_ordered > 0;

% 1. Fond noir, cellules hors assemblée en blanc
img = image_raster(actif, [colors_noir; 1 1 1], cell_id_sorted, [0 0 0]);
exporter_raster(img, strcat(namefull, 'rasterall.png'), 'noir', 150, 'sorted SCE #');

% 2. Fond blanc, cellules hors assemblée en gris clair
couleurs_gris = [colors; 0.75 0.75 0.75];
img = image_raster(actif, couleurs_gris, cell_id_sorted, [1 1 1]);
exporter_raster(img, strcat(namefull, 'rasterall_blanc.png'), 'blanc', 300, 'sorted SCE #');

% 3. « core » : couleur seulement quand l'assemblée de la cellule est
%    significativement active dans ce SCE (PCl_final), gris sinon
PCl_ordered = PCl_final(:, RList);
id_point = repmat(cell_id_sorted(:), 1, size(actif, 2));
dans_asm = cell_id_sorted <= NCl;
asm_active = false(size(actif));
asm_active(dans_asm, :) = PCl_ordered(cell_id_sorted(dans_asm), :) > 0;
id_point(~asm_active) = NCl + 1;
img = epaissir(image_raster(actif, couleurs_gris, id_point, [1 1 1]));
exporter_raster(img, strcat(namefull, 'rasterall_blanc_core.png'), 'blanc', 300, 'sorted SCE #');

% 4. Toutes les décharges colorées par assemblée, points épaissis
img = epaissir(image_raster(actif, couleurs_gris, cell_id_sorted, [1 1 1]));
exporter_raster(img, strcat(namefull, 'rasterall_blanc_thick.png'), 'blanc', 300, 'sorted SCE #');

% 5. Seulement les SCE où au moins une assemblée est active
sce_asm = sum(PCl_ordered, 1) > 0;
if any(sce_asm)
    img = epaissir(image_raster(actif(:, sce_asm), couleurs_gris, cell_id_sorted, [1 1 1]));
    exporter_raster(img, strcat(namefull, 'raster_assemblies_only.png'), 'blanc', 300, ' SCE # ');
end

% 6. Raster noir et blanc
fig = figure('visible', 'off', 'Color', 'k');
ax = axes('Parent', fig);
imagesc(ax, actif);
colormap(ax, gray);
axis(ax, 'ij'); axis(ax, 'tight');
set(ax, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
xlabel(ax, 'sorted SCE #'); ylabel(ax, 'sorted Cell #');
exportgraphics(fig, strcat(namefull, 'rasterallblack.png'), 'Resolution', 150, 'BackgroundColor', 'black');
close(fig)
end


function img = image_raster(actif, palette, id, fond)
% Image RGB : chaque point actif prend la couleur palette(id) ; les autres
% points prennent la couleur de fond. id : une valeur par ligne, ou par point.
if ~isequal(size(id), size(actif))
    id = repmat(id(:), 1, size(actif, 2));
end
img = zeros([size(actif) 3]);
for c = 1:3
    canal = reshape(palette(id, c), size(actif));
    canal(~actif) = fond(c);
    img(:, :, c) = canal;
end
end


function img_thick = epaissir(img)
% Dilatation 2x2 des points colorés sur fond blanc (meilleure lisibilité)
se = strel('rectangle', [2 2]);
img_thick = img;
for c = 1:3
    img_thick(:, :, c) = 1 - imdilate(1 - img(:, :, c), se);
end
end


function exporter_raster(img, fichier, fond, resolution, legende_x)
if strcmp(fond, 'noir'), c_fig = 'k'; c_axes = 'w'; c_export = 'black';
else, c_fig = 'w'; c_axes = 'k'; c_export = 'white'; end
fig = figure('visible', 'off', 'Color', c_fig);
imshow(img);
axis on;
xlabel(legende_x);
ylabel('sorted Cell #');
set(gca, 'XColor', c_axes, 'YColor', c_axes, 'FontSize', 12);
box off
exportgraphics(fig, fichier, 'Resolution', resolution, 'BackgroundColor', c_export, 'ContentType', 'image');
close(fig)
end
