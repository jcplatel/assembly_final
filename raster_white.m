%new raster (fond BLANC, version scatter - points visibles)
assembly=assemblyortho;
% assembly=assemblyraw;
colors = [
    1 0 0;     % Red 
    0 1 0;     % Green
    0 0 1;     % Blue
    1 1 0;     % Yellow
    1 0 1;     % Magenta
    0 1 1;     % Cyan
    0.5 0 0;   % Maroon
    0 0.5 0;   % Olive
    0 0 0.5;   % Navy
    0.5 0.5 0; % Olive Green
    0.5 0 0.5; % Purple
    0 0.5 0.5; % Teal
    0.75 0 0;  % Dark Red
    0 0.75 0;  % Lime Green
    0 0 0.75;  % Dark Blue
    0.75 0.75 0; % Dark Yellow
    0.75 0 0.75; % Dark Magenta
    0 0.75 0.75; % Dark Cyan
    ];
% Cl_order = [ Cl1  Cl2 Cl3 Cl4];
Cl_order = [Cl0 Cl1 Cl2 Cl3];
Race_ordered=Race(x1,Cl_order);

% Couleur des cellules NON assignées à une assemblée :
% on utilise un GRIS FONCE (pas clair !) pour rester bien visible sur fond blanc
colors(NCl+1,:) = [0.4 0.4 0.4]; % Gris foncé au lieu de gris clair

cell_id=(NCl+1)*ones(1,NCell);
for n=1:NCl
    cell_id(1,cell2mat(assembly(n))) = n;
end
cell_id=cell_id(x1);

%% Construction du scatter plot (points de taille controlable)
[row_idx, col_idx] = find(Race_ordered > 0); % row = cellule, col = SCE
point_colors = colors(cell_id(row_idx), :);   % couleur associée à chaque point

fig = figure('visible','off');
set(fig, 'Color', 'w'); % Fond de la figure en BLANC
ax = axes('Parent', fig);
set(ax, 'Color', 'w');
hold(ax, 'on');

% MarkerSize : à ajuster si besoin (4-8 selon la densité de points souhaitée)
markerSize = 6;
scatter(ax, col_idx, row_idx, markerSize, point_colors, 'filled');

set(ax, 'YDir', 'reverse'); % Pour garder cellule #1 en haut, comme imshow
xlim(ax, [1, size(Race_ordered,2)]);
ylim(ax, [1, NCell]);
xlabel(ax, 'sorted SCE #');
ylabel(ax, 'sorted Cell #');
set(ax, 'XColor', 'k', 'YColor', 'k'); % Axes en noir, lisibles sur fond blanc
box(ax, 'on');

namegraph=strcat(namefull,['rasterall_blanc_v2' , '.png']);
if isfolder(namefull)
    exportgraphics(gcf,namegraph,'Resolution',150,'BackgroundColor','white', 'ContentType', 'image')
    close gcf
end