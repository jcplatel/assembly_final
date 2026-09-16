%new raster 
assembly=assemblyortho;
% assembly=assemblyraw;

colors = [
    1 0 0;     % Red 
    0 1 0;     % Green
    0 0 1;     % Blue
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

% fig = figure('visible','on');
% % fig=figure;
% set(fig, 'Color', 'k'); % Set the figure background color to black
% 
% ax = axes('Parent', fig);
% set(ax, 'Color', 'k'); % Set the axis background color to black
colors(NCl+1,:)=1;

cell_id=(NCl+1)*ones(1,NCell);

for n=1:NCl
    cell_id(1,cell2mat(assembly(n))) = n;
end
cell_id=cell_id(x1);

%%
% Construire une matrice de couleurs RGB selon Race_ordered et cell_id
img = zeros(NCell, size(Race_ordered,2), 3);

for n = 1:NCell
    for y = 1:size(Race_ordered,2)
        if Race_ordered(n,y) > 0
            img(n,y,:) = colors(cell_id(n),:);
        else
            img(n,y,:) = [0 0 0]; % Noir pour les zéros
        end
    end
end
fig = figure('visible','off');
set(fig, 'Color', 'k');
% figure('Color','k');
imshow(img);
axis on;
xlabel('sorted SCE #');
ylabel('sorted Cell #');
set(gca,'XColor','w','YColor','w');

namegraph=strcat(namefull,['rasterall' , '.png']);

if isfolder(namefull)
    exportgraphics(gcf,namegraph,'Resolution',150,'BackgroundColor','black', 'ContentType', 'image')

    close gcf
end

%%%%


% Couleur des cellules NON assignées à une assemblée :
% sur fond blanc, on ne peut plus utiliser du blanc [1 1 1] (invisible).
% On utilise un gris clair pour rester visible tout en restant discret.
colors(NCl+1,:) = [0.75 0.75 0.75]; % Gris clair au lieu de blanc

cell_id=(NCl+1)*ones(1,NCell);
for n=1:NCl
    cell_id(1,cell2mat(assembly(n))) = n;
end
cell_id=cell_id(x1);

%%
% Construire une matrice de couleurs RGB selon Race_ordered et cell_id
img = ones(NCell, size(Race_ordered,2), 3); % <-- ones() au lieu de zeros() : fond blanc par défaut
for n = 1:NCell
    for y = 1:size(Race_ordered,2)
        if Race_ordered(n,y) > 0
            img(n,y,:) = colors(cell_id(n),:);
        else
            img(n,y,:) = [1 1 1]; % Blanc pour les zéros (au lieu de noir)
        end
    end
end

fig = figure('visible','off');
set(fig, 'Color', 'w'); % Fond de la figure en BLANC (au lieu de 'k')
imshow(img);
axis on;
xlabel('sorted SCE #');
ylabel('sorted Cell #');
set(gca,'XColor','k','YColor','k'); % Texte des axes en NOIR (au lieu de 'w') pour rester lisible sur fond blanc
box off
namegraph=strcat(namefull,['rasterall_blanc' , '.png']);
if isfolder(namefull)
    exportgraphics(gcf,namegraph,'Resolution',300,'BackgroundColor','white', 'ContentType', 'image')
    % close gcf
end

%%%nouveau code
%% 1. Préparation des matrices et assignations
NSCE_ord = size(Race_ordered, 2);

% Couleur grise pour le bruit de fond / décharges hors assemblée
color_noise = [0.75, 0.75, 0.75]; 

% Assignation des cellules aux clusters
cell_id = (NCl + 1) * ones(1, NCell);
for n = 1:NCl
    cell_id(1, cell2mat(assembly(n))) = n;
end
cell_id_sorted = cell_id(x1); % Alignement sur le tri vertical des cellules

% Alignement horizontal de la matrice d'activation des assemblées
% RList contient l'ordre des colonnes utilisé pour trier Race_ordered
PCl_ordered = PCl_final(:, RList); 

%% 2. Construction de l'image RGB (Fond Blanc)
img = ones(NCell, NSCE_ord, 3); % Fond blanc par défaut

for n = 1:NCell
    cid = cell_id_sorted(n); % ID du cluster pour la cellule n

    for y = 1:NSCE_ord
        if Race_ordered(n, y) > 0
            % Condition clé : la cellule appartient à une assemblée (cid <= NCl)
            % ET cette assemblée est active lors de ce SCE (PCl_ordered == 1)
            if (cid <= NCl) && (PCl_ordered(cid, y) == 1)
                img(n, y, :) = colors(cid, :); % Couleur de l'assemblée
            else
                img(n, y, :) = color_noise;    % Gris clair (bruit / hors assemblée)
            end
        end
    end
end

%% 3. Épaississement des points pour améliorer la lisibilité
% Dilatation morphologique (3x3 pixels) pour grossir chaque point
% Réduit le fond blanc pour donner plus d'épaisseur aux couleurs
se = strel("rectangle",[2 2]); % Épaisseur : testez 'square', 2 ou 'disk', 1 selon le rendu
img_thick = img;
for c = 1:3
    % On dilate les pixels sombres/colorés par rapport au fond blanc (valeur 1)
    img_thick(:, :, c) = 1 - imdilate(1 - img(:, :, c), se);
end

%% 4. Affichage et Exportation
fig = figure('visible', 'off');
set(fig, 'Color', 'w');

imshow(img_thick);
axis on;
xlabel('sorted SCE #');
ylabel('sorted Cell #');
set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', 12);
box off;

namegraph = strcat(namefull, 'rasterall_blanc_core.png');
if isfolder(namefull)
    exportgraphics(gcf, namegraph, 'Resolution', 300, 'BackgroundColor', 'white', 'ContentType', 'image');
    % close gcf;
end

%%%
%% 1. Définition des couleurs et assignation des cellules
% Couleur pour les cellules non assignées (Gris clair discret sur fond blanc)
colors(NCl + 1, :) = [0.75, 0.75, 0.75];

% Attribution de chaque cellule à son assemblée
cell_id = (NCl + 1) * ones(1, NCell);
for n = 1:NCl
    cell_id(1, cell2mat(assembly(n))) = n;
end
cell_id_sorted = cell_id(x1); % Réordonné selon le tri vertical

%% 2. Construction de la matrice RGB de base (Fond Blanc)
NSCE_ord = size(Race_ordered, 2);
img = ones(NCell, NSCE_ord, 3); % Fond blanc par défaut

for n = 1:NCell
    cid = cell_id_sorted(n);
    cell_color = colors(cid, :);

    for y = 1:NSCE_ord
        if Race_ordered(n, y) > 0
            img(n, y, :) = cell_color; % Colore toutes les décharges de la cellule
        end
    end
end

%% 3. Épaississement des points (Dilatation morphologique)
% 'disk', 1 donne un épaississement doux (croix de 3x3)
% Pour des points encore plus marqués, testez strel('square', 2) ou strel('disk', 2)
% se = strel('line', 2,90);
se = strel("rectangle",[2 2]);

img_thick = img;
for c = 1:3
    % Inversion (1 - img) pour dilater les points actifs sur le fond blanc
    img_thick(:, :, c) = 1 - imdilate(1 - img(:, :, c), se);
end

%% 4. Affichage et Exportation Haute Définition
fig = figure('visible', 'off');
set(fig, 'Color', 'w');

imshow(img_thick);
axis on;
xlabel('sorted SCE #');
ylabel('sorted Cell #');
set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', 12);
box off;

namegraph = strcat(namefull, 'rasterall_blanc_thick.png');
if isfolder(namefull)
    exportgraphics(gcf, namegraph, 'Resolution', 300, 'BackgroundColor', 'white', 'ContentType', 'image');
    % close gcf;
end


%%%
%% 1. Sélection et tri exclusif des SCE appartenant à des assemblées
% On ne garde que les SCE associés à au moins une assemblée (on exclut Cl0)
RList_asm = [Cl1, Cl2, Cl3]; 
NSCE_asm  = length(RList_asm);

if NSCE_asm == 0
    warning('Aucun SCE n''appartient à une assemblée.');
    return;
end

% Matrice Race sous-échantillonnée et réordonnée (Cellules triées x SCE assemblées)
Race_asm_ordered = Race(x1, RList_asm);

%% 2. Définition des couleurs et assignation des cellules
% Couleur pour les cellules non assignées (Gris clair discret)
colors(NCl + 1, :) = [0.75, 0.75, 0.75];

% Attribution de chaque cellule à son assemblée
cell_id = (NCl + 1) * ones(1, NCell);
for n = 1:NCl
    cell_id(1, cell2mat(assembly(n))) = n;
end
cell_id_sorted = cell_id(x1); % Alignement sur le tri vertical (x1)

%% 3. Construction de l'image RGB (Fond Blanc)
img = ones(NCell, NSCE_asm, 3); % Matrice fond blanc de dimension [NCell x NSCE_asm]

for n = 1:NCell
    cid = cell_id_sorted(n);
    cell_color = colors(cid, :);

    for y = 1:NSCE_asm
        if Race_asm_ordered(n, y) > 0
            img(n, y, :) = cell_color; % Colore toutes les décharges de la cellule
        end
    end
end

%% 4. Épaississement des points (Dilatation morphologique)
% se = strel('disk', 1); % Épaisseur : ajustez à strel('disk', 2) si besoin
se = strel("rectangle",[2 2]); 
img_thick = img;
for c = 1:3
    % Dilatation des points colorés par rapport au fond blanc
    img_thick(:, :, c) = 1 - imdilate(1 - img(:, :, c), se);
end

%% 5. Affichage et Exportation
fig = figure('visible', 'off');
set(fig, 'Color', 'w');

imshow(img_thick);
axis on;
xlabel(' SCE # ');
ylabel('sorted Cell #');
set(gca, 'XColor', 'k', 'YColor', 'k', 'FontSize', 12);
box off;

namegraph = strcat(namefull, 'raster_assemblies_only.png');
if isfolder(namefull)
    exportgraphics(gcf, namegraph, 'Resolution', 300, 'BackgroundColor', 'white', 'ContentType', 'image');
    % close gcf;
end

%%%
%new raster 
% assembly=assemblyortho;
% 
% colors = [
%     1 0 0;     % Red 
%     0 1 0;     % Green
%     0 0 1;     % Blue
%     1 1 0;     % Yellow
%     1 0 1;     % Magenta
%     0 1 1;     % Cyan
%     0.5 0 0;   % Maroon
%     0 0.5 0;   % Olive
%     0 0 0.5;   % Navy
%     0.5 0.5 0; % Olive Green
%     0.5 0 0.5; % Purple
%     0 0.5 0.5; % Teal
%     0.75 0 0;  % Dark Red
%     0 0.75 0;  % Lime Green
%     0 0 0.75;  % Dark Blue
%     0.75 0.75 0; % Dark Yellow
%     0.75 0 0.75; % Dark Magenta
%     0 0.75 0.75; % Dark Cyan
% ];
% Cl_order = [ Cl1  Cl2 Cl3 Cl4];
% % Cl_order = [Cl0 Cl1 Cl2 Cl3 Cl4];
% Race_ordered=Race(x1,Cl_order);

% fig = figure('visible','off');
% % fig=figure;
% set(fig, 'Color', 'k'); % Set the figure background color to black
% 
% ax = axes('Parent', fig);
% set(ax, 'Color', 'k'); % Set the axis background color to black
% colors(NCl+1,:)=1;
% 
% cell_id=(NCl+1)*ones(1,NCell);
% 
% for n=1:NCl
%     cell_id(1,cell2mat(assembly(n))) = n;
% end
% cell_id=cell_id(x1);
% 
% hold on
% for n = 1:NCell
%     for y=1:size(Race_ordered,2)
%         if Race_ordered(n,y)>0
%             rectangle('Position',[y,n, 1,1],'FaceColor','w','EdgeColor','w','LineWidth',0.5)
%         end
%     end
% end
% 
% axis ij
% axis tight
% set(ax, 'XColor', 'w');
% set(ax, 'YColor', 'w');
% set(ax, 'Color', 'k');
% xlabel('sorted SCE #')     %was RACE
% ylabel('sorted Cell #')
% 
% hold off
% namegraph=strcat(namefull,['rasterallblack' , '.png']);
% if isfolder(namefull)
%     exportgraphics(gcf,namegraph,'Resolution',150,'BackgroundColor','black')
%     close gcf
% end
mask = Race_ordered > 0;

fig = figure('visible','off');
set(fig, 'Color', 'k');
ax = axes('Parent', fig);
set(ax, 'Color', 'k');

% Affichage avec imagesc, blanc = 1, noir = 0
imagesc(mask);
colormap(gray); % Couleurs noir/blanc
axis ij;
axis tight;
set(ax, 'XColor', 'w', 'YColor', 'w');

xlabel('sorted SCE #');
ylabel('sorted Cell #');
hold off;

namegraph = fullfile(namefull, 'rasterallblack.png');
if isfolder(namefull)
    exportgraphics(fig, namegraph, 'Resolution', 150, 'BackgroundColor', 'black');
    close(fig);
end
