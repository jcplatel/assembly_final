%new raster 
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
%%
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
