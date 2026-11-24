function triangle_maxwell (colorcell,Donnees_Clustering,write_number,id_low_intensity)

if nargin < 4
    id_low_intensity = [];
end

X_tri = Donnees_Clustering(:,1);
Y_tri = Donnees_Clustering(:,2);

% % --- 2. Définition des labels ---
% labels_str = {'Red', 'Green', 'Blue', 'Yellow', 'Magenta', 'Cyan', 'White', 'Black'};
% n_classes = length(labels_str);
% all_classes = 1:8; 
%  
% same_class = (colorcell_method1(:) == colorcell_method2(:));
% different_class = (colorcell_method1(:) ~= colorcell_method2(:));
% dif_class=find (different_class==1)
% Palette pour method2
palette = [
    1.0 0.0 0.0;   % 1 red
    0.0 0.8 0.0;   % 2 green
    0.0 0.0 1.0;   % 3 blue
    0.9 0.9 0.0;   % 4 yellow
    1.0 0.0 1.0;   % 5 magenta
    0.0 0.8 0.8;   % 6 cyan
    0.92 0.92 0.92;% 7 white
    0.1 0.1 0.1    % 8 black
];
%si on donne id_low_intensity
if nargin == 4
    colorcell(id_low_intensity) = 8;
end

% Couleur de chaque cellule selon method2
C_points = palette(colorcell(:), :);


figure('Name', 'Triangle de Maxwell - couleur = method2', 'Color', 'w');
hold on;

% Triangle
plot([0 1 0.5 0], [0 0 sqrt(3)/2 0], 'k-', 'LineWidth', 1.5);

% Labels
text(1.02, 0, 'Red', 'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'left', 'Color', [0.8 0 0]);
text(-0.02, 0, 'Green', 'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'right', 'Color', [0 0.6 0]);
text(0.5, sqrt(3)/2 + 0.03, 'Blue', 'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'center', 'Color', [0 0 0.8]);

% Graduations
pas = 0.25;
for frac = pas:pas:1-pas
    text(frac, -0.04, sprintf('%.2f', frac), 'FontSize', 9, 'HorizontalAlignment', 'center', 'Color', [0.6 0 0]);
    x_droit = 1 - (frac/2);
    y_droit = frac * sqrt(3)/2;
    text(x_droit + 0.03, y_droit, sprintf('%.2f', frac), 'FontSize', 9, 'HorizontalAlignment', 'left', 'Color', [0 0 0.6]);
    x_gauche = 0.5 - ((1-frac)/2);
    y_gauche = frac * sqrt(3)/2;
    text(x_gauche - 0.03, y_gauche, sprintf('%.2f', 1-frac), 'FontSize', 9, 'HorizontalAlignment', 'right', 'Color', [0 0.5 0]);
end

% % 1) Tous les points coloriés avec method2
% scatter(X_tri, Y_tri, 35, C_points, ...
%     'filled', 'MarkerEdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.5);

% % % 2) Les points où method1 differents method2 sont marqués en blanc par-dessus
% % scatter(X_tri(different_class), Y_tri(different_class), 55, [1 1 1],  'MarkerEdgeColor', 'k', 'LineWidth', 1.0);
% scatter(X_tri(different_class), Y_tri(different_class), 35, C_points(different_class,:),'filled','LineWidth', 1.0);
% % 1) Tous les points coloriés avec method2 - SAUVEGARDE DU HANDLE hScatter
hScatter = scatter(X_tri, Y_tri, 35, C_points, 'filled', 'MarkerEdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.5);

% =========================================================================
% AJOUT DES INFOBULLES (DATATIPS) AU SURVOL
% 1. Calcul inverse des composantes R, V, B à partir de X et Y
B_coord = Y_tri / (sqrt(3)/2);
R_coord = X_tri - 0.5 .* B_coord;
V_coord = 1 - R_coord - B_coord;

% 2. Ajout des lignes dans l'infobulle
hScatter.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('R', R_coord, '%.3f');
hScatter.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('V', V_coord, '%.3f');
hScatter.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('B', B_coord, '%.3f');

% (Optionnel) Ajout du numéro du point (Index)
idx = 1:length(X_tri);
hScatter.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Index', idx, '%d');



if write_number==1
    for k = 1:length(X_tri)
        text(X_tri(k) + 0.015, Y_tri(k), num2str(k), ...
        'FontSize', 7, 'Color', [0.4 0.4 0.4]); % Gris pour ne pas masquer les couleurs
    end
end

axis equal;
axis off;
hold off;

%%%

end