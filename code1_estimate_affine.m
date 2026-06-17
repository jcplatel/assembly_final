%% 1. Calcul de la transformation affine (Estimation)
clear; clc; close all;

% Chemins vers vos CSV
% Le fichier in vivo ne contient que 2 axes XY pertinents (car Z=0 pour tous)
% Le fichier ex vivo contient les axes Z, Y, X (axis-0=Z, axis-1=Y, axis-2=X ou similaire)
% path_csv_invivo = "E:\Data\Aurelie\data\chroms\582\230331plane0\old\Pointage_Manuel_invivo.csv";
% path_csv_exvivo = "E:\Data\Aurelie\data\chroms\582\230331plane0\old\Pointage_Manuel_exvivo.csv";
path_csv_invivo = "E:\Data\Aurelie\data\chroms\175\221125plane0\Pointage_Manuel_invivo.csv";
path_csv_exvivo = "E:\Data\Aurelie\data\chroms\175\221125plane0\Pointage_Manuel_exvivo.csv";
[path,name,ext] = fileparts(path_csv_invivo);

% Lecture des tables
T_invivo = readtable(path_csv_invivo);
T_exvivo = readtable(path_csv_exvivo);

% On s'assure qu'on utilise les mêmes labels (les points correspondants)
[~, idx_in, idx_ex] = intersect(T_invivo.label, T_exvivo.label);

% Extraction des coordonnées
% Attention à l'ordre dans le CSV !
% Souvent napari/BiaPy sauve en (Z, Y, X) -> axis-0=Z, axis-1=Y, axis-2=X
Z_ex = T_exvivo.axis_0(idx_ex);
Y_ex = T_exvivo.axis_1(idx_ex);
X_ex = T_exvivo.axis_2(idx_ex);

% In vivo est 2D, donc Z=0 (axis-0), Y=axis-1, X=axis-2
Y_in = T_invivo.axis_1(idx_in);
X_in = T_invivo.axis_2(idx_in);

% Nos points de contrôle
movingPoints = [X_ex, Y_ex];  % Ex-vivo (à transformer)
fixedPoints  = [X_in, Y_in];  % In-vivo (référence cible)

% Calcul de la transformation affine
tform = fitgeotform2d(movingPoints, fixedPoints, "affine");

disp('Matrice Affine 2D estimée (Ex-vivo -> In-vivo) :');
disp(tform.A);

% Petit plot de vérification
figure;
scatter(fixedPoints(:,1), fixedPoints(:,2), 50, 'r', 'filled'); hold on;
% Application de la transfo aux points ex-vivo pour vérifier l'alignement
moving_transformed = transformPointsForward(tform, movingPoints);
scatter(moving_transformed(:,1), moving_transformed(:,2), 30, 'b', 'o', 'LineWidth', 2);
legend('In vivo (Cible)', 'Ex vivo (Transformé)');
title('Vérification de l''alignement affine 2D');
axis ij; % car ce sont des images
axis equal;

% Sauvegarde de la transformation
save(strcat(path,'\registration_affine_2D.mat'), 'tform', 'movingPoints', 'fixedPoints', 'Z_ex');
disp('Transformation sauvée dans registration_affine_2D.mat');
