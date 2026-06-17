%% 2. Application de la matrice au volume Ex-vivo
clear; clc; close all;

% Paramètres
path_tif_exvivo = 'E:\Data\Aurelie\data\chroms\119\exvivo.tif'; % REMPLACER par le vrai chemin
radius_z = 10; % nombre de plans à extraire autour du plan moyen

% Chargement de la matrice et des profondeurs Z
load('registration_affine_2D.mat', 'tform', 'Z_ex');

% Calcul du plan Z de référence (plan moyen où se trouvent les points pointés)
z_ref = round(median(Z_ex));
fprintf('Plan Z de référence ex-vivo : %d\n', z_ref);

% Lecture des infos du TIF
info = imfinfo(path_tif_exvivo);
num_channels = 4;
total_pages = numel(info);
num_Z = total_pages / num_channels;

imgW = info(1).Width;
imgH = info(1).Height;

% Bornes Z à extraire
z_min = max(1, z_ref - radius_z);
z_max = min(num_Z, z_ref + radius_z);
nb_z_utiles = z_max - z_min + 1;

fprintf('Extraction du sous-volume de Z=%d à Z=%d (%d plans)\n', z_min, z_max, nb_z_utiles);

% Pré-allocation du sous-volume original
Volume_Sub_Exvivo = zeros(imgH, imgW, nb_z_utiles, num_channels, 'uint16');

% Lecture
idx_z = 0;
for z = z_min:z_max
    idx_z = idx_z + 1;
    for c = 1:num_channels
        page = (z - 1) * num_channels + c;
        Volume_Sub_Exvivo(:, :, idx_z, c) = imread(path_tif_exvivo, 'Index', page);
    end
end

% --- Déformation spatiale (Warping) vers le repère In-vivo ---
% On doit définir la taille de l'image de sortie (qui doit matcher le 512x512 de l'in-vivo)
outW = 512;
outH = 512;
Rout = imref2d([outH, outW]); % Référentiel cible

Volume_Warped = zeros(outH, outW, nb_z_utiles, num_channels, 'uint16');

fprintf('Warping du sous-volume vers 512x512...\n');
for z = 1:nb_z_utiles
    for c = 1:num_channels
        slice_in = Volume_Sub_Exvivo(:, :, z, c);
        % imwarp applique la transformation spatiale
        slice_out = imwarp(slice_in, tform, 'OutputView', Rout);
        Volume_Warped(:, :, z, c) = slice_out;
    end
end

disp('Warping terminé !');
% Volume_Warped contient maintenant ton ex-vivo aligné spatialement sur l'in-vivo
% (taille 512x512x21x4)

save('exvivo_warped_subvolume.mat', 'Volume_Warped', 'z_min', 'z_max', 'z_ref', '-v7.3');
