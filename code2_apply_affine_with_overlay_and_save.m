%% 2. Application de la matrice au volume Ex-vivo, Superposition et Sauvegarde TIF
% clear; clc; close all;
% --- Paramètres ---
path_chroms = "E:\Data\Aurelie\data\chroms\175\221125plane0\";
path_tif_exvivo = "E:\Data\Aurelie\data\chroms\175\exvivo_rotated.tif"; % REMPLACER par le vrai chemin
path_tif_invivo = strcat(path_chroms,"invivo.tif"); % REMPLACER par le vrai chemin
radius_z = 10; % nombre de plans à extraire autour du plan moyen
canal_a_afficher = 2; % 2 pour le vert (GCaMP)
out_tif_volume = strcat(path_chroms,"exvivo_warped_subvolume.tif");
out_tif_overlay = strcat(path_chroms,"overlay_invivo_exvivo.tif");
out_tif_volume_sum2 = strcat(path_chroms,"exvivo_warped_subvolume_sum2.tif"); % NOUVEAU

% Chargement de la matrice et des profondeurs Z
% load('registration_affine_2D.mat', 'tform', 'Z_ex');
z_ref = round(median(Z_ex));
fprintf('Plan Z de référence ex-vivo : %d\n', z_ref);

info = imfinfo(path_tif_exvivo);
num_channels = 4;
total_pages = numel(info);
num_Z = total_pages / num_channels;
imgW = info(1).Width;
imgH = info(1).Height;

z_min = max(1, z_ref - radius_z);
z_max = min(num_Z, z_ref + radius_z);
nb_z_utiles = z_max - z_min + 1;

fprintf('Extraction du sous-volume de Z=%d à Z=%d (%d plans)\n', z_min, z_max, nb_z_utiles);

Volume_Sub_Exvivo = zeros(imgH, imgW, nb_z_utiles, num_channels, 'uint16');
idx_z = 0;
for z = z_min:z_max
    idx_z = idx_z + 1;
    for c = 1:num_channels
        page = (z - 1) * num_channels + c;
        Volume_Sub_Exvivo(:, :, idx_z, c) = imread(path_tif_exvivo, 'Index', page);
    end
end

% --- Déformation spatiale (Warping) vers le repère In-vivo ---
outW = 512;
outH = 512;
Rout = imref2d([outH, outW]); 
Volume_Warped = zeros(outH, outW, nb_z_utiles, num_channels, 'uint16');

fprintf('Warping du sous-volume vers 512x512...\n');
for z = 1:nb_z_utiles
    for c = 1:num_channels
        slice_in = Volume_Sub_Exvivo(:, :, z, c);
        Volume_Warped(:, :, z, c) = imwarp(slice_in, tform, 'OutputView', Rout);
    end
end
disp('Warping terminé !');

% --- Chargement de l'image In vivo ---
Img_Invivo = zeros(outH, outW, num_channels, 'uint16');
for c = 1:num_channels
    Img_Invivo(:, :, c) = imread(path_tif_invivo, 'Index', c);
end

% --- Projection MIP ---
Exvivo_MIP = squeeze(max(Volume_Warped, [], 3)); 

% --- Création et Sauvegarde de l'Overlay (RGB) ---
img_in_adj = imadjust(Img_Invivo(:,:,canal_a_afficher));
img_ex_adj = imadjust(Exvivo_MIP(:,:,canal_a_afficher));

% On crée une image RGB : Vert = in vivo, Magenta = ex vivo (Rouge + Bleu)
overlay_rgb = cat(3, img_ex_adj, img_in_adj, img_ex_adj);
imwrite(overlay_rgb, out_tif_overlay);
fprintf('Image superposée sauvée sous : %s\n', out_tif_overlay);

% --- Affichage ---
figure('Name', 'Superposition In vivo vs Ex vivo', 'Position', [100, 100, 1400, 500]);
subplot(1,3,1); imshow(img_in_adj); title('In vivo');
subplot(1,3,2); imshow(img_ex_adj); title('Ex vivo Warped (MIP)');
subplot(1,3,3); imshow(overlay_rgb); title('Overlay (Vert=In, Magenta=Ex)');

% --- Sauvegarde du Volume 3D en Multi-page TIF ---
fprintf('Sauvegarde du volume 3D aligné sous : %s\n', out_tif_volume);
if isfile(out_tif_volume)
    delete(out_tif_volume); % On supprime l'ancien pour ne pas append à l'infini
end

for z = 1:nb_z_utiles
    for c = 1:num_channels
        imwrite(Volume_Warped(:, :, z, c), out_tif_volume, 'WriteMode', 'append', 'Compression', 'none');
    end
end

disp('Volume 3D original sauvegardé !');
save('exvivo_warped_subvolume.mat', 'Volume_Warped', 'Exvivo_MIP', 'z_min', 'z_max', 'z_ref', '-v7.3');

%% === NOUVEAU : Volume avec somme de 2 plans consécutifs ===
fprintf('\n=== Création du volume avec somme de 2 plans consécutifs ===\n');

% Nombre de plans après somme de 2: si nb_z_utiles est impair, on garde le dernier plan seul
nb_z_sum2 = floor(nb_z_utiles / 2);
a_garde_dernier = (nb_z_utiles % 2) == 1;

% On crée un volume uint16 pour la somme (pas de overflow avec uint16 pour 2 plans)
Volume_Warped_Sum2 = zeros(outH, outW, nb_z_sum2 + a_garde_dernier, num_channels, 'uint16');

fprintf('Somme de 2 plans consécutifs : %d plans originaux → %d plans sommé2\n', ...
    nb_z_utiles, nb_z_sum2 + a_garde_dernier);

for z_idx = 1:nb_z_sum2
    z1 = 2 * z_idx - 1;  % plan 1, 3, 5, ...
    z2 = 2 * z_idx;      % plan 2, 4, 6, ...
    
    for c = 1:num_channels
        slice1 = Volume_Warped(:, :, z1, c);
        slice2 = Volume_Warped(:, :, z2, c);
        
        % Somme de 2 plans (uint16 → pas de overflow pour 2 plans)
        slice_sum = slice1 + slice2;
        
        Volume_Warped_Sum2(:, :, z_idx, c) = slice_sum;
    end
end

% Si nombre de plans impair, on garde le dernier plan seul
if a_garde_dernier
    z_last = nb_z_utiles;
    for c = 1:num_channels
        Volume_Warped_Sum2(:, :, nb_z_sum2 + 1, c) = Volume_Warped(:, :, z_last, c);
    end
end

disp('Volume sommé2 créé !');

% --- Sauvegarde du Volume 3D sommé2 en Multi-page TIF ---
fprintf('Sauvegarde du volume 3D sommé2 sous : %s\n', out_tif_volume_sum2);

if isfile(out_tif_volume_sum2)
    delete(out_tif_volume_sum2);
end

nb_z_final = nb_z_sum2 + a_garde_dernier;

for z = 1:nb_z_final
    for c = 1:num_channels
        imwrite(Volume_Warped_Sum2(:, :, z, c), out_tif_volume_sum2, 'WriteMode', 'append', 'Compression', 'none');
    end
end

disp('Volume 3D sommé2 sauvegardé !');

% Sauvegarde du MAT pour le volume sommé2
save('exvivo_warped_subvolume_sum2.mat', 'Volume_Warped_Sum2', 'nb_z_sum2', 'a_garde_dernier', '-v7.3');

fprintf('=== Tout est terminé et sauvegardé ! ===\n');
fprintf('- Volume original : %s\n', out_tif_volume);
fprintf('- Volume sommé2   : %s\n', out_tif_volume_sum2);
fprintf('- Overlay         : %s\n', out_tif_overlay);