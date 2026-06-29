%% 2. Application de la matrice au volume Ex-vivo, Superposition et Sauvegarde TIF
% clear; clc; close all;

% --- Paramètres ---
path_chroms = "E:\Data\Aurelie\data\chroms\175\221125plane0\";
path_tif_exvivo = "E:\Data\Aurelie\data\chroms\175\exvivo_rotated_dimcorrected.tif"; 
path_tif_invivo = strcat(path_chroms,"invivo.tif"); 

radius_z = 14; % nombre de plans à extraire autour du plan moyen
canal_a_afficher = 2; % 2 pour le vert (GCaMP)
out_tif_volume = strcat(path_chroms, "exvivo_warped_subvolumebis.tif");
out_tif_overlay = strcat(path_chroms, "overlay_invivo_exvivo.tif");
out_tif_volume_sum2 = strcat(path_chroms, "exvivo_warped_subvolume_sum2bis.tif");

% Assurez-vous d'avoir load('registration_affine_2D.mat', 'tform', 'Z_ex') 
% dans le workspace avant de lancer ce code !
z_ref = round(median(Z_ex));
fprintf('Plan Z de référence ex-vivo : %d\n', z_ref);

info = imfinfo(path_tif_exvivo);
total_pages = numel(info);
num_channels = 4;
num_Z = total_pages / num_channels;

imgW = info(1).Width;
imgH = info(1).Height;
z_min = max(1, z_ref - radius_z);
z_max = min(num_Z, z_ref + radius_z);
nb_z_utiles = z_max - z_min + 1;
fprintf('Extraction du sous-volume de Z=%d à Z=%d (%d plans)\n', z_min, z_max, nb_z_utiles);

% --- 1. DÉTECTION DU FORMAT (8-bits ou 16-bits) ---
temp_img = imread(path_tif_exvivo, 'Index', 1);
img_class = class(temp_img); % Détectera 'uint8' ou 'uint16'
fprintf('Format de l''image détecté : %s\n', img_class);

Volume_Sub_Exvivo = zeros(imgH, imgW, nb_z_utiles, num_channels, img_class);

% --- 2. DÉTECTION DE L'ORDRE DES PAGES (C/Z ou Z/C) ---
if isfield(info(1), 'ImageDescription')
    image_desc = info(1).ImageDescription;
    est_ordre_CZ = contains(image_desc, 'order=cz', 'IgnoreCase', true);
else
    est_ordre_CZ = false;
end

idx_z = 0;
for z = z_min:z_max
    idx_z = idx_z + 1;
    for c = 1:num_channels
        if est_ordre_CZ
            % Ordre ImageJ (C1Z1, C1Z2... puis C2Z1...)
            page = (c - 1) * num_Z + z;
        else
            % Ordre Classique (Z1C1, Z1C2... puis Z2C1...)
            page = (z - 1) * num_channels + c;
        end
        Volume_Sub_Exvivo(:, :, idx_z, c) = imread(path_tif_exvivo, 'Index', page);
    end
end

% --- Déformation spatiale (Warping) vers le repère In-vivo ---
outW = 512;
outH = 512;
Rout = imref2d([outH, outW]); 
Volume_Warped = zeros(outH, outW, nb_z_utiles, num_channels, img_class);
fprintf('Warping du sous-volume vers 512x512...\n');
for z = 1:nb_z_utiles
    for c = 1:num_channels
        slice_in = Volume_Sub_Exvivo(:, :, z, c);
        Volume_Warped(:, :, z, c) = imwarp(slice_in, tform, 'OutputView', Rout);
    end
end
disp('Warping terminé !');

% --- Chargement de l'image In vivo ---
Img_Invivo = zeros(outH, outW, num_channels, img_class);
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


%% --- Sauvegarde du Volume 3D en Multi-page TIF (Format ImageJ) avec IMWRITE ---
fprintf('Sauvegarde du volume 3D aligné sous : %s\n', out_tif_volume);
if isfile(out_tif_volume)
    delete(out_tif_volume); % On supprime l'ancien
end

total_images = nb_z_utiles * num_channels;
ij_desc = sprintf('ImageJ=1.53c\nimages=%d\nchannels=%d\nslices=%d\nframes=1\nhyperstack=true\nmode=composite\nloop=false', ...
    total_images, num_channels, nb_z_utiles);

page = 1;
for z = 1:nb_z_utiles
    for c = 1:num_channels
        img_slice = Volume_Warped(:, :, z, c);
        
        if page == 1
            % Le premier plan crée le fichier et inclut les tags ImageJ
            imwrite(img_slice, out_tif_volume, 'tif', ...
                'Description', ij_desc, ...
                'WriteMode', 'overwrite');
        else
            % Les plans suivants sont ajoutés à la suite (pas de décalage X possible)
            imwrite(img_slice, out_tif_volume, 'tif', ...
                'WriteMode', 'append');
        end
        page = page + 1;
    end
end
disp('Volume 3D original sauvegardé avec imwrite !');
save('exvivo_warped_subvolume.mat', 'Volume_Warped', 'Exvivo_MIP', 'z_min', 'z_max', 'z_ref', '-v7.3');


%% === NOUVEAU : Volume avec somme de 2 plans consécutifs ===
fprintf('\n=== Création du volume avec somme de 2 plans consécutifs ===\n');
nb_z_sum2 = floor(nb_z_utiles / 2);
a_garde_dernier = mod(nb_z_utiles, 2) == 1;
nb_z_final = nb_z_sum2 + a_garde_dernier;

Volume_Warped_Sum2 = zeros(outH, outW, nb_z_final, num_channels, img_class);
fprintf('Somme de 2 plans consécutifs : %d plans originaux → %d plans sommé2\n', ...
    nb_z_utiles, nb_z_final);

for z_idx = 1:nb_z_sum2
    z1 = 2 * z_idx - 1;  % plan 1, 3, 5, ...
    z2 = 2 * z_idx;      % plan 2, 4, 6, ...
    
    for c = 1:num_channels
        slice1 = Volume_Warped(:, :, z1, c);
        slice2 = Volume_Warped(:, :, z2, c);
        
        % --- CORRECTION SATURATION : Somme en double ---
        slice_sum = double(slice1) + double(slice2);
        
        % On rabat proprement selon le format détecté
        if strcmp(img_class, 'uint8')
            Volume_Warped_Sum2(:, :, z_idx, c) = uint8(min(slice_sum, 255));
        else
            Volume_Warped_Sum2(:, :, z_idx, c) = uint16(min(slice_sum, 65535));
        end
    end
end

% Si nombre de plans impair, on garde le dernier plan seul
if a_garde_dernier
    z_last = nb_z_utiles;
    for c = 1:num_channels
        Volume_Warped_Sum2(:, :, nb_z_final, c) = Volume_Warped(:, :, z_last, c);
    end
end
disp('Volume sommé2 créé !');

% --- Sauvegarde du Volume 3D sommé2 en Multi-page TIF (Format ImageJ) avec IMWRITE ---
fprintf('Sauvegarde du volume 3D sommé2 sous : %s\n', out_tif_volume_sum2);
if isfile(out_tif_volume_sum2)
    delete(out_tif_volume_sum2);
end

total_images_sum2 = nb_z_final * num_channels;
ij_desc_sum2 = sprintf('ImageJ=1.53c\nimages=%d\nchannels=%d\nslices=%d\nframes=1\nhyperstack=true\nmode=composite\nloop=false', ...
    total_images_sum2, num_channels, nb_z_final);

page = 1;
for z = 1:nb_z_final
    for c = 1:num_channels
        img_slice2 = Volume_Warped_Sum2(:, :, z, c);
        
        if page == 1
            imwrite(img_slice2, out_tif_volume_sum2, 'tif', ...
                'Description', ij_desc_sum2, ...
                'WriteMode', 'overwrite');
        else
            imwrite(img_slice2, out_tif_volume_sum2, 'tif', ...
                'WriteMode', 'append');
        end
        page = page + 1;
    end
end
disp('Volume 3D sommé2 sauvegardé avec imwrite !');

% Sauvegarde du MAT pour le volume sommé2
save('exvivo_warped_subvolume_sum.mat', 'Volume_Warped_Sum2', 'nb_z_sum2', 'a_garde_dernier', '-v7.3');

fprintf('\n=== Tout est terminé et sauvegardé ! ===\n');
fprintf('- Volume original : %s\n', out_tif_volume);
fprintf('- Volume sommé2   : %s\n', out_tif_volume_sum2);
fprintf('- Overlay         : %s\n', out_tif_overlay);