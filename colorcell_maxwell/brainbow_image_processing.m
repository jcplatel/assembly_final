function [images_norm] = brainbow_image_processing(path,red, red890,green, blue,calcium)

%%filter images (background substraction rolling ball): 
radius = 50;
se = strel('disk', radius);

% Perform the rolling ball background subtraction
green_filt = imtophat(green, se);
red1040_filt  = imtophat(red, se);
blue_filt  = imtophat(blue, se);
red890_filt  = imtophat(red890, se);

%% Min Max
images_input.green = green_filt;
images_input.blue = blue_filt;
images_input.red1040 = red1040_filt;
images_input.calcium = calcium;
images_input.red890 = red890_filt;

channel_names = fieldnames(images_input);

% Structure pour stocker les résultats
images_norm = struct();

% Paramètres
target_size = [512, 512];
n_pixels = prod(target_size);

% figure('Name', 'Comparaison des canaux normalisés', Visible="on");

for i = 1:length(channel_names)

    % a. Récupération du nom et de l'image courante
    name = channel_names{i};
    img_raw = images_input.(name);

    % b. Conversion et Nettoyage (Pixel en trop)
    vec_double = double(img_raw(:));

    if numel(vec_double) > n_pixels
        vec_double = vec_double(1:n_pixels);
    elseif numel(vec_double) < n_pixels
        warning('Image %s trop petite ! Remplie avec des zéros.', name);
        vec_double(n_pixels) = 0; % Padding basique si besoin
    end

    % c. Calcul des bornes (Robust Min-Max)
    val_min = min(vec_double);
    val_max_99 = prctile(vec_double, 99.99);
    % val_min = prctile(vec_double, 1);
    % val_max_99 = prctile(vec_double, 99);

    % d. Normalisation
    if val_max_99 > val_min
        vec_norm = (vec_double - val_min) / (val_max_99 - val_min);
    else
        vec_norm = zeros(size(vec_double));
    end

    % e. Clipping (Saturation du top 1%)
    vec_norm(vec_norm > 1) = 1;
    vec_norm(vec_norm < 0) = 0;

    % f. Reshape et Stockage
    img_final = reshape(vec_norm, target_size);
    images_norm.(name) = img_final;

    % g. Affichage immédiat pour vérifier
    % subplot(2, 3, i);
    % imshow(img_final);
    % title([name ' (Norm 99%)']);
end

% for i = 1:length(channel_names)
%     name = channel_names{i};
%     images_norm.(name) = images_input.(name);
% end



end
