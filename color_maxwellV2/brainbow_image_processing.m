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
    % CORRECTION du 11 nov. : avant, on gardait les n_pixels premières valeurs
    % du vecteur puis on remettait en 512 x 512. Avec une LIGNE en trop
    % (513 x 512), chaque colonne se décalait d'un pixel de plus que la
    % précédente (image cisaillée). Maintenant l'image est rognée (ou
    % complétée par des zéros) ligne par ligne et colonne par colonne.
    img_d = double(img_raw(:, :, 1));
    if ~isequal(size(img_d), target_size)
        warning('Image %s : %d x %d au lieu de %d x %d, rognée / complétée.', name, size(img_d, 1), ...
            size(img_d, 2), target_size(1), target_size(2));
        tmp = zeros(target_size);
        nl = min(size(img_d, 1), target_size(1));
        nc = min(size(img_d, 2), target_size(2));
        tmp(1:nl, 1:nc) = img_d(1:nl, 1:nc);
        img_d = tmp;
    end
    vec_double = img_d(:);

    % c. Calcul des bornes (Robust Min-Max)
    val_min = min(vec_double);
    val_max_99 = prctile(vec_double, 99.9);
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
