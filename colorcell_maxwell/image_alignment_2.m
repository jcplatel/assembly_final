function [aligned_red, aligned_green, aligned_blue, aligned_red890] = image_alignment_2(images_norm)


% ÉTAPE 1 : Alignement de base rigide (imregcorr trouve la translation et rotation d'un coup, sans tâtonner)
tform_init = imregcorr(images_norm.red1040, images_norm.calcium, 'transformtype', 'rigid'); 

% ÉTAPE 2 : On configure l'optimiseur multimodal
[optimizer, metric] = imregconfig('multimodal');

% ÉTAPE 3 : On affine l'alignement pour trouver le cisaillement et le scaling anisotrope
% En utilisant 'InitialTransformation', imregtform commence son tâtonnement très près du but !
tform_finale = imregtform(images_norm.red1040, images_norm.calcium, 'affine', optimizer, metric, ...
    'InitialTransformation', tform_init);

% ÉTAPE 4 : On déforme l'image finale
aligned_red = imwarp(images_norm.red1040, tform_finale, 'OutputView', imref2d(size(images_norm.calcium)));
aligned_green = imwarp(images_norm.green, tform_finale, 'OutputView', imref2d(size(images_norm.calcium)));

tform_finale = imregtform(images_norm.red890, aligned_red, 'affine', optimizer, metric, ...
    'InitialTransformation', tform_init);

aligned_red890 = imwarp(images_norm.red890, tform_finale, 'OutputView', imref2d(size(images_norm.calcium)));
aligned_blue = imwarp(images_norm.blue, tform_finale, 'OutputView', imref2d(size(images_norm.calcium)));

end