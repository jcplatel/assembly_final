function [aligned_red, aligned_green, aligned_blue, calcium_norm,aligned_red890] = image_alignment(images_norm)

%registration images
% %red 1040 aligned to calcium 920
% images_norm.
reg_obj = imregcorr(images_norm.red1040, images_norm.calcium, 'Method', 'phasecorr');
% reg_obj = imregcorr(green, calcium);
T = reg_obj.T;
aligned_red = imwarp(images_norm.red1040, affine2d(T), 'OutputView', imref2d(size(images_norm.red1040)));
% imwrite(uint16(aligned_red), ([path,'aligned_red_new.tif']), 'tif');
% 
aligned_green = imwarp(images_norm.green, affine2d(T), 'OutputView', imref2d(size(images_norm.green)));
% imwrite(uint16(aligned_green), ([path,'aligned_green_new.tif']), 'tif');

% %red 890 aligned to red1040 and used to align Blue890 
reg_obj = imregcorr(images_norm.red890, aligned_red, 'Method', 'phasecorr');
T = reg_obj.T;
aligned_red890 = imwarp(images_norm.red890, affine2d(T), 'OutputView', imref2d(size(images_norm.red890)));
% imwrite(aligned_red890, ([path,'aligned_red890.tif']), 'tif');

aligned_blue = imwarp(images_norm.blue, affine2d(T), 'OutputView', imref2d(size(images_norm.blue)));
% imwrite(uint16(aligned_blue), ([path,'aligned_blue_new.tif']), 'tif');
calcium_norm = images_norm.calcium;

aligned_red890(aligned_red890 < 0) = 0;
aligned_red(aligned_red < 0) = 0;
aligned_green(aligned_green < 0) = 0;
aligned_blue(aligned_blue < 0) = 0;

end