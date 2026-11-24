function [red, red890,green, blue, calcium] = brainbow_image_read (path)

% read original brainbow§calcium tif from path

% blue=double(imread([path,'blue_new.tif']));
% red=double(imread([path,'red_new.tif']));
% red890=double(imread([path,'red890_new.tif']));
% green=double(imread([path,'green_new.tif']));
% % green=double(imread([path,'deconv_green.tif']));

blue=double(imread(strcat(path,'blue.tif')));
red=double(imread(strcat(path,'red.tif')));
red890=double(imread(strcat(path,'red890.tif')));
green=double(imread(strcat(path,'green.tif')));
calcium=double(imread(strcat(path,'calcium.tif')));
end