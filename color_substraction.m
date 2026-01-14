function [Tr1b,colorcell] = color_substraction (Tr1b,path) 

fileExists_colorcellnew = isfile(strcat(path ,'colorcellnew.mat'));
fileExists_colorcell = isfile(strcat(path ,'colorcell.mat'));
if fileExists_colorcellnew
    load (strcat(path ,'colorcellnew.mat'))
elseif fileExists_colorcell
    load (strcat(path ,'colorcell.mat'))
end
load("E:\Data\Aurelie\data\nocues\444119\220919_plane0\colorcell_registration.mat");%colorcell from Solene;%for444119 220919plane0
Tr1b = double(F);

if exist('colorcell','var')
    mask_to_keep = colorcell < 7;
    Tr1b = Tr1b(mask_to_keep, :);
    colorcell = colorcell(mask_to_keep);
end