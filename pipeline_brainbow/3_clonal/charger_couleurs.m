function [colorcell, path_colorcell, source_couleur] = charger_couleurs(path, fichier_colorcell)
% Couleur de chaque cellule. SEULE fonction du pipeline qui lit les couleurs.
% Priorité : fichier indiqué dans la table (colonne fichier_colorcell),
% sinon colorcell_GT, puis colorcell_Maxwellnew, puis colorcell_Maxwell
% (même ordre qu'avant). source_couleur garde la trace du fichier utilisé
% (pour une analyse de sensibilité par source).
%   1=R 2=G 3=B 4=RG 5=RB 6=GB 7=All (blanc) 8=none

if nargin >= 2 && ~isempty(fichier_colorcell)
    path_colorcell = char(fichier_colorcell);
    source_couleur = 'table';
elseif isfile([path 'colorcell_GT.mat'])
    path_colorcell = [path 'colorcell_GT.mat'];
    source_couleur = 'GT';
elseif isfile([path 'colorcell_Maxwellnew.mat'])
    path_colorcell = [path 'colorcell_Maxwellnew.mat'];
    source_couleur = 'Maxwellnew';
elseif isfile([path 'colorcell_Maxwell.mat'])
    path_colorcell = [path 'colorcell_Maxwell.mat'];
    source_couleur = 'Maxwell';
else
    error('Aucun fichier colorcell_GT / _Maxwellnew / _Maxwell dans %s', path);
end
if ~isfile(path_colorcell)
    error('Fichier couleur introuvable : %s', path_colorcell);
end
S = load(path_colorcell, 'colorcell');
colorcell = double(S.colorcell(:));
end
