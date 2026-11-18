function f = trouver_composite(dossier)
%TROUVER_COMPOSITE  Image composite d'une session (11 nov.) : d'abord
%   compositecalcium.tif (rouge, vert, bleu, calcium, écrit par
%   color_predictor_model : l'image qui a servi aux couleurs), sinon un autre
%   fichier dont le nom contient « composite », dans le dossier des données
%   puis dans son dossier parent. '' si rien.
f = '';
d = regexprep(char(dossier), '[\\/]+$', '');
for niveau = 1:2
    for motif = {'compositecalcium.tif', 'compositecalcium.tiff', '*omposite*.tif', '*omposite*.tiff'}
        L = dir(fullfile(d, motif{1}));
        if ~isempty(L)
            f = fullfile(L(1).folder, L(1).name);
            return
        end
    end
    d = fileparts(d);
end
end
