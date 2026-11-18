function v = est_v2(chemin)
%EST_V2  true si le chemin est dans une analyse « version 2 » des assemblées
%   (11 nov.) : <dossier de l'étape 1>\analyses\<nom>\..., avec <nom> = 'v2',
%   ou qui commence par 'v2_', ou qui finit par '_v2'. La version est une
%   propriété de l'analyse : version 1 (RACEORTHO_FDR2, une assemblée par
%   cellule) et version 2 (RACEORTHO_V2, cellules partagées) ne se mélangent
%   jamais dans un même dossier. Voir fichier_assemblees.
v = false;
if isempty(chemin), return; end
tok = regexp(char(chemin), '[\\/]analyses[\\/]([^\\/]+)', 'tokens', 'once');
if isempty(tok), return; end
v = nom_v2(tok{1});
end


function v = nom_v2(nom)
v = ~isempty(regexpi(nom, '(^v2$|^v2_|_v2$)', 'once'));
end
