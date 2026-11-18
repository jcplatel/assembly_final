function executer_scripts_anciens(noms_scripts, V)
% Lance d'anciens scripts (que je n'ai pas relus : raster_rastermap, ...)
% dans un espace de travail qui
% contient les champs de V sous leurs noms habituels (Race, IDX2,
% assemblyortho, namefull, colorcell, ...).
%
% Une erreur n'arrête pas le pipeline : elle s'affiche en avertissement.
% Un script qui a besoin de colorcell pendant l'étape 1 échouera : c'est
% voulu, aucune couleur n'existe avant l'étape 3. Il faut alors le passer
% dans config.scripts_anciens_clonal.

liste_scripts__ = noms_scripts;
noms_variables__ = fieldnames(V);
for i_var__ = 1:numel(noms_variables__)
    eval([noms_variables__{i_var__} ' = V.(noms_variables__{i_var__});']);
end
clear V i_var__ noms_variables__ noms_scripts

for i_script__ = 1:numel(liste_scripts__)
    nom_script__ = liste_scripts__{i_script__};
    if ~exist(nom_script__, 'file')
        warning('Ancien script %s introuvable : ignoré.', nom_script__);
        continue
    end
    try
        eval(nom_script__);
    catch err__
        warning('Ancien script %s non exécuté : %s', nom_script__, err__.message);
    end
    close all
end
end
