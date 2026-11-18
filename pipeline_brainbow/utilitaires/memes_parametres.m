function ok = memes_parametres(fichier, nom_variable, valeur)
% Vrai si le fichier existe, contient nom_variable, et si sa valeur est
% identique à valeur. Sert à ne pas refaire un calcul déjà à jour.
ok = false;
try
    contenu = whos('-file', fichier);
    if ~ismember(nom_variable, {contenu.name}), return; end
    s = load(fichier, nom_variable);
    ok = isequaln(s.(nom_variable), valeur);
catch
end
end
