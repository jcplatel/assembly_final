function chemin = avec_separateur(chemin)
% Ajoute \ (ou /) à la fin d'un chemin de dossier s'il n'y est pas déjà.
% Tes scripts écrivent strcat(path, 'Fall.mat') : il faut le séparateur.
chemin = char(chemin);
if ~isempty(chemin) && ~any(chemin(end) == '/\')
    chemin = [chemin filesep];
end
end
