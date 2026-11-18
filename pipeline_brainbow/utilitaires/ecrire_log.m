function historique = ecrire_log(path_file, historique, nouveau_message)
% Journal de traitement : affiché dans la console et ajouté au fichier texte.
% (Même fonction qu'avant, sans la fenêtre uifigure, qui était invisible.)
historique{end+1} = nouveau_message;
fprintf('%s\n', nouveau_message);
fid = fopen(path_file, 'a', 'n', 'UTF-8');
if fid ~= -1
    fprintf(fid, '%s\r\n', nouveau_message);
    fclose(fid);
end
end
