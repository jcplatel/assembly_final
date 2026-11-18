function verifier_alignement(NCell, colorcell, outline_gcampx, nom)
% Couleurs, contours et activité doivent décrire EXACTEMENT les mêmes
% cellules. Sinon : erreur, jamais de décalage silencieux.
% (Le même nombre ne garantit pas le même ordre : il faut que colorcell ait
% été calculé sur les cellules iscell == 1 de Fall.mat, dans leur ordre.)
if numel(colorcell) ~= NCell
    error(['%s : colorcell a %d valeurs mais l''activité a %d cellules. Couleurs et ' ...
        'activité ne décrivent pas la même liste de cellules (curation suite2p modifiée, ' ...
        'ou colorcell calculé sur toutes les ROI).'], nom, numel(colorcell), NCell);
end
if numel(outline_gcampx) ~= NCell
    error('%s : Fall.mat contient %d cellules (iscell = 1) mais l''activité en a %d.', ...
        nom, numel(outline_gcampx), NCell);
end
end
