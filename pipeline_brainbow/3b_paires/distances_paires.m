function [bornes, cumulatif, etiquettes] = distances_paires(opts)
%DISTANCES_PAIRES  Distances de l'étape 3b.
%   opts.distances = 'disques' (défaut depuis le 5 nov.) : toutes les paires
%     à moins de R µm, pour chaque R de opts.rayons_um (cumulatif).
%     bornes = [0, rayons], cumulatif = true.
%   opts.distances = 'tranches' : tranches [bornes(i), bornes(i+1)) de
%     opts.bornes_um, cumulatif = false.
%   etiquettes : texte de chaque distance ('<= 200' ou '20-40').
if isfield(opts, 'distances') && strcmp(opts.distances, 'disques')
    bornes = [0, opts.rayons_um(:)'];
    cumulatif = true;
    etiquettes = arrayfun(@(r) sprintf('<= %g', r), bornes(2:end), 'UniformOutput', false);
    etiquettes = strrep(etiquettes, '<= Inf', 'all');
else
    bornes = opts.bornes_um(:)';
    cumulatif = false;
    etiquettes = arrayfun(@(i) sprintf('%g-%g', bornes(i), bornes(i + 1)), 1:numel(bornes) - 1, ...
        'UniformOutput', false);
    etiquettes = strrep(etiquettes, '-Inf', '+');
end
end
