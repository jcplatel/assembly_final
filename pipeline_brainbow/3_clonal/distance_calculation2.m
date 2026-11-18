function dist = distance_calculation2(gcamp_centroid, assemblyortho, colorcell, size_pixel, nshuffle, graine)
% = distance_calculation (seules les parties actives : la fin du script,
% après « return », n'était jamais exécutée).
%   size_pixel : µm par pixel (705/512 dans ton script)
% Changements :
%  - toutes les distances en µm (avant : distance_assembly en PIXELS et
%    meandistance_gcamp en µm, exportées côte à côte)
%  - distance_assembly recalculée à chaque appel (avant, une valeur restée
%    d'un K précédent avec plus d'assemblées entrait dans la moyenne)
%  - test « chaque couleur est-elle regroupée dans l'espace ? » : le mélange
%    utilise maintenant les positions des BONNES cellules. Avant,
%    colorcell_sh était plus court que NCell : gcamp_centroid(colorcell_sh==n,:)
%    prenait les premières cellules du champ, pas les cellules colorées.
%  - seules les couleurs 1 à 6 entrent dans le mélange (avant : colorcell < 7,
%    donc aussi un éventuel code 0)
%  - p-value = (1 + nb de mélanges aussi regroupés) / (nshuffle + 1) ;
%    NaN pour une couleur qui a moins de 2 cellules

%% Distances entre toutes les cellules (µm)
distances_gcamp = pdist(gcamp_centroid);
meandistance_gcamp = mean(distances_gcamp) * size_pixel;

%% Distances entre cellules d'une même assemblée (µm)
nbassemblies = numel(assemblyortho);
distance_assembly = nan(1, nbassemblies);
for num = 1:nbassemblies
    distance_assembly(num) = mean(pdist(gcamp_centroid(assemblyortho{num}, :))) * size_pixel;
end
mean_distance_assembly = mean(distance_assembly, 'omitnan');

%% Chaque couleur est-elle regroupée dans le champ ?
colorcell = colorcell(:);
real_color = find(colorcell >= 1 & colorcell <= 6);   % cellules colorées (ni blanches ni noires)
couleurs_reelles = colorcell(real_color);
mean_distances_color = nan(1, 6);
for n = 1:6
    mean_distances_color(n) = mean(pdist(gcamp_centroid(colorcell == n, :))) * size_pixel;
end
rng(graine);
mean_distances_color_sh = nan(nshuffle, 6);
for shuffle = 1:nshuffle
    colorcell_sh = couleurs_reelles(randperm(numel(real_color)));    % couleurs mélangées
    for n = 1:6
        % positions des cellules colorées qui reçoivent la couleur n
        mean_distances_color_sh(shuffle, n) = mean(pdist(gcamp_centroid(real_color(colorcell_sh == n), :))) * size_pixel;
    end
end
pvalue_distance = (1 + sum(mean_distances_color_sh <= mean_distances_color, 1)) / (nshuffle + 1);
pvalue_distance(isnan(mean_distances_color)) = NaN;

dist = struct('size_pixel', size_pixel, 'meandistance_gcamp', meandistance_gcamp, ...
    'distance_assembly', distance_assembly, 'mean_distance_assembly', mean_distance_assembly, ...
    'mean_distances_color', mean_distances_color, 'pvalue_distance', pvalue_distance);
end
