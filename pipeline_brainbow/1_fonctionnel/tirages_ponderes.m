function comptes = tirages_ponderes(poids, n, NShuf, groupe, nGroupes)
% n éléments tirés SANS remise, avec une probabilité proportionnelle à
% poids (méthode d'Efraimidis et Spirakis 2006), répété NShuf fois.
%   poids  : 1 x nElements (0 = jamais tiré)
%   groupe : 1 x nElements, numéro de groupe de chaque élément (0 = aucun)
% Sortie : comptes(g, t) = nombre d'éléments du groupe g au tirage t.

cles = rand(NShuf, numel(poids)) .^ (1 ./ poids);   % les n plus grandes clés sont tirées
[~, idx] = maxk(cles, n, 2);                         % NShuf x n
g_tires = reshape(groupe(idx), size(idx));          % groupe de chaque élément tiré

comptes = zeros(nGroupes, NShuf);
for g = 1:nGroupes
    comptes(g, :) = sum(g_tires == g, 2)';
end
end
