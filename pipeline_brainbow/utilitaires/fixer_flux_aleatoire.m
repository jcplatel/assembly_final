function fixer_flux_aleatoire(graine, indice, nombre)
% Donne au générateur aléatoire un flux propre à (graine, indice).
% Les tirages sont reproductibles, y compris dans un parfor, et
% indépendants d'une itération à l'autre (flux mrg32k3a séparés).
flux = RandStream.create('mrg32k3a', 'NumStreams', max(nombre, indice), ...
    'StreamIndices', indice, 'Seed', graine);
RandStream.setGlobalStream(flux);
end
