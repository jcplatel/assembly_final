%% demoRaceHierarchicalClustering.m
% Script de démonstration : génère un jeu de données RACE synthétique,
% appelle raceHierarchicalClustering, et affiche les graphiques :
%   1) dendrogramme des SCE avec la coupe statistique
%   2) matrice RACE réordonnée selon les clusters
%   3) histogramme des p-values de biais de couleur le long de la hiérarchie
%
% Remplace le bloc "DONNEES SYNTHETIQUES" par ta vraie matrice RACE et
% ton vecteur colorLabel pour l'utiliser sur tes données.

RACE=Race';
fprintf('Seuil statistique de coupe (distance cosinus) : %.4f\n', results.threshold);
fprintf('Nombre d''assemblages significatifs detectes : %d\n', numel(results.assemblies));
for a = 1:numel(results.assemblies)
    fprintf('  Assemblage %d : %d SCE, %d cellules coeur -> [%s]\n', a, ...
        numel(results.assemblies(a).sceIdx), numel(results.assemblies(a).coreCells), ...
        num2str(results.assemblies(a).coreCells(:)'));
end

%% ---------------- FIGURE 1 : dendrogramme + coupe ----------------
figure('Name','Dendrogramme des SCE','Color','w');
[H, T_, outperm] = dendrogram(results.Z, 0, 'ColorThreshold', results.threshold);
hold on;
yline(results.threshold, 'r--', 'Seuil statistique', 'LineWidth', 1.5);
xlabel('SCE (feuilles reordonnees)');
ylabel('Distance cosinus');
title('Clustering hierarchique des evenements synchrones (SCE)');

%% ---------------- FIGURE 2 : matrice RACE reordonnee ----------------
figure('Name','Matrice RACE reordonnee','Color','w');
[~, sceOrder] = sort(results.sceClusterID);
imagesc(RACE(sceOrder, :)');
colormap(flipud(gray));
ylabel('Cellules'); xlabel('SCE (triees par cluster)');
title('RACE trie par assemblage');
colorbar;

%% ---------------- FIGURE 3 : biais de couleur le long de la hierarchie ----------------
if ~isempty(results.colorBias)
    figure('Name','Biais de couleur dans la hierarchie','Color','w');
    scatter(results.colorBias.MergeDistance, -log10(results.colorBias.pAdj), ...
        40, results.colorBias.nSCE, 'filled');
    colorbar; colormap(parula);
    xlabel('Distance de fusion (niveau dans la hierarchie)');
    ylabel('-log_{10}(p ajustee, FDR)');
    yline(-log10(0.05), 'r--', 'seuil p=0.05');
    title('Biais de couleur a chaque noeud du dendrogramme');
    disp(results.colorBias(1:min(10,end),:));
end
