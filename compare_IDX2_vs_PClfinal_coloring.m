%% Comparaison directe : coloration par IDX2 (K-means, meme espace que l'UMAP)
%% vs coloration par PCl_final (test de permutation sur cellules d'assemblee)
%
% Objectif : verifier si les "points isoles" observes avec PCl_final
% disparaissent quand on colore avec IDX2 -- ce qui confirmerait que le
% decalage vient du critere de coloration, pas d'un probleme de clustering.
%
% Necessite : Y_umap, IDX2 (sortie de kmeansoptnew), NClOK, PCl_final, NCl

NRace = size(Y_umap, 1);

%% 1. Label IDX2 restreint aux clusters valides (meme convention que votre pipeline)
label_idx2 = IDX2;
label_idx2(IDX2 > NClOK) = 0;   % 0 = cluster non retenu / bruit, comme Cl0

%% 2. Label PCl_final (ce qui a ete affiche jusqu'ici)
label_pclfinal = zeros(1, NRace);
for j = 1:NRace
    a = find(PCl_final(:, j) == 1, 1, 'first');
    if ~isempty(a), label_pclfinal(j) = a; end
end

%% 3. Accord entre les deux criteres (devrait etre imparfait par construction)
ari_idx2_vs_pcl = rand_index(label_idx2', label_pclfinal, 'adjusted');
fprintf('ARI entre IDX2 (K-means) et PCl_final (permutation) : %.3f\n', ari_idx2_vs_pcl);

%% 4. Figure comparative
palette = lines(max(max(label_idx2), NCl) + 1);

fig = figure('visible', 'off', 'Color', 'w', 'Position', [100 100 1200 500]);

subplot(1,2,1);
gscatter(Y_umap(:,1), Y_umap(:,2), label_idx2, [], [], 12);
title('Coloration par IDX2 (K-means, meme espace que UMAP)');
xlabel('UMAP 1'); ylabel('UMAP 2'); legend off; box on;

subplot(1,2,2);
gscatter(Y_umap(:,1), Y_umap(:,2), label_pclfinal, [], [], 12);
title('Coloration par PCl\_final (test permutation, cellules assemblee)');
xlabel('UMAP 1'); ylabel('UMAP 2'); legend off; box on;

sgtitle(sprintf('Meme embedding, deux criteres de coloration differents (ARI=%.2f)', ari_idx2_vs_pcl));

namegraph = strcat(namefull, 'compare_IDX2_vs_PClfinal.png');
if isfolder(namefull)
    exportgraphics(fig, namegraph, 'Resolution', 200, 'BackgroundColor', 'white');
end
close(fig);

%% 5. Verification specifique : les "points isoles" en PCl_final sont-ils
%% bien groupes en IDX2 ?
% (reprend le calcul d'outliers du script diagnose_outlier_events.m)
dist_to_centroid = NaN(1, NRace);
for a = 1:NCl
    idx_a = find(label_pclfinal == a);
    if numel(idx_a) < 2, continue; end
    centroid_a = median(Y_umap(idx_a, :), 1);
    dist_to_centroid(idx_a) = vecnorm(Y_umap(idx_a, :) - centroid_a, 2, 2);
end

is_outlier_pcl = false(1, NRace);
for a = 1:NCl
    idx_a = find(label_pclfinal == a);
    if numel(idx_a) < 5, continue; end
    thresh = prctile(dist_to_centroid(idx_a), 95);
    is_outlier_pcl(idx_a(dist_to_centroid(idx_a) > thresh)) = true;
end

% Pour ces outliers PCl_final, quel est leur label IDX2 ?
fprintf('\nPour les %d events "isoles" selon PCl_final :\n', sum(is_outlier_pcl));
tabulate(label_idx2(is_outlier_pcl));
% Si ces outliers appartiennent presque tous a un MEME cluster IDX2 (souvent
% differt de leur label PCl_final), ca confirme qu'ils sont mal classes par
% le test de permutation mais bien regroupes par le K-means sur la
% similarite de population complete.