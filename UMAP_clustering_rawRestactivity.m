%% =================================================================
% CLUSTERING DES CELLULES A PARTIR DE LEUR ACTIVITE AU REPOS (WinRest)
% =================================================================
% Objectif : au lieu de clusteriser les evenements RACE (SCE), on
% clusterise directement les CELLULES selon leur profil d'activite
% pendant les periodes de repos, pour voir si une structure fonctionnelle
% (assemblees) existe deja hors des evenements synchrones detectes.
%
% HYPOTHESES SUR VOS VARIABLES (a adapter aux noms reels) :
%   Activity : matrice NCell x NFrames, activite continue de chaque
%              cellule sur toute la session (dF/F, activite deconvoluee,
%              ou trace binaire framewise -- PAS la matrice Race qui est
%              restreinte aux SCE)
%   WinRest  : soit un vecteur logique 1 x NFrames (true = frame de repos),
%              soit un vecteur d'indices de frames de repos
%
% Si vos donnees sont deja sous forme "Race"-like (binaire, cellule x
% evenement) mais pour les frames de repos, remplacez simplement
% ActivityRest par cette matrice directement (sautez l'etape 1).

%% 1. Extraction de l'activite pendant le repos
Activity= Tr1b;
if islogical(WinRest)
    idxRest = find(WinRest);
else
    idxRest = WinRest;   % deja une liste d'indices
end

ActivityRest = Activity(:, idxRest);      % NCell x NFrames_rest
% ActivityRest = Activity; 
NCell = size(ActivityRest, 1);
NFramesRest = size(ActivityRest, 2);
fprintf('Repos : %d frames sur %d cellules.\n', NFramesRest, NCell);

% Optionnel : z-score de chaque cellule pour ne pas laisser dominer les
% cellules a forte amplitude/frequence de decharge
ActivityRest_z = zscore(ActivityRest, 0, 2);

%% 2. Matrice de correlation entre cellules (repos uniquement)
% C(i,j) = correlation entre l'activite des cellules i et j pendant le repos
C = corr(ActivityRest_z');            % NCell x NCell
D = 1 - C;                             % distance simple (1 - correlation)
D = (D + D') / 2;                      % forcer la symetrie (erreurs numeriques)
D(1:NCell+1:end) = 0;                  % diagonale a 0

%% 3. Reduction de dimension : PCA -> t-SNE / UMAP
nComp = min(50, NFramesRest);
[~, score] = pca(ActivityRest_z, 'NumComponents', nComp);   % NCell x nComp

rng(0);
Y_tsne_cells = tsne(score, 'NumDimensions', 2, 'Perplexity', 30);

% UMAP (plugin run_umap deja installe)
Y_umap_cells = run_umap(score, 'n_neighbors', 20, 'min_dist', 0.1, ...
    'n_components', 2, 'randomize', false, 'verbose', 'none');

%% 4. Clustering direct sur l'embedding (DBSCAN)
minpts = 8;

D_emb = pdist2(Y_umap_cells, Y_umap_cells);
D_emb_sorted = sort(D_emb, 2);
kdist = sort(D_emb_sorted(:, minpts), 'descend');
epsilon_umap = prctile(kdist, 90);      % a affiner visuellement (k-distance plot)
idx_dbscan_cells_umap = dbscan(Y_umap_cells, epsilon_umap, minpts);

D_emb2 = pdist2(Y_tsne_cells, Y_tsne_cells);
D_emb2_sorted = sort(D_emb2, 2);
kdist2 = sort(D_emb2_sorted(:, minpts), 'descend');
epsilon_tsne = prctile(kdist2, 90);
idx_dbscan_cells_tsne = dbscan(Y_tsne_cells, epsilon_tsne, minpts);

% -1 = cellules "isolees" (pas d'assemblee au repos)
valid_clusters = unique(idx_dbscan_cells_umap(idx_dbscan_cells_umap > 0));
NCl_rest = numel(valid_clusters);
fprintf('Clustering repos (UMAP+DBSCAN) : %d assemblees fonctionnelles trouvees.\n', NCl_rest);

%% 5. Construction des assemblees "repos" (meme format que assemblyortho)
assembly_rest = cell(1, NCl_rest);
for k = 1:NCl_rest
    assembly_rest{k} = find(idx_dbscan_cells_umap == valid_clusters(k))';
end

%% 6. Comparaison avec les assemblees issues des SCE (assemblyortho)
% Construire un label par cellule pour chaque methode (0 = non assignee)
label_sce = zeros(1, NCell);
for n = 1:numel(assemblyortho)
    label_sce(cell2mat(assemblyortho(n))) = n;
end

label_rest = zeros(1, NCell);
for n = 1:NCl_rest
    label_rest(assembly_rest{n}) = n;
end

ari_sce_vs_rest = rand_index(label_sce, label_rest,'adjusted');
fprintf('ARI (assemblees SCE vs assemblees repos) : %.3f\n', ari_sce_vs_rest);
% ARI proche de 0 -> les deux structures sont independantes
% ARI proche de 1 -> les assemblees actives pendant les SCE sont deja
%                     visibles dans la correlation d'activite au repos

% Matrice de confusion brute (utile si adjustedRandIndex indisponible)
% crosstab(label_sce, label_rest)

%% 7. Figure : embedding des cellules colore par cluster "repos" vs assemblee SCE
colors_rest = lines(max(NCl_rest, 1));   % palette generique si NCl_rest > taille de "colors"

fig = figure('visible', 'off', 'Color', 'w', 'Position', [100 100 1300 500]);

subplot(1,2,1);
gscatter(Y_umap_cells(:,1), Y_umap_cells(:,2), idx_dbscan_cells_umap, [], [], 15);
title(sprintf('UMAP cellules (repos) - DBSCAN : %d clusters', NCl_rest));
xlabel('UMAP 1'); ylabel('UMAP 2'); legend off; box on;

subplot(1,2,2);
hold on;
cell_color_sce = repmat([0.8 0.8 0.8], NCell, 1);
for n = 1:numel(assemblyortho)
    cell_color_sce(cell2mat(assemblyortho(n)), :) = repmat(colors(min(n,size(colors,1)), :), ...
        numel(cell2mat(assemblyortho(n))), 1);
end
scatter(Y_umap_cells(:,1), Y_umap_cells(:,2), 20, cell_color_sce, 'filled');
title('Meme embedding, colore par assemblee SCE (reference)');
xlabel('UMAP 1'); ylabel('UMAP 2'); box on; hold off;

sgtitle('Cellules projetees selon leur activite au repos : cluster repos vs assemblee SCE');

namegraph = strcat(namefull, 'umap_cells_resting_vs_SCE.png');
if isfolder(namefull)
    exportgraphics(fig, namegraph, 'Resolution', 200, 'BackgroundColor', 'white');
    close(fig);
end

%% 8. (Optionnel) Silhouette de la partition "repos"
sil_rest = silhouette(Y_umap_cells, idx_dbscan_cells_umap);
fprintf('Silhouette moyenne (clusters repos, UMAP) : %.3f\n', mean(sil_rest));