%% =================================================================
% UMAP_assembly_Raster.m -- Version sur la matrice Raster (Cellules x Temps)
% =================================================================
% Difference conceptuelle avec UMAP_assembly.m (sur Race) :
%   - Avant : chaque point UMAP = un evenement RACE, decrit par les
%     cellules actives pendant cet evenement (X = Race' -> NRace x NCell)
%   - Maintenant : chaque point UMAP = une CELLULE, decrite par son
%     profil d'activite binaire au cours du temps (X = Raster -> NCell x Nz)
%
% Ceci permet de verifier si les cellules appartenant a la meme assemblee
% (definie par votre methode de permutation, CellCl_final/assemblyortho)
% se regroupent naturellement dans un espace de reduction de dimension
% construit uniquement a partir de leur activite temporelle brute --
% une validation independante de vos assemblees, au niveau cellulaire.
%
% Necessite : Raster (NCell x Nz), assemblyortho (cell array), CellCl_final
% (NCl x NCell, sortie de RaceOrthoJC), colors, namefull, NClini.

[NCell, Nz] = size(Raster);

%% 1. Embedding UMAP + t-SNE (points = cellules, variables = temps)
X = double(Raster); % NCell x Nz : une ligne = une cellule

% Nombre de composantes PCA limite par le nombre d'observations (NCell),
% pas par le nombre de frames (Nz), contrairement a la version Race.
nComp = min(50, NCell - 1);
[~, score] = pca(X, 'NumComponents', nComp);

rng(0);
Y_tsne = tsne(score, 'NumDimensions', 2, 'Perplexity', min(30, floor((NCell-1)/3)));
rng(0, 'twister')
Y_umap = run_umap(score, 'n_neighbors', 5, 'min_dist', 0.01, ...
    'n_components', 2, 'randomize', false, 'verbose', 'none');

%% 2. Detection d'epsilon (coude du k-distance plot)
minpts = min(8, max(3, round(NCell/20))); % adapte si peu de cellules
epsilon_umap = local_select_epsilon_knee(Y_umap, minpts, strcat(namefull, 'kdist_umap_raster.png'));
epsilon_tsne = local_select_epsilon_knee(Y_tsne, minpts, strcat(namefull, 'kdist_tsne_raster.png'));

%% 3. Clustering DBSCAN
idx_dbscan_umap = dbscan(Y_umap, epsilon_umap, minpts);
idx_dbscan_tsne = dbscan(Y_tsne, epsilon_tsne, minpts);

% Garde-fou : si le coude donne un resultat degenere (tout fusionne)
if numel(unique(idx_dbscan_umap(idx_dbscan_umap > 0))) <= 1
    D_fb = pdist2(Y_umap, Y_umap);
    D_fb_sorted = sort(D_fb, 2);
    kdist_fb = sort(D_fb_sorted(:, minpts), 'descend');
    epsilon_umap = prctile(kdist_fb, 50);
    idx_dbscan_umap = dbscan(Y_umap, epsilon_umap, minpts);
end

fprintf('[NClini=%d] UMAP (Raster) : eps=%.3f, %d clusters, %.1f%% bruit\n', ...
    NClini, epsilon_umap, numel(unique(idx_dbscan_umap(idx_dbscan_umap>0))), ...
    100*mean(idx_dbscan_umap==-1));

%% 4. Comparaison avec votre resultat officiel (CellCl_final, niveau cellule)
% CellCl_final : [NCl x NCell], 1 si la cellule appartient a l'assemblee
label_perm = zeros(1, NCell);
for c = 1:NCell
    a = find(CellCl_final(:, c) == 1, 1, 'first');
    if ~isempty(a), label_perm(c) = a; end
end

ari_umap = rand_index(label_perm, idx_dbscan_umap', 'adjusted');
ari_tsne = rand_index(label_perm, idx_dbscan_tsne', 'adjusted');
fprintf('[NClini=%d] ARI vs DBSCAN (Raster) : UMAP=%.3f | t-SNE=%.3f\n', NClini, ari_umap, ari_tsne);

%% 5. Stockage pour le graphique recapitulatif de fin de boucle
ari_umap_raster_all(NClini) = ari_umap;
ari_tsne_raster_all(NClini) = ari_tsne;
epsilon_umap_raster_all(NClini) = epsilon_umap;
NCl_dbscan_raster_all(NClini) = numel(unique(idx_dbscan_umap(idx_dbscan_umap > 0)));

%% 6. Figure comparative : projection des assemblees + DBSCAN
NCl = size(CellCl_final, 1);
colors_asm = colors(1:NCl, :);
color_noise = [0.8 0.8 0.8];

% Une cellule appartient a au plus une assemblee (deja resolu par le
% traitement "promiscuous" en amont) -> pas de couleur mixte necessaire
cell_color = repmat(color_noise, NCell, 1);
for a = 1:NCl
    cells_a = find(CellCl_final(a, :) == 1);
    cell_color(cells_a, :) = repmat(colors_asm(a, :), numel(cells_a), 1);
end
in_assembly = label_perm > 0;

fig = figure('visible', 'on', 'Color', 'w', 'Position', [100 100 1500 500]);

subplot(1,3,1); hold on;
scatter(Y_umap(~in_assembly,1), Y_umap(~in_assembly,2), 18, cell_color(~in_assembly,:), 'filled', 'MarkerFaceAlpha', 0.5);
scatter(Y_umap(in_assembly,1), Y_umap(in_assembly,2), 40, cell_color(in_assembly,:), 'filled', 'MarkerEdgeColor','k');
title('UMAP (Raster) - assemblees (CellCl\_final)'); xlabel('UMAP 1'); ylabel('UMAP 2'); box on; hold off;

subplot(1,3,2);
gscatter(Y_umap(:,1), Y_umap(:,2), idx_dbscan_umap, [], [], 15);
title(sprintf('UMAP (Raster) - DBSCAN (eps=%.2f, ARI=%.2f)', epsilon_umap, ari_umap));
xlabel('UMAP 1'); ylabel('UMAP 2'); legend off; box on;

subplot(1,3,3);
gscatter(Y_tsne(:,1), Y_tsne(:,2), idx_dbscan_tsne, [], [], 15);
title(sprintf('t-SNE (Raster) - DBSCAN (eps=%.2f, ARI=%.2f)', epsilon_tsne, ari_tsne));
xlabel('t-SNE 1'); ylabel('t-SNE 2'); legend off; box on;

sgtitle(sprintf('NClini=%d : assemblees (permutation) vs clustering direct sur activite temporelle', NClini));

if isfolder(namefull)
    exportgraphics(fig, strcat(namefull, 'umap_dbscan_comparison_raster.png'), ...
        'Resolution', 200, 'BackgroundColor', 'white');
end
% close(fig);

%% =================================================================
% FONCTION LOCALE (identique a la version Race, doit rester en fin de fichier)
% =================================================================

function epsilon = local_select_epsilon_knee(Y, minpts, save_path)
    if nargin < 3, save_path = ''; end

    D = pdist2(Y, Y);
    D_sorted = sort(D, 2);
    kdist = sort(D_sorted(:, minpts), 'descend');
    n = numel(kdist);
    x = (1:n)';

    xn = (x - x(1)) / (x(end) - x(1));
    yn = (kdist - min(kdist)) / (max(kdist) - min(kdist));

    p1 = [xn(1), yn(1)];
    p2 = [xn(end), yn(end)];
    line_vec = (p2 - p1) / norm(p2 - p1);
    pts = [xn, yn] - p1;
    proj = (pts * line_vec') * line_vec;
    perp_dist = vecnorm(pts - proj, 2, 2);
    [~, knee_idx] = max(perp_dist);
    epsilon = kdist(knee_idx);

    % if ~isempty(save_path)
        fig = figure('visible', 'on');
        plot(x, kdist, 'b-', 'LineWidth', 1.2); hold on;
        plot(knee_idx, kdist(knee_idx), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
        xlabel(sprintf('Points tries (%d-eme voisin decroissant)', minpts));
        ylabel('Distance');
        title(sprintf('k-distance plot -- epsilon retenu = %.4f', epsilon));
        legend('k-distance', 'Coude detecte'); grid on; hold off;
        % exportgraphics(fig, save_path, 'Resolution', 150, 'BackgroundColor', 'white');
        % close(fig);
    % end
end