%% =================================================================
% UMAP_assembly.m -- VERSION PROPRE ET AUTONOME (corrigee)
% =================================================================
% Appele dans la boucle "for nanalysis = 4:16" de votre script principal,
% juste apres SCE_clustering (qui fournit Race, PCl_final, NCl, colors,
% namefull, NCell, P_VAL_THRESHOLD, MIN_CELL_SIZE, NShuf, NClini).
%
% Ce script :
%   1. Calcule l'embedding UMAP + t-SNE des events RACE
%   2. Detecte epsilon automatiquement (coude du k-distance plot)
%   3. Clusterise avec DBSCAN (independant de tout NClini/K-means)
%   4. Compare ce clustering a votre resultat officiel (PCl_final) via ARI
%   5. Sauvegarde les metriques dans des tableaux indexes par nanalysis,
%      pour un graphique recapitulatif a la fin de la boucle complete
%
% Necessite : run_umap.m (plugin), Statistics and Machine Learning Toolbox
% (dbscan, tsne, silhouette), et rand_index.m (FileExchange, deja installe)

NRace = size(Race, 2);

%% 1. Embedding UMAP + t-SNE
X = double(Race)';

nComp = min(50, size(X, 2));
[~, score] = pca(X, 'NumComponents', nComp);

rng(0);
Y_tsne = tsne(score, 'NumDimensions', 2, 'Perplexity', 30);
rng(0, 'twister')
Y_umap = run_umap(score, 'n_neighbors', 5, 'min_dist', 0.01, ...
    'n_components', 2, 'randomize', false, 'verbose', 'none');

%% 2. Detection d'epsilon (coude du k-distance plot) -- fonction locale en fin de fichier
minpts = 8;
% epsilon_umap = local_select_epsilon_knee(Y_umap, minpts);
% epsilon_tsne = local_select_epsilon_knee(Y_tsne, minpts);

epsilon_umap = local_select_epsilon_knee(Y_umap, minpts, strcat(namefull, 'kdist_umap.png'));
epsilon_tsne = local_select_epsilon_knee(Y_tsne, minpts, strcat(namefull, 'kdist_tsne.png'));

%% 3. Clustering DBSCAN
idx_dbscan_umap = dbscan(Y_umap, epsilon_umap, minpts);
idx_dbscan_tsne = dbscan(Y_tsne, epsilon_tsne, minpts);

% Garde-fou : si le coude donne un resultat degenere (tout fusionne en un
% seul cluster), on relache epsilon en prenant un percentile plus bas
if numel(unique(idx_dbscan_umap(idx_dbscan_umap > 0))) <= 1
    D_fb = pdist2(Y_umap, Y_umap);
    D_fb_sorted = sort(D_fb, 2);
    kdist_fb = sort(D_fb_sorted(:, minpts), 'descend');
    epsilon_umap = prctile(kdist_fb, 50);
    idx_dbscan_umap = dbscan(Y_umap, epsilon_umap, minpts);
end

fprintf('[NClini=%d] UMAP : eps=%.3f, %d clusters, %.1f%% bruit\n', ...
    NClini, epsilon_umap, numel(unique(idx_dbscan_umap(idx_dbscan_umap>0))), ...
    100*mean(idx_dbscan_umap==-1));

%% 4. Comparaison avec votre resultat officiel (PCl_final)
label_perm = zeros(1, NRace);
for j = 1:NRace
    a = find(PCl_final(:, j) == 1, 1, 'first');
    if ~isempty(a), label_perm(j) = a; end
end

ari_umap = rand_index(label_perm, idx_dbscan_umap', 'adjusted');
ari_tsne = rand_index(label_perm, idx_dbscan_tsne', 'adjusted');
fprintf('[NClini=%d] ARI vs DBSCAN : UMAP=%.3f | t-SNE=%.3f\n', NClini, ari_umap, ari_tsne);

%% 5. Stockage pour le graphique recapitulatif de fin de boucle
% (variables persistantes dans le workspace du script appelant, comme inertia_all)
ari_umap_all(NClini) = ari_umap;
ari_tsne_all(NClini) = ari_tsne;
epsilon_umap_all(NClini) = epsilon_umap;
NCl_dbscan_all(NClini) = numel(unique(idx_dbscan_umap(idx_dbscan_umap > 0)));

%% 6. Figure comparative (sauvegardee a chaque NClini)
colors_asm = colors(1:NCl, :);
color_noise = [0.8 0.8 0.8];
event_color = zeros(NRace, 3);
is_complex = false(1, NRace);
for j = 1:NRace
    active = find(PCl_final(:, j) == 1);
    if isempty(active)
        event_color(j, :) = color_noise;
    elseif numel(active) == 1
        event_color(j, :) = colors_asm(active, :);
    else
        event_color(j, :) = mean(colors_asm(active, :), 1);
        is_complex(j) = true;
    end
end

fig = figure('visible', 'off', 'Color', 'w', 'Position', [100 100 1500 500]);

subplot(1,3,1); hold on;
scatter(Y_umap(~is_complex,1), Y_umap(~is_complex,2), 18, event_color(~is_complex,:), 'filled', 'MarkerFaceAlpha', 0.7);
scatter(Y_umap(is_complex,1), Y_umap(is_complex,2), 40, event_color(is_complex,:), 'd', 'filled', 'MarkerEdgeColor','k');
title('UMAP - assemblee (PCl\_final)'); xlabel('UMAP 1'); ylabel('UMAP 2'); box on; hold off;

subplot(1,3,2);
gscatter(Y_umap(:,1), Y_umap(:,2), idx_dbscan_umap, [], [], 15);
title(sprintf('UMAP - DBSCAN (eps=%.2f, ARI=%.2f)', epsilon_umap, ari_umap));
xlabel('UMAP 1'); ylabel('UMAP 2'); legend off; box on;

subplot(1,3,3);
gscatter(Y_tsne(:,1), Y_tsne(:,2), idx_dbscan_tsne, [], [], 15);
title(sprintf('t-SNE - DBSCAN (eps=%.2f, ARI=%.2f)', epsilon_tsne, ari_tsne));
xlabel('t-SNE 1'); ylabel('t-SNE 2'); legend off; box on;

sgtitle(sprintf('NClini=%d : permutation vs clustering direct sur embedding', NClini));

if isfolder(namefull)
    exportgraphics(fig, strcat(namefull, 'umap_dbscan_comparison.png'), ...
        'Resolution', 200, 'BackgroundColor', 'white');
end
close(fig);

%% =================================================================
% FONCTION LOCALE (doit rester en fin de fichier)
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

    % --- Trace et sauvegarde (pas d'affichage interactif) ---
    if ~isempty(save_path)
        fig = figure('visible', 'off');
        plot(x, kdist, 'b-', 'LineWidth', 1.2); hold on;
        plot(knee_idx, kdist(knee_idx), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
        xlabel(sprintf('Points tries (%d-eme voisin decroissant)', minpts));
        ylabel('Distance');
        title(sprintf('k-distance plot -- epsilon retenu = %.4f', epsilon));
        legend('k-distance', 'Coude detecte'); grid on; hold off;
        exportgraphics(fig, save_path, 'Resolution', 150, 'BackgroundColor', 'white');
        close(fig);
    end
end