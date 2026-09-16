%% Selection automatique d'epsilon par detection de coude (knee/elbow)
% Remplace le prctile(kdist,90) trop permissif par une vraie detection
% du coude du k-distance plot (methode "kneedle" simplifiee : point le
% plus eloigne de la droite reliant le debut et la fin de la courbe).

function epsilon = select_epsilon_knee(Y, minpts, do_plot)
if nargin < 3, do_plot = false; end

D = pdist2(Y, Y);
D_sorted = sort(D, 2);
kdist = sort(D_sorted(:, minpts), 'descend');   % courbe decroissante
n = numel(kdist);
x = (1:n)';

% Normaliser x et y entre 0 et 1 pour rendre la geometrie comparable
xn = (x - x(1)) / (x(end) - x(1));
yn = (kdist - min(kdist)) / (max(kdist) - min(kdist));

% Droite reliant le premier et le dernier point
p1 = [xn(1), yn(1)];
p2 = [xn(end), yn(end)];
line_vec = p2 - p1;
line_vec = line_vec / norm(line_vec);

% Distance perpendiculaire de chaque point a cette droite
pts = [xn, yn] - p1;
proj_len = pts * line_vec';
proj = proj_len * line_vec;
perp_dist = vecnorm(pts - proj, 2, 2);

[~, knee_idx] = max(perp_dist);
epsilon = kdist(knee_idx);

if do_plot
    figure('visible','on');
    plot(x, kdist, 'b-', 'LineWidth', 1.2); hold on;
    plot(knee_idx, kdist(knee_idx), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    xlabel(sprintf('Points tries (%d-eme voisin decroissant)', minpts));
    ylabel('Distance');
    title(sprintf('k-distance plot -- coude detecte : epsilon = %.4f', epsilon));
    legend('k-distance', 'Coude detecte (epsilon retenu)');
    grid on; hold off;
end
end

%% Utilisation dans votre script (remplace le bloc "Estimation semi-automatique")
minpts = 8;
epsilon_umap = select_epsilon_knee(Y_umap, minpts, true);   % true = affiche le plot
epsilon_tsne = select_epsilon_knee(Y_tsne, minpts, true);

fprintf('Epsilon retenu (UMAP, via coude) : %.4f\n', epsilon_umap);
fprintf('Epsilon retenu (t-SNE, via coude) : %.4f\n', epsilon_tsne);

idx_dbscan_umap = dbscan(Y_umap, epsilon_umap, minpts);
idx_dbscan_tsne = dbscan(Y_tsne, epsilon_tsne, minpts);

fprintf('UMAP  : %d clusters trouves (hors bruit), %.1f%% de bruit\n', ...
    numel(unique(idx_dbscan_umap(idx_dbscan_umap>0))), 100*mean(idx_dbscan_umap==-1));
fprintf('t-SNE : %d clusters trouves (hors bruit), %.1f%% de bruit\n', ...
    numel(unique(idx_dbscan_tsne(idx_dbscan_tsne>0))), 100*mean(idx_dbscan_tsne==-1));

%% Verification complementaire : balayage d'epsilon autour du coude
% Utile si le coude automatique donne encore un resultat degenere
% (1 seul cluster geant, ou tout en bruit)
eps_range = linspace(epsilon_umap*0.3, epsilon_umap*2, 15);
n_clusters_vs_eps = zeros(size(eps_range));
pct_noise_vs_eps  = zeros(size(eps_range));

for e = 1:numel(eps_range)
    idx_tmp = dbscan(Y_umap, eps_range(e), minpts);
    n_clusters_vs_eps(e) = numel(unique(idx_tmp(idx_tmp > 0)));
    pct_noise_vs_eps(e)  = 100 * mean(idx_tmp == -1);
end

figure('visible','on');
yyaxis left
plot(eps_range, n_clusters_vs_eps, '-o', 'LineWidth', 1.5);
ylabel('Nombre de clusters');
yyaxis right
plot(eps_range, pct_noise_vs_eps, '-s', 'LineWidth', 1.5);
ylabel('% de points bruit');
xlabel('epsilon'); grid on;
title('Sensibilite du clustering DBSCAN a epsilon (UMAP)');
xline(epsilon_umap, '--k', 'epsilon retenu (coude)');