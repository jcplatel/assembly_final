%% Projection t-SNE des RACE events, coloree par assemblee(s) active(s)
% A placer juste apres RACEORTHOJC_NEW, une fois obtenus :
% Race, assemblyortho, PCl_final, CellCl_final, NCl, x1, RList, namefull
%
% IMPORTANT : PCl_final(i,j)==1 signifie que l'assemblee i est
% statistiquement active lors de l'event j (test de permutation).
% sum(PCl_final(:,j)) donne le NOMBRE d'assemblees co-actives pour cet
% event (0 = bruit/Cl0, 1 = Cl1, 2 = Cl2, >=3 = Cl3), independant de
% l'identite de l'assemblee -> on ne peut pas colorer par Cl0/Cl1/Cl2/Cl3
% directement, il faut regarder quelle(s) ligne(s) de PCl_final sont a 1.
colors = [
    1 0 0;     % Red 
    0 1 0;     % Green
    0 0 1;     % Blue
    1 0 1;     % Magenta
    0 1 1;     % Cyan
    0.5 0 0;   % Maroon
    0 0.5 0;   % Olive
    0 0 0.5;   % Navy
    0.5 0.5 0; % Olive Green
    0.5 0 0.5; % Purple
    0 0.5 0.5; % Teal
    0.75 0 0;  % Dark Red
    0 0.75 0;  % Lime Green
    0 0 0.75;  % Dark Blue
    0.75 0.75 0; % Dark Yellow
    0.75 0 0.75; % Dark Magenta
    0 0.75 0.75; % Dark Cyan
];
%% 1. Couleur / complexite de chaque event RACE
NRace = size(Race, 2);
n_active = sum(PCl_final, 1);           % nb d'assemblees actives par event

colors_asm = colors(1:NCl, :);          % votre palette existante
color_noise = [0.8 0.8 0.8];

event_color = zeros(NRace, 3);
is_complex  = false(1, NRace);

for j = 1:NRace
    active = find(PCl_final(:, j) == 1);
    if isempty(active)
        event_color(j, :) = color_noise;
    elseif numel(active) == 1
        event_color(j, :) = colors_asm(active, :);
    else
        event_color(j, :) = mean(colors_asm(active, :), 1);  % couleur melangee
        is_complex(j) = true;
    end
end

%% 2. Reduction de dimension (PCA -> t-SNE) sur le profil cellulaire complet
X = double(Race)';                       % NRace x NCell
nComp = min(50, size(X, 2));
[~, score] = pca(X, 'NumComponents', nComp);

rng(0);
Y_tsne = tsne(score, 'NumDimensions', 2, 'Perplexity', 30);

%% 3. Figure : points simples en rond, co-activations en losange
fig = figure('visible', 'off', 'Color', 'w', 'Position', [100 100 1200 550]);

simple_mask  = n_active <= 1;
complex_mask = n_active >= 2;

subplot(1,2,1); hold on;
scatter(score(simple_mask,1), score(simple_mask,2), 20, event_color(simple_mask,:), 'filled', 'MarkerFaceAlpha', 0.7);
scatter(score(complex_mask,1), score(complex_mask,2), 40, event_color(complex_mask,:), 'd', 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
xlabel('PC1'); ylabel('PC2'); title('PCA des RACE events'); box on; hold off;

subplot(1,2,2); hold on;
scatter(Y_tsne(simple_mask,1), Y_tsne(simple_mask,2), 20, event_color(simple_mask,:), 'filled', 'MarkerFaceAlpha', 0.7);
scatter(Y_tsne(complex_mask,1), Y_tsne(complex_mask,2), 40, event_color(complex_mask,:), 'd', 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
xlabel('t-SNE 1'); ylabel('t-SNE 2'); title('t-SNE des RACE events'); box on; hold off;

sgtitle('Rond = 0/1 assemblee active, Losange = co-activation (>=2 assemblees)');

%% 4. Export
namegraph = strcat(namefull, 'tsne_coactivation.png');
if isfolder(namefull)
    exportgraphics(fig, namegraph, 'Resolution', 200, 'BackgroundColor', 'white');
    close(fig);
end

%% 5. (Optionnel) Meme logique pour les CELLULES plutot que les events
% Utile pour voir si des cellules "promiscuous" (rejetees en 2.1 de votre
% script car actives dans 2 clusters) se placent bien entre deux
% assemblees dans l'espace reduit.
%
% Xc = double(Race);                     % NCell x NRace
% [~, scoreC] = pca(Xc, 'NumComponents', min(50, size(Xc,2)));
% Yc_tsne = tsne(scoreC, 'NumDimensions', 2, 'Perplexity', 30);
% cell_color = repmat(color_noise, NCell, 1);
% for i = 1:NCl
%     cell_color(assemblyortho{i}, :) = repmat(colors_asm(i,:), numel(assemblyortho{i}), 1);
% end
% figure; scatter(Yc_tsne(:,1), Yc_tsne(:,2), 15, cell_color, 'filled');
% title('t-SNE des cellules, colore par assemblee assignee');