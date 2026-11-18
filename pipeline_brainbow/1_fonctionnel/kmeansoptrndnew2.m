function sCl = kmeansoptrndnew2(E, N, NCl, MinSCE, null_kmeans, online)
% Version 2 de kmeansoptrndnew : même procédure que kmeansoptnew2, sur des
% données mélangées. Changements :
%  - mêmes options k-means que pour les vraies données (Replicates = N,
%    OnlinePhase) : sinon les vraies données sont mieux optimisées et
%    obtiennent mécaniquement de meilleures silhouettes
%  - même règle MinSCE (sans elle, un singleton à s = 1 fixerait le seuil à 1)
%  - modèle nul 'curveball' : garde à la fois le nombre de SCE de chaque
%    cellule et la taille de chaque SCE. 'permutation_colonnes' = ancien
%    modèle (garde seulement la taille des SCE).
% Le générateur aléatoire est fixé par l'appelant (fixer_flux_aleatoire).
% online = false (facultatif) : sans phase en ligne, comme l'ancien
% kmeansoptrndnew (sert à retrouver_ancien).

if nargin < 6, online = true; end
if online, phase = 'on'; else, phase = 'off'; end

[Np, Ne] = size(E);

%% Randomization
switch null_kmeans
    case 'curveball'
        ERnd = double(curveball(E > 0));
    case 'permutation_colonnes'
        ERnd = zeros(Np, Ne);
        for i = 1:Ne
            ERnd(:, i) = E(randperm(Np), i);
        end
    otherwise
        error('null_kmeans doit valoir ''curveball'' ou ''permutation_colonnes''.');
end

%% Covariance matrix
M_shuff = CovarMjc(ERnd);

%% k-means (UseParallel off : cette fonction tourne déjà dans un parfor)
options = statset('UseParallel', false);
IDX = kmeans(M_shuff, NCl, 'Options', options, 'MaxIter', 1000, ...
    'OnlinePhase', phase, 'Replicates', N);
IDX = IDX';

[~, sCl] = silhouette_clusters(M_shuff, IDX, NCl, MinSCE);
end
