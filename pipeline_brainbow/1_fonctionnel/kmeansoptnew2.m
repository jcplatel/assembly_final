function [IDXs,sCl,M,S,NClini,sumd,nSCE_cluster] = kmeansoptnew2(E,N,type,NClini,Mcons,MinSCE,graine,online)
% Version 2 de kmeansoptnew : on clusterise les SCE sur la base des
% cellules qui y participent.
%   E : cellules x SCE ; N : nombre de réplicats k-means
% Changements :
%  - MinSCE : un cluster de moins de MinSCE SCE reçoit la silhouette -1 et
%    n'est donc jamais retenu. (MATLAB donne s = 1 à un SCE seul dans son
%    cluster : avant, les singletons passaient toujours le seuil.)
%  - générateur aléatoire fixé (graine) : résultat reproductible
%  - la silhouette (distance 1 - M, moyenne par cluster, silhouette_clusters)
%    est LA définition utilisée partout dans le pipeline
%  - sorties en plus : nSCE_cluster (nombre de SCE par cluster)
%  - online (facultatif, défaut true) : phase en ligne de kmeans (lente ;
%    doit être la même pour les données et les surrogates)

if nargin < 6, MinSCE = 3; end
if nargin < 7, graine = 1; end
if nargin < 8, online = true; end
if online, phase = 'on'; else, phase = 'off'; end
Ne = size(E, 2);

%% Covariance matrix
if strcmp(type, 'var')
    if isempty(Mcons)
        M = CovarMjc(E);
    else
        M = Mcons;
    end
else
    error('kmeansoptnew2 : type doit valoir ''var''.');
end

NCl = NClini;
options = statset('UseParallel', true, 'UseSubstreams', true, ...
    'Streams', RandStream('mlfg6331_64', 'Seed', graine));
[IDX, ~, sumd] = kmeans(M, NCl, 'Options', options, 'MaxIter', 1000, ...
    'OnlinePhase', phase, 'Replicates', N);   % k-means sur les lignes de la matrice de covariance
IDX = IDX';

[s, sCl, nSCE_cluster] = silhouette_clusters(M, IDX, NCl, MinSCE);
S = mean(s);

%% sort RACE/silhouette of best cluster
[sCl, xCl] = sort(sCl, 'descend');
IDXs = zeros(1, Ne);
for i = 1:NCl
    IDXs = IDXs + (IDX == xCl(i)) * i;
end
nSCE_cluster = nSCE_cluster(xCl);
sumd = sumd(xCl);
end
