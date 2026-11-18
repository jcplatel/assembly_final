function [s, sCl, nSCE_cluster] = silhouette_clusters(M, IDX, NCl, MinSCE)
% Silhouette de chaque SCE sur la distance 1 - M, puis moyenne par cluster.
% C'est la SEULE définition de la silhouette dans le pipeline (données,
% surrogates, clusters finaux).
% Un cluster de moins de MinSCE SCE reçoit -1 (valeur minimale possible) :
% il ne peut jamais passer le seuil.

D_square = 1 - M;
D_square(logical(eye(size(D_square)))) = 0;   % diagonale à 0
D_square = (D_square + D_square') / 2;         % symétrie parfaite
D_vec = squareform(D_square);
s = silhouette([], IDX, D_vec);

nSCE_cluster = accumarray(IDX(:), 1, [NCl 1])';
sCl = -ones(1, NCl);
for i = 1:NCl
    if nSCE_cluster(i) >= MinSCE
        sCl(i) = mean(s(IDX == i));
    end
end
end
