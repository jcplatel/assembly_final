function figure_clusters(A, Race, namefull)
% Figure clusters.png de SCE_clustering : matrice de covariance des SCE
% triés par cluster, et raster Race trié (cellules par assemblée, SCE par
% assemblées actives).
if A.NCl < 1, return; end
MSort = A.M(A.x2, A.x2);
fig = figure('visible', 'off');
subplot(1, 2, 1)
imagesc(MSort)
colormap jet
axis image
xlabel('sorted SCE #')
ylabel('sorted SCE #')
subplot(1, 2, 2)
imagesc(Race(A.x1, A.RList), [-1 1.2])
axis image
xlabel('sorted SCE #')
ylabel('sorted Cell #')
exportgraphics(fig, strcat(namefull, 'clusters.png'), 'Resolution', 300)
close(fig)
end
