function dessiner_camembert(ax, X, labels, rgb)
% Camembert dont chaque part garde SA couleur, même quand certaines parts
% valent 0 (elles sont retirées avant le tracé).
%   labels = []  -> pourcentages par défaut de pie
nz = find(X > 0);
if isempty(nz), return; end
if isempty(labels)
    h = pie(ax, X(nz), ones(1, numel(nz)));
else
    h = pie(ax, X(nz), ones(1, numel(nz)), labels(nz));
end
for k = 1:numel(nz)
    h(2*k-1).FaceColor = rgb(nz(k), :);          % part (patch)
    h(2*k).FontSize    = 8;                      % étiquette (texte)
    h(2*k).Position    = h(2*k).Position * 1.1;  % éloigne un peu le texte
end
end
