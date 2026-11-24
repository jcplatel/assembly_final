function [Donnees_Clustering, r, g, b, info] = maxwell_proportions(R_corr, G_corr, B_corr, methode)
% CORRECTION du 11 nov. (demandée par l'utilisateur). Avant, l'échelle de
% chaque canal venait de brainbow_image_processing : le 99,9e percentile des
% PIXELS de l'image (262 pixels, soit 2 ou 3 cellules). Une seule cellule
% très brillante écrasait tout son canal et changeait la classe des autres.
% Maintenant ('roi95', par défaut) : après le retrait du fond, chaque canal
% est divisé par le 95e percentile de l'intensité des CELLULES qui ont du
% signal dans ce canal (une vingtaine de cellules sur 500 : une cellule
% isolée ne change plus rien). Cette division annule exactement l'ancienne
% échelle, qui était un simple facteur.
%   methode : 'roi95' (défaut) ou 'ancien' (calcul d'avant, pour comparer)
%   info    : fond et échelle de chaque canal [R G B], somme des trois
%             canaux par cellule (intensité totale : faible = couleur peu sûre)
% Attention : les proportions changent, donc aussi les classes ; comparer
% sur une session avant de tout refaire (confusion_matrix).

    if nargin < 4 || isempty(methode), methode = 'roi95'; end

    % 1. Soustraction du Bruit de Fond Constant 

    bg_R = prctile(R_corr(R_corr > 0), 5);
    bg_G = prctile(G_corr(G_corr > 0), 5);
    bg_B = prctile(B_corr(B_corr > 0), 5);

    % On soustrait ce fond constant et on coupe strictement à 0
    R_clean = max(0, R_corr - bg_R);
    G_clean = max(0, G_corr - bg_G);
    B_clean = max(0, B_corr - bg_B);

    % 1 bis. Échelle de chaque canal
    echelle = [1 1 1];
    if strcmpi(methode, 'roi95')
        echelle = [ref95(R_clean), ref95(G_clean), ref95(B_clean)];
        R_clean = R_clean / echelle(1);
        G_clean = G_clean / echelle(2);
        B_clean = B_clean / echelle(3);
    end

    % 2. Calcul des Proportions de Maxwell
    Somme = R_clean + G_clean + B_clean + 1e-6; 
    
    r = R_clean ./ Somme;
    g = G_clean ./ Somme;
    b = B_clean ./ Somme;

    % 3. Projection géométrique de Maxwell
    Xtri = r + b * 0.5;
    Ytri = b * (sqrt(3)/2);
    
    Donnees_Clustering = [Xtri, Ytri];
    info = struct('methode', methode, 'fond', [bg_R, bg_G, bg_B], 'echelle', echelle, 'somme', Somme);

end

function v = ref95(x)
    % 95e percentile des cellules qui ont du signal ; 1 s'il n'y en a pas assez
    x = x(x > 0 & isfinite(x));
    if numel(x) < 20
        v = 1;
    else
        v = prctile(x, 95);
        if ~(v > 0), v = 1; end
    end
end
