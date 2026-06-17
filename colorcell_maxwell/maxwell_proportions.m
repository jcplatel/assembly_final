function [Donnees_Clustering, r, g, b] = maxwell_proportions(R_corr, G_corr, B_corr)

    % 1. Soustraction du Bruit de Fond Constant 

    bg_R = prctile(R_corr(R_corr > 0), 5);
    bg_G = prctile(G_corr(G_corr > 0), 5);
    bg_B = prctile(B_corr(B_corr > 0), 5);

    % On soustrait ce fond constant et on coupe strictement à 0
    R_clean = max(0, R_corr - bg_R);
    G_clean = max(0, G_corr - bg_G);
    B_clean = max(0, B_corr - bg_B);

    % 2. Calcul des Proportions de Maxwell
    Somme = R_clean + G_clean + B_clean + 1e-6; 
    
    r = R_clean ./ Somme;
    g = G_clean ./ Somme;
    b = B_clean ./ Somme;

    % 3. Projection géométrique de Maxwell
    Xtri = r + b * 0.5;
    Ytri = b * (sqrt(3)/2);
    
    Donnees_Clustering = [Xtri, Ytri];

end