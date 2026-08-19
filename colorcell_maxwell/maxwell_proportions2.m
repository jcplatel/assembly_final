function [Donnees_Clustering, r, g, b] = maxwell_proportions2(R_corr, G_corr, B_corr)

bg_R = prctile(R_corr(R_corr > 0), 5);
bg_G = prctile(G_corr(G_corr > 0), 5);
bg_B = prctile(B_corr(B_corr > 0), 5);

R_clean = max(0, R_corr - bg_R);
G_clean = max(0, G_corr - bg_G);
B_clean = max(0, B_corr - bg_B);

Rmax = prctile(R_clean(R_clean > 0), 99.9);
Gmax = prctile(G_clean(G_clean > 0), 99.9);
Bmax = prctile(B_clean(B_clean > 0), 99.9);

R_clean = R_clean ./ max(Rmax, eps);
G_clean = G_clean ./ max(Gmax, eps);
B_clean = B_clean ./ max(Bmax, eps);

R_clean = min(max(R_clean, 0), 1);
G_clean = min(max(G_clean, 0), 1);
B_clean = min(max(B_clean, 0), 1);

Somme = R_clean + G_clean + B_clean + eps;

r = R_clean ./ Somme;
g = G_clean ./ Somme;
b = B_clean ./ Somme;

Xtri = r + 0.5*b;
Ytri = (sqrt(3)/2) * b;

Donnees_Clustering = [Xtri, Ytri];

end