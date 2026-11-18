function E = residus_paires(C, V, dist_um, pas_um, multiplicatif)
%RESIDUS_PAIRES  Retire de la similarité des paires ce qui ne dépend pas de
%   la lignée. AUCUNE couleur n'entre ici.
%
%     C_ij = f(distance_ij) + u_i + u_j [+ a_i a_j] + E_ij
%
%   f       : moyenne par tranche fine de distance (pas_um)
%   u_i     : niveau propre à chaque cellule (une cellule très active ou très
%             propre est plus synchrone avec tout le monde)
%   a_i a_j : couplage de chaque cellule au signal commun du champ
%             (multiplicatif ; utile pour la corrélation du ΔF/F)
%   Ajusté par moindres carrés sur les paires valides V (matrice logique).
%   E garde ce qui reste : c'est sur E qu'on compare les couleurs.

N = size(C, 1);
nbf = floor(400 / pas_um) + 1;
Bf = min(floor(dist_um / pas_um), nbf - 1) + 1;     % tranche fine de chaque paire (1..nbf)
Vf = double(V);
deg = sum(Vf, 2);
Cz = C;
Cz(~V) = 0;
u = zeros(N, 1);
a = zeros(N, 1);
if multiplicatif, a(:) = 0.1; end
A_u = diag(deg) + Vf + 1e-9 * eye(N);
nb_par_tranche = accumarray(Bf(V), 1, [nbf 1]);
for it = 1:15
    AA = a * a';
    R0 = Cz - u - u' - AA;
    f = accumarray(Bf(V), R0(V), [nbf 1]) ./ max(nb_par_tranche, 1);
    F = f(Bf);
    R1 = Cz - F - AA;
    R1(~V) = 0;
    u = A_u \ sum(R1, 2);
    if multiplicatif
        M = Cz - F - u - u';
        M(~V) = 0;
        for k = 1:3
            a = (M * a) ./ max(Vf * (a .^ 2), 1e-12);
        end
    end
end
E = C - f(Bf) - u - u' - a * a';
E(~V) = NaN;
end
