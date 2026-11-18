function P = coparticipation_sce(Race)
%COPARTICIPATION_SCE  Co-participation de chaque paire de cellules aux SCE
%   (8 nov., demandé par l'utilisateur : « est-ce que les cellules de même
%   couleur participent plus aux mêmes SCE ? »).
%   Coefficient phi = corrélation de Pearson des deux lignes 0/1 de Race
%   (cellules x SCE). Tient compte du nombre de SCE de chaque cellule : deux
%   cellules très actives n'ont pas un phi élevé du seul fait de leur
%   activité. Même mesure que l'arbre 'sce' de l'étape 3c (distance = 1 - phi).
%   Cellule dans aucun SCE (ou dans tous) : NaN ; diagonale : NaN.
%   AUCUNE couleur n'entre ici.
X = double(Race > 0);
X = X - mean(X, 2);
s = sqrt(sum(X .^ 2, 2));
P = (X * X') ./ (s * s');
P(s == 0, :) = NaN;
P(:, s == 0) = NaN;
P(logical(eye(size(P, 1)))) = NaN;
end
