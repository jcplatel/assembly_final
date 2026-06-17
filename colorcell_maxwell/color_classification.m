function idx_cluster = color_classification (Donnees_Clustering, NCell, r, g, b,ep)

X_tri = Donnees_Clustering(:,1);
Y_tri = Donnees_Clustering(:,2);
noms_theoriques = ["Rouge", "Vert", "Bleu", "Jaune", "Magenta", "Cyan"];
idx_cluster = zeros(NCell, 1);

for i = 1:NCell
    pR = r(i);
    pG = g(i);
    pB = b(i);
    
    if pR > pG && pR > pB
        % Le Rouge domine.
        if pG > pB
            % Tendance Jaune
            diff = pR - pG;
            if abs(diff - 0.4) < ep
                idx_cluster(i) = 8; % Frontière Noire
            elseif diff > 0.4 
                idx_cluster(i) = 1; % Rouge
            else
                idx_cluster(i) = 4; % Jaune
            end
        else
            % Tendance Magenta
            diff = pR - pB;
            if abs(diff - 0.4) < ep
                idx_cluster(i) = 8; % Frontière Noire
            elseif diff > 0.4 
                idx_cluster(i) = 1; % Rouge
            else
                idx_cluster(i) = 5; % Magenta
            end
        end
        
    elseif pG > pR && pG > pB
        % Le Vert domine.
        if pR > pB
            % Tendance Jaune
            diff = pG - pR;
            if abs(diff - 0.3) < ep
                idx_cluster(i) = 8; % Frontière Noire
            elseif diff > 0.3
                idx_cluster(i) = 2; % Vert
            else
                idx_cluster(i) = 4; % Jaune
            end
        else
            % Tendance Cyan
            diff = pG - pB;
            if abs(diff - 0.3) < ep
                idx_cluster(i) = 8; % Frontière Noire
            elseif diff > 0.3
                idx_cluster(i) = 2; % Vert
            else
                idx_cluster(i) = 6; % Cyan
            end
        end
        
    else
        % Le Bleu domine.
        if pR > pG
            % Tendance Magenta
            diff = pB - pR;
            if abs(diff - 0.3) < ep
                idx_cluster(i) = 8; % Frontière Noire
            elseif diff > 0.3
                idx_cluster(i) = 3; % Bleu
            else
                idx_cluster(i) = 5; % Magenta
            end
        else
            % Tendance Cyan
            diff = pB - pG;
            if abs(diff - 0.3) < ep
                idx_cluster(i) = 8; % Frontière Noire
            elseif diff > 0.3
                idx_cluster(i) = 3; % Bleu
            else
                idx_cluster(i) = 6; % Cyan
            end
        end
    end
    
        % =======================================================
        % 2. FORÇAGE DES FRONTIÈRES SUPPLÉMENTAIRES (NOIR)
        % =======================================================
        
        valeurs = sort([pR, pG, pB], 'descend');
        
        % A. Les frontières du "Y" (séparation entre primaires)
        % On veut cette ligne UNIQUEMENT quand on n'est PAS dans une couleur secondaire claire.
        % Dans une couleur secondaire, la 3ème couleur (la plus faible) est très basse (< 0.2).
        % Donc on n'applique cette frontière noire que si la couleur la plus faible est > 0.2
        % (ce qui correspond au milieu gris/blanc du triangle).
        
        if (valeurs(1) - valeurs(2)) < ep && valeurs(3) > 0.20
            idx_cluster(i) = 8;
        end
        
        % B. La ligne "horizontale" (séparation Jaune / Cyan-Magenta)
        % Elle sépare le "monde du bas" du "monde du haut"
        if (valeurs(2) - valeurs(3)) < ep
            if valeurs(1) < 0.7 
                idx_cluster(i) = 8;
            end
        end
    
end

end