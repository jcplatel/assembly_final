function M = CovarMjc(E)
    % E doit être la matrice (Cellules x SCE) ou (SCE x Cellules) selon ce que vous clusterisez
    
    % 1. Calcul vectorisé (équivalent à vos boucles mais instantané)
    % 'Rows','pairwise' gère les NaNs automatiquement s'il y en a
    M = corr(E, 'Rows', 'pairwise'); 
    
    % 2. SUPPRESSION DE LA DIAGONALE (Crucial pour vos shuffles !)
    % On met la diagonale à 0 pour que le clustering ne soit basé que 
    % sur les interactions ENTRE éléments, pas sur leur propre variance.
    n = size(M, 1);
    M(logical(eye(n))) = 0;
    
    % 3. Nettoyage final
    M(isnan(M)) = 0;
end