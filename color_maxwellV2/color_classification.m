function [idx_cluster, marge] = color_classification (Donnees_Clustering, NCell, r, g, b, ep, seuil_blanc, seuils)
    % CORRECTION du 11 nov. (demandée par l'utilisateur) : classe « blanc »
    % (7) pour les cellules qui ont les TROIS couleurs : la plus faible des
    % trois proportions dépasse seuil_blanc (défaut 0,20). Avant, ces
    % cellules étaient forcées dans une classe à deux couleurs (ex. 0,34 /
    % 0,33 / 0,33 -> magenta). seuil_blanc = Inf : comme avant (pas de blanc).
    if nargin < 7 || isempty(seuil_blanc), seuil_blanc = 0.20; end
    % 8e argument facultatif (17 nov.) : frontières « couleur pure / mélange » [rouge vert bleu]
    % (défaut [0.4 0.3 0.3] = comme avant ; MaxwellAuto utilise [0.4 0.3 0.4])
    if nargin < 8 || isempty(seuils), seuils = [0.4, 0.3, 0.3]; end
    % Ajout du 11 nov. (ne change RIEN à idx_cluster) : 2e sortie facultative
    % 'marge' = distance de la cellule à la frontière la plus proche de sa
    % classe, en proportion (0 = sur la frontière ; plus c'est grand, plus la
    % couleur est sûre). Frontières : seuil pur / mixte (0,4 pour le rouge,
    % 0,3 pour le vert et le bleu) et, pour une couleur mixte, l'écart entre
    % la 2e et la 3e proportion (jaune ou magenta ? etc.). NaN pour les
    % cellules sans signal (classe 8).
    marge = nan(NCell, 1);
    X_tri = Donnees_Clustering(:,1);
    Y_tri = Donnees_Clustering(:,2);
    idx_cluster = zeros(NCell, 1);
    
    for i = 1:NCell
        pR = r(i);
        pG = g(i);
        pB = b(i);
        
        % =======================================================
        % 0. GESTION DE SÉCURITÉ DES CELLULES REJETÉES (NaN)
        % =======================================================
        % Si une des valeurs est NaN, on classe immédiatement en Noir (8)
        if pR==0 && pG==0 && pB==0
            idx_cluster(i) = 8;
            continue; % Passe directement à l'itération 'i' suivante
        end
        % =======================================================
        
        if pR > pG && pR > pB
            % Le Rouge domine.
            if pG > pB
                % Tendance Jaune
                diff = pR - pG;
                if abs(diff - seuils(1)) < ep
                    idx_cluster(i) = 8; % Frontière Noire
                elseif diff > seuils(1) 
                    idx_cluster(i) = 1; % Rouge
                else
                    idx_cluster(i) = 4; % Jaune
                end
            else
                % Tendance Magenta
                diff = pR - pB;
                if abs(diff - seuils(1)) < ep
                    idx_cluster(i) = 8; % Frontière Noire
                elseif diff > seuils(1) 
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
                if abs(diff - seuils(2)) < ep
                    idx_cluster(i) = 8; % Frontière Noire
                elseif diff > seuils(2)
                    idx_cluster(i) = 2; % Vert
                else
                    idx_cluster(i) = 4; % Jaune
                end
            else
                % Tendance Cyan
                diff = pG - pB;
                if abs(diff - seuils(2)) < ep
                    idx_cluster(i) = 8; % Frontière Noire
                elseif diff > seuils(2)
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
                if abs(diff - seuils(3)) < ep
                    idx_cluster(i) = 8; % Frontière Noire
                elseif diff > seuils(3)
                    idx_cluster(i) = 3; % Bleu
                else
                    idx_cluster(i) = 5; % Magenta
                end
            else
                % Tendance Cyan
                diff = pB - pG;
                if abs(diff - seuils(3)) < ep
                    idx_cluster(i) = 8; % Frontière Noire
                elseif diff > seuils(3)
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
        if pR > pG && pR > pB, seuil_m = seuils(1); elseif pG > pR && pG > pB, seuil_m = seuils(2); else, seuil_m = seuils(3); end
        ecart = valeurs(1) - valeurs(2);
        if ecart > seuil_m
            marge(i) = ecart - seuil_m;
        else
            marge(i) = min(seuil_m - ecart, valeurs(2) - valeurs(3));
        end

        % Blanc : les trois couleurs présentes
        if valeurs(3) > seuil_blanc
            idx_cluster(i) = 7;
            marge(i) = valeurs(3) - seuil_blanc;
        elseif isfinite(seuil_blanc)
            marge(i) = min(marge(i), seuil_blanc - valeurs(3));
        end
        
        % A. Les frontières du "Y" (séparation entre primaires)
        if (valeurs(1) - valeurs(2)) < ep && valeurs(3) > 0.20
            idx_cluster(i) = 8;
        end
        
        % B. La ligne "horizontale" (séparation Jaune / Cyan-Magenta)
        if (valeurs(2) - valeurs(3)) < ep
            if valeurs(1) < 0.7 
                idx_cluster(i) = 8;
            end
        end
        
    end
end