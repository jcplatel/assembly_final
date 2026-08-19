function [roi_corrected, k_estime, b_estime, image_corrected] =...
    bleedthrough_correction2(roi_color, roi_ca, limit_max_k, image_color, image_ca, ax_target)
% BLEEDTHROUGH_CORRECTION Estime et corrige la fuite spectrale (GCaMP -> couleur)
% 
% Arguments d'entrée :
%   - roi_color, roi_ca : Vecteurs (Nx1) des intensités des cellules
%   - k                 : (Optionnel) Forcer la valeur de pente manuellement. Vide [] pour auto.
%   - graph             : (Optionnel) 1 pour afficher le graphique, 0 sinon.
%   - limit_max_k       : (Optionnel) Limite maximale biologique pour la pente k (ex: 0.40)
%   - image_color       : (Optionnel) L'image entière (matrice 2D) du canal couleur
%   - image_ca          : (Optionnel) L'image entière (matrice 2D) du canal GCaMP
%
% Arguments de sortie :
%   - roi_corrected     : Vecteur (Nx1) des intensités de cellules corrigées
%   - k_estime, b_estime: Les paramètres du fit (y = kx + b)
%   - image_corrected   : (Optionnel) L'image entière corrigée pixel par pixel

    % Gestion des arguments optionnels
    if nargin < 4
        image_color = []; 
        image_ca = [];
    end

    % 1. Filtrer les zéros et les valeurs aberrantes (NaN)
    idx_valid = (roi_ca > 0) & (roi_color > 0) & isfinite(roi_ca) & isfinite(roi_color);
    x_val = roi_ca(idx_valid);
    y_val = roi_color(idx_valid);
    
    if length(x_val) < 10
        roi_corrected = roi_color; k_estime = 0; b_estime = 0; 
        image_corrected = image_color;
        return;
    end
    
    % =====================================================================
    % 2. MÉTHODE GÉOMÉTRIQUE : L'ÉLASTIQUE INFÉRIEUR (Lower Convex Hull)
    % =====================================================================
    x_aug = [x_val; max(x_val)];
    y_aug = [y_val; max(y_val)*10]; 
    
    try
        k_hull = convhull(x_aug, y_aug);
        x_hull = x_aug(k_hull);
        y_hull = y_aug(k_hull);
        
        [~, idx_min_x] = min(x_hull);
        idx_max_x = max(x_val); 
        
        x_lower = [];
        y_lower = [];
        for i = 1:length(x_hull)
            if x_hull(i) >= x_hull(idx_min_x) && x_hull(i) <= idx_max_x && y_hull(i) <= max(y_val)
                x_lower = [x_lower; x_hull(i)];
                y_lower = [y_lower; y_hull(i)];
            end
        end
        
        if length(x_lower) < 2
            x_lower = [min(x_val); max(x_val)];
            y_lower = [min(y_val); min(y_val)];
        end
    catch
        x_lower = [min(x_val); max(x_val)];
        y_lower = [min(y_val); min(y_val)];
    end
    
    % =====================================================================
    % 3. FIT ROBUSTE SUR L'ÉLASTIQUE INFÉRIEUR
    % =====================================================================
    if length(x_lower) >= 2
        p = robustfit(x_lower, y_lower);
        b_estime = p(1);
        k_estime = p(2);
    else
        k_estime = 0; b_estime = min(y_val);
    end
    
    % =====================================================================
    % 4. SÉCURITÉS PHYSIQUES ABSOLUES
    % =====================================================================
    % if k_estime < 0
    %     k_estime = 0;
    % end
    % 
    % if k_estime > limit_max_k
    %     fprintf('Attention: k (%.3f) dépasse le max. Bridage à %.3f.\n', k_estime, limit_max_k);
    %     k_estime = limit_max_k;
    % end
    % 
    % y_virtuel = y_val - k_estime * x_val;
    % b_estime = prctile(y_virtuel, 5); 
    % 
    % if b_estime < 0, b_estime = 0; end
    % if b_estime > min(y_val), b_estime = min(y_val); end
    % 
    % fprintf('Fuite GCaMP (k) : %.3f | Bruit (b) : %.3f\n', k_estime, b_estime);
    if k_estime < 0
        k_estime = 0;
    end
    
    % On sauvegarde le 'b' trouvé par le fit naturel avant tout bridage
    b_fit_original = b_estime; 
    
    % Bridage de la pente (k)
    k_bridage_actif = false;
    if k_estime > limit_max_k
        fprintf('Attention: k (%.3f) dépasse le max. Bridage à %.3f.\n', k_estime, limit_max_k);
        k_estime = limit_max_k;
        k_bridage_actif = true;
    end
    
    % Recalcul du 'b' basé sur les points réels
    y_virtuel = y_val - k_estime * x_val;
    b_estime = prctile(y_virtuel, 5); 
    
    % BRIDAGE DE SÉCURITÉ SUR LE 'b' :
    % Si on a bridé 'k', le 'b' va avoir tendance à exploser pour compenser.
    % On s'assure qu'il ne dépasse pas le 'b' d'origine, ou un maximum absolu (ex: 5% du max de Y).
    if k_bridage_actif
        b_max_autorise = min([b_fit_original, max(y_val)*0.05]); 
        if b_estime > b_max_autorise
            b_estime = b_max_autorise;
        end
    end
    
    % Limites standards de sécurité basse
    if b_estime < 0, b_estime = 0; end
    if b_estime > min(y_val), b_estime = min(y_val); end
    % =====================================================================
    % 5. AFFICHAGE
    % =====================================================================
  if isempty(ax_target)
        % Comportement d'origine : on crée une nouvelle figure
        figure('Name', 'Correction Bleedthrough', 'Color', 'w');
        ax = gca;
    else
        % On utilise l'axe fourni par l'interface UI
        ax = ax_target;
    end
    
    % On trace tout en spécifiant explicitement 'ax'
    scatter(ax, roi_ca, roi_color, 30, [0.7 0.7 0.7], 'o'); hold(ax, 'on');
    
    x_line = linspace(0, max(roi_ca), 100);
    y_line = k_estime * x_line + b_estime;
    plot(ax, x_line, y_line, 'g-', 'LineWidth', 3);
    
    scatter(ax, x_lower, y_lower, 50, 'k', 'filled');
    
    % title(ax, 'Correction GCaMP (Lower Convex Hull)'); % Commenté car géré par l'UI principale
    legend(ax, 'Cellules', sprintf('Fit: y = %.2fx + %.2f', k_estime, b_estime), 'Points Convex Hull', 'Location', 'best');
    
    xlim(ax, [0 max(x_val)*1.05]);
    ylim(ax, [0 max(y_val)*1.05]);
    
    hold(ax, 'off');
    
    % =====================================================================
    % 6. CORRECTION FINALE (ROIs)
    % =====================================================================
    roi_corrected = roi_color - (k_estime * roi_ca + b_estime);
    roi_corrected(roi_corrected < 0) = 0;
    
    % =====================================================================
    % 7. CORRECTION FINALE (IMAGE ENTIÈRE)
    % =====================================================================
    if ~isempty(image_color) && ~isempty(image_ca)
        % On applique exactement la même équation mathématique aux pixels
        image_corrected = image_color - (k_estime * image_ca + b_estime);
        % Coup de hache pixel-wise pour éviter les valeurs de lumière négatives [cite:108]
        image_corrected(image_corrected < 0) = 0; 
    else
        image_corrected = [];
    end
end