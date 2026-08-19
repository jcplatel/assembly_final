function [k_estime, b_estime, x_val, y_val, x_lower, y_lower] = ...
    estimate_bleedthrough_params(roi_color, roi_ca, limit_max_k)
    % 1. Filtrer les zéros et les valeurs aberrantes (NaN)
    idx_valid = (roi_ca > 0) & (roi_color > 0) & isfinite(roi_ca) & isfinite(roi_color);
    x_val = roi_ca(idx_valid);
    y_val = roi_color(idx_valid);
    
    if length(x_val) < 10
        k_estime = 0; b_estime = 0; x_lower = []; y_lower = []; return;
    end
    
    % 2. MÉTHODE GÉOMÉTRIQUE : L'ÉLASTIQUE INFÉRIEUR
    x_aug = [x_val; max(x_val)];
    y_aug = [y_val; max(y_val)*10]; 
    try
        k_hull = convhull(x_aug, y_aug);
        x_hull = x_aug(k_hull);
        y_hull = y_aug(k_hull);
        [~, idx_min_x] = min(x_hull);
        idx_max_x = max(x_val); 
        x_lower = []; y_lower = [];
        for i = 1:length(x_hull)
            if x_hull(i) >= x_hull(idx_min_x) && x_hull(i) <= idx_max_x && y_hull(i) <= max(y_val)
                x_lower = [x_lower; x_hull(i)]; y_lower = [y_lower; y_hull(i)];
            end
        end
        if length(x_lower) < 2
            x_lower = [min(x_val); max(x_val)]; y_lower = [min(y_val); min(y_val)];
        end
    catch
        x_lower = [min(x_val); max(x_val)]; y_lower = [min(y_val); min(y_val)];
    end
    
    % 3. FIT ROBUSTE SUR L'ÉLASTIQUE INFÉRIEUR
    if length(x_lower) >= 2
        p = robustfit(x_lower, y_lower);
        b_estime = p(1); k_estime = p(2);
    else
        k_estime = 0; b_estime = min(y_val);
    end
    
    % 4. SÉCURITÉS PHYSIQUES ABSOLUES
    if k_estime < 0, k_estime = 0; end
    b_fit_original = b_estime; 
    k_bridage_actif = false;
    
    if k_estime > limit_max_k
        k_estime = limit_max_k;
        k_bridage_actif = true;
    end
    
    y_virtuel = y_val - k_estime * x_val;
    b_estime = prctile(y_virtuel, 5); 
    
    if k_bridage_actif
        b_max_autorise = min([b_fit_original, max(y_val)*0.05]); 
        if b_estime > b_max_autorise, b_estime = b_max_autorise; end
    end
    
    if b_estime < 0, b_estime = 0; end
    if b_estime > min(y_val), b_estime = min(y_val); end
end