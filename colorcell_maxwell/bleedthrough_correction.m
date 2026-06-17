function [roi_corrected, k_estime, b_estime] = bleedthrough_correction(roi_color, roi_ca, k,graph)


    % 1. Filtrer les zéros
    idx_valid = (roi_ca > 0) & (roi_color > 0);
    x_val = roi_ca(idx_valid);
    y_val = roi_color(idx_valid);

    % Sécurité si canal vide
    if length(x_val) < 10
        roi_corrected = roi_color; k_estime = 0; b_estime = 0; return;
    end

    % 2. Découper l'axe X en tranches plus larges pour garantir des points
    num_bins = 12; 
    edges = linspace(0, max(x_val), num_bins + 1);
    x_min = []; y_min = [];
    
    % 3. Chercher le MINIMUM ABSOLU dans chaque tranche
    for i = 1:num_bins
        in_bin = (x_val >= edges(i)) & (x_val < edges(i+1));
        if sum(in_bin) > 3 % Au moins 3 points pour avoir un "bord"
            y_in_bin = y_val(in_bin);
            x_in_bin = x_val(in_bin);
            
            % On prend le minimum absolu
            [min_y, idx_min] = min(y_in_bin);
            
            x_min = [x_min; x_in_bin(idx_min)];
            y_min = [y_min; min_y];
        end
    end

    % 4. Nettoyage des minimums trouvés
    % Si un minimum est en fait au milieu du nuage (ex: Y élevé alors qu'on est en bas),
    % on l'ignore car la tranche manquait de cellules "pures".
    % On garde uniquement la moitié inférieure des minimums trouvés.
    seuil_haut = median(y_min);
    idx_vrais_bas = y_min <= seuil_haut;
    
    x_fit = x_min(idx_vrais_bas);
    y_fit = y_min(idx_vrais_bas);

    % 5. Fit linéaire robuste (y = k*x + b)
    if length(x_fit) >= 2
        b_k = robustfit(x_fit, y_fit); 
        b_estime = b_k(1); 
        k_estime = b_k(2);
    else
        k_estime = 0; b_estime = 0;
    end
    
    % Option de repli : Si le fit donne un k négatif (impossible physiquement), on le met à 0
    if k_estime < 0
        k_estime = 0; 
        b_estime = min(y_min); % Le fond devient juste le min global
    end
    
    fprintf('Fuite GCaMP (k) : %.3f | Bruit (b) : %.3f\n', k_estime, b_estime);

    % 6. Affichage
    if graph == 1
        figure('Name', 'Correction Bleedthrough', 'Color', 'w');
        scatter(roi_ca, roi_color, 30, [0.7 0.7 0.7], 'o'); hold on;
        scatter(x_fit, y_fit, 80, 'k', 'filled'); % Affiche uniquement les points gardés
        
        x_line = linspace(0, max(roi_ca), 100);
        plot(x_line, k_estime * x_line + b_estime, 'g-', 'LineWidth', 3);
        title('Correction GCaMP ' );
        legend('Données', 'Points du plancher', sprintf('Fit: y = %.2fx + %.2f', k_estime, b_estime), 'Location', 'best');
        hold off;
    end
    % 7. Correction
    if nargin==3 
        k_estime=k;
        roi_corrected = roi_color - (k_estime * roi_ca + b_estime);
        roi_corrected(roi_corrected < 0) = 0;
    else
        roi_corrected = roi_color - (k_estime * roi_ca + b_estime);
        roi_corrected(roi_corrected < 0) = 0;
    end
end