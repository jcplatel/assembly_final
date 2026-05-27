function draw_image_mask(img1, img2, img3, img4, Coomasqx, Coomasqy, colorcell, prctile_low, prctile_high, image_to_show, roi_show)

figure
hold on

% --- Création de colormaps sur mesure ---
map_length = 256;
zeros_col = zeros(map_length, 1);
lin_col = linspace(0, 1, map_length)';
cmap_R = [lin_col, zeros_col, zeros_col]; 
cmap_G = [zeros_col, lin_col, zeros_col]; 
cmap_B = [zeros_col, zeros_col, lin_col]; 

% --- Affichage des images ---
if image_to_show == 0 
    img1_norm = mat2gray(img1, [prctile(img1, prctile_low, 'all'), prctile(img1, prctile_high, 'all')]); 
    img2_norm = mat2gray(img2, [prctile(img2, prctile_low, 'all'), prctile(img2, prctile_high, 'all')]);
    img3_norm = mat2gray(img3, [prctile(img3, prctile_low, 'all'), prctile(img3, prctile_high, 'all')]);
    img_rgb = cat(3, img1_norm, img2_norm, img3_norm);
    imshow(img_rgb);
    title('RGB Image (1=R, 2=G, 3=B)');
elseif image_to_show == 1
    imagesc(img1); colormap(cmap_R); title('Image 1 (Red)');
    clim([prctile(img1, prctile_low, "all"), prctile(img1, prctile_high, "all")]);
elseif image_to_show == 2
    imagesc(img2); colormap(cmap_G); title('Image 2 (Green)');
    clim([prctile(img2, prctile_low, "all"), prctile(img2, prctile_high, "all")]);
elseif image_to_show == 3
    imagesc(img3); colormap(cmap_B); title('Image 3 (Blue)');
    clim([prctile(img3, prctile_low, "all"), prctile(img3, prctile_high, "all")]);
elseif image_to_show == 4
    imagesc(img4); colormap("gray"); title('Image 4 (Gray)');
    clim([prctile(img4, prctile_low, "all"), prctile(img4, prctile_high, "all")]);
end

% --- Dessin Rapide et Exact basé sur les pixels de Suite2p ---
if roi_show
    NCell = length(Coomasqx);
    
    bounds_x = cell(NCell, 1);
    bounds_y = cell(NCell, 1);
    max_pts = 0; 
    
    for n = 1:NCell
        % Récupération de tous les pixels de la ROI depuis stat.xpix / stat.ypix
        x_all = Coomasqx{n}(:); 
        y_all = Coomasqy{n}(:);
        
        if length(x_all) >= 3 
            % boundary(..., 1) donne l'enveloppe exacte des pixels de Suite2p sans dépasser.
            % C'est très naturel car ça suit les vrais contours d'intensité !
            k = boundary(x_all, y_all, 1);
            
            % On récupère le contour fermé
            xs = x_all(k);
            ys = y_all(k);
            
            bounds_x{n} = xs;
            bounds_y{n} = ys;
            
            if length(xs) > max_pts
                max_pts = length(xs);
            end
        end
    end
    
    % Création des matrices pour l'affichage vectorisé instantané
    X_mat = NaN(max_pts, NCell);
    Y_mat = NaN(max_pts, NCell);
    
    for n = 1:NCell
        if ~isempty(bounds_x{n})
            len = length(bounds_x{n});
            X_mat(1:len, n) = bounds_x{n};
            Y_mat(1:len, n) = bounds_y{n};
        end
    end
    
    % Dessiner toutes les ROIs
    patch(X_mat, Y_mat, 'w', 'FaceAlpha', 0.1, 'EdgeColor', 'w', 'LineWidth', 1);
end

axis off
axis image
hold off
end



% figure
% hold on
% 
% % --- Création de colormaps sur mesure ---
% map_length = 256;
% zeros_col = zeros(map_length, 1);
% lin_col = linspace(0, 1, map_length)';
% cmap_R = [lin_col, zeros_col, zeros_col]; 
% cmap_G = [zeros_col, lin_col, zeros_col]; 
% cmap_B = [zeros_col, zeros_col, lin_col]; 
% 
% % --- Affichage des images ---
% if image_to_show == 0 
%     img1_norm = mat2gray(img1, [prctile(img1, prctile_low, 'all'), prctile(img1, prctile_high, 'all')]); 
%     img2_norm = mat2gray(img2, [prctile(img2, prctile_low, 'all'), prctile(img2, prctile_high, 'all')]);
%     img3_norm = mat2gray(img3, [prctile(img3, prctile_low, 'all'), prctile(img3, prctile_high, 'all')]);
%     img_rgb = cat(3, img1_norm, img2_norm, img3_norm);
%     imshow(img_rgb);
%     title('RGB Image (1=R, 2=G, 3=B)');
% elseif image_to_show == 1
%     imagesc(img1); colormap(cmap_R); title('Image 1 (Red)');
%     clim([prctile(img1, prctile_low, "all"), prctile(img1, prctile_high, "all")]);
% elseif image_to_show == 2
%     imagesc(img2); colormap(cmap_G); title('Image 2 (Green)');
%     clim([prctile(img2, prctile_low, "all"), prctile(img2, prctile_high, "all")]);
% elseif image_to_show == 3
%     imagesc(img3); colormap(cmap_B); title('Image 3 (Blue)');
%     clim([prctile(img3, prctile_low, "all"), prctile(img3, prctile_high, "all")]);
% elseif image_to_show == 4
%     imagesc(img4); colormap("gray"); title('Image 4 (Gray)');
%     clim([prctile(img4, prctile_low, "all"), prctile(img4, prctile_high, "all")]);
% end
% 
% % --- Dessin Rapide et Lisse des contours ---
% if roi_show
%     NCell = length(maskx);
% 
%     % On prépare des cellules pour stocker les contours de chaque ROI
%     bounds_x = cell(NCell, 1);
%     bounds_y = cell(NCell, 1);
%     max_pts = 0; % Pour déterminer la taille de la matrice finale
% 
%     for n = 1:NCell
%         x = maskx{n}(:); 
%         y = masky{n}(:);
% 
%         if length(x) >= 4
%             % 1. Trouver les bords serrés
%             k = boundary(x, y, 0.8);
%             xb = x(k);
%             yb = y(k);
% 
%             % 2. Lissage périodique robuste (sans spline)
%             xb_ext = [xb(end-2:end-1); xb; xb(2:3)]; 
%             yb_ext = [yb(end-2:end-1); yb; yb(2:3)];
% 
%             % Interpolation linéaire pour densifier les points
%             t = 1:length(xb_ext);
%             t_iq = linspace(1, length(xb_ext), length(xb_ext)*5);
%             xb_dense = interp1(t, xb_ext, t_iq, 'linear');
%             yb_dense = interp1(t, yb_ext, t_iq, 'linear');
% 
%             % Lissage Gaussien
%             xs_smooth = smoothdata(xb_dense, 'gaussian', 15);
%             ys_smooth = smoothdata(yb_dense, 'gaussian', 15);
% 
%             % 3. Découpage propre de la boucle centrale
%             N_dense = length(xb_dense);
%             idx_start = floor(N_dense * (2 / length(xb_ext))) + 1;
%             idx_end = ceil(N_dense * ((length(xb_ext)-2) / length(xb_ext)));
% 
%             xs = xs_smooth(idx_start:idx_end)';
%             ys = ys_smooth(idx_start:idx_end)';
% 
%             % --- CORRECTION CRUCIALE ICI ---
%             % On force explicitement le dernier point à être identique au premier
%             % pour garantir que le polygone dessiné par patch() soit parfaitement fermé
%             xs = [xs; xs(1)];
%             ys = [ys; ys(1)];
% 
%             bounds_x{n} = xs;
%             bounds_y{n} = ys;
% 
%             if length(xs) > max_pts
%                 max_pts = length(xs);
%             end
%         end
%     end
% 
%     % 3. Création des matrices pour un dessin vectorisé (remplissage avec des NaN)
%     X_mat = NaN(max_pts, NCell);
%     Y_mat = NaN(max_pts, NCell);
% 
%     for n = 1:NCell
%         if ~isempty(bounds_x{n})
%             len = length(bounds_x{n});
%             X_mat(1:len, n) = bounds_x{n};
%             Y_mat(1:len, n) = bounds_y{n};
%         end
%     end
% 
%     % 4. Dessiner TOUTES les cellules en UN SEUL appel (Instantané)
%     patch(X_mat, Y_mat, 'w', 'FaceAlpha', 0.1, 'EdgeColor', 'w', 'LineWidth', 1);
% end
% 
% axis off
% axis image
% hold off
% end
% 
% % 
% % figure
% hold on
% 
% % --- Création de colormaps sur mesure (noir vers rouge, vert, ou bleu) ---
% map_length = 256;
% zeros_col = zeros(map_length, 1);
% lin_col = linspace(0, 1, map_length)';
% cmap_R = [lin_col, zeros_col, zeros_col]; % Noir -> Rouge
% cmap_G = [zeros_col, lin_col, zeros_col]; % Noir -> Vert
% cmap_B = [zeros_col, zeros_col, lin_col]; % Noir -> Bleu
% 
% if image_to_show == 0 
%     % Pour l'image RGB, on ajuste le contraste direct dans mat2gray
%     img1_norm = mat2gray(img1, [prctile(img1, prctile_low, 'all'), prctile(img1, prctile_high, 'all')]); 
%     img2_norm = mat2gray(img2, [prctile(img2, prctile_low, 'all'), prctile(img2, prctile_high, 'all')]);
%     img3_norm = mat2gray(img3, [prctile(img3, prctile_low, 'all'), prctile(img3, prctile_high, 'all')]);
% 
%     img_rgb = cat(3, img1_norm, img2_norm, img3_norm);
%     imshow(img_rgb);
%     title('RGB Image (1=R, 2=G, 3=B)');
% 
% elseif image_to_show == 1
%     imagesc(img1)
%     colormap(cmap_R)
%     clim([prctile(img1, prctile_low, "all"), prctile(img1, prctile_high, "all")]);
%     title('Image 1 (Red)');
% 
% elseif image_to_show == 2
%     imagesc(img2)
%     colormap(cmap_G)
%     clim([prctile(img2, prctile_low, "all"), prctile(img2, prctile_high, "all")]);
%     title('Image 2 (Green)');
% 
% elseif image_to_show == 3
%     imagesc(img3)
%     colormap(cmap_B)
%     clim([prctile(img3, prctile_low, "all"), prctile(img3, prctile_high, "all")]);
%     title('Image 3 (Blue)');
% 
% elseif image_to_show == 4
%     imagesc(img4)
%     colormap("gray")
%     clim([prctile(img4, prctile_low, "all"), prctile(img4, prctile_high, "all")]);
%     title('Image 4 (Gray)');
% end
% 
% % --- Dessin des contours (Patch) ---
% if roi_show
%     NCell = length(maskx);
%     for n = 1:NCell
%         x = cell2mat(maskx(n));
%         y = cell2mat(masky(n));
% 
%         % Optionnel : sécurité si un mask est vide ou n'a pas assez de points
%         if length(x) >= 3 
%             sorted_indices = convhull(x, y); 
%             x_hull = x(sorted_indices);
%             y_hull = y(sorted_indices);
%             patch(x_hull, y_hull, 'w', 'FaceAlpha', 0.1, 'EdgeColor', 'w', 'LineWidth', 1);
%         end
%     end
% end
% axis off
% axis image
% hold off
%end