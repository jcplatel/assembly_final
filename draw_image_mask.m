function draw_image_mask(img1, img2, img3, img4, Coomasqx, Coomasqy, colorcell, prctile_low, prctile_high, image_to_show, roi_show)
    figure('Visible', 'off')
    
    % --- Création de colormaps sur mesure ---
    map_length = 256;
    zeros_col = zeros(map_length, 1);
    lin_col = linspace(0, 1, map_length)';
    cmap_R = [lin_col, zeros_col, zeros_col]; 
    cmap_G = [zeros_col, lin_col, zeros_col]; 
    cmap_B = [zeros_col, zeros_col, lin_col]; 
    
    % --- Affichage UNIFIÉ des images ---
    if image_to_show == 0 
        img1_norm = mat2gray(img1, [prctile(img1, prctile_low, 'all'), prctile(img1, prctile_high, 'all')]); 
        img2_norm = mat2gray(img2, [prctile(img2, prctile_low, 'all'), prctile(img2, prctile_high, 'all')]);
        img3_norm = mat2gray(img3, [prctile(img3, prctile_low, 'all'), prctile(img3, prctile_high, 'all')]);
        img_rgb = cat(3, img1_norm, img2_norm, img3_norm);
        
        % Remplacement de imshow par l'affichage d'objet image natif
        h = image(img_rgb); 
        % title('RGB Image (1=R, 2=G, 3=B)');
    elseif image_to_show == 1
        imagesc(img1); colormap(cmap_R); % title('Image 1 (Red)');
        clim([prctile(img1, prctile_low, "all"), prctile(img1, prctile_high, "all")]);
    elseif image_to_show == 2
        imagesc(img2); colormap(cmap_G); % title('Image 2 (Green)');
        clim([prctile(img2, prctile_low, "all"), prctile(img2, prctile_high, "all")]);
    elseif image_to_show == 3
        imagesc(img3); colormap(cmap_B); % title('Image 3 (Blue)');
        clim([prctile(img3, prctile_low, "all"), prctile(img3, prctile_high, "all")]);
    elseif image_to_show == 4
        imagesc(img4); colormap("gray"); % title('Image 4 (Gray)');
        clim([prctile(img4, prctile_low, "all"), prctile(img4, prctile_high, "all")]);
    end
    
    hold on; % On active hold ON une fois que le fond est dessiné
    
    % --- Dessin Rapide et Exact basé sur les pixels de Suite2p ---
    if roi_show
        NCell = length(Coomasqx);
        
        bounds_x = cell(NCell, 1);
        bounds_y = cell(NCell, 1);
        max_pts = 0; 
        
        for n = 1:NCell
            x_all = Coomasqx{n}(:); 
            y_all = Coomasqy{n}(:);
            
            if length(x_all) >= 3 
                k = boundary(x_all, y_all, 1);
                xs = x_all(k);
                ys = y_all(k);
                
                bounds_x{n} = xs;
                bounds_y{n} = ys;
                
                if length(xs) > max_pts
                    max_pts = length(xs);
                end
            end
        end
        
        X_mat = NaN(max_pts, NCell);
        Y_mat = NaN(max_pts, NCell);
        
        for n = 1:NCell
            if ~isempty(bounds_x{n})
                len = length(bounds_x{n});
                X_mat(1:len, n) = bounds_x{n};
                Y_mat(1:len, n) = bounds_y{n};
            end
        end
        
        patch(X_mat, Y_mat, 'w', 'FaceAlpha', 0.1, 'EdgeColor', 'w', 'LineWidth', 1);

    end
    
    % --- VERROUILLAGE DÉFINITIF DES AXES ---
    axis image; % Force des pixels carrés
    axis off; % Enlève les graduations
    set(gca, 'YDir', 'reverse'); % Force l'axe Y vers le bas pour matcher avec les matrices
    
    hold off;
end