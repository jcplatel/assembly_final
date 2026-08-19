function colorcell_updated = verify_colorcell_gui(aligned_red, aligned_green, aligned_blue, calcium_norm, ...
                                                  outline_gcampx, outline_gcampy, colorcell)
% VERIFY_COLORCELL_GUI Inspecte visuellement et corrige les couleurs des ROIs
%
% Entrées : Images, Contours, Vecteur colorcell (1 à 8)
% Sortie  : colorcell_updated (mis à jour avec vos clics)

    NCell_total = length(colorcell);
    zero_img = zeros(size(aligned_red));
    colorcell_updated = colorcell;
    
    % Palette de couleurs (1 à 8)
    cmap = {
        [1 0 0],       % 1 = Rouge
        [0 1 0],       % 2 = Vert
        [0.2 0.5 1],   % 3 = Bleu
        [1 1 0],       % 4 = Jaune
        [1 0 1],       % 5 = Magenta
        [0 1 1],       % 6 = Cyan
        [1 1 1],       % 7 = Blanc
        [0.5 0.5 0.5]  % 8 = Noir (Affiché en gris)
    };
    noms_couleurs = {'Rouge','Vert','Bleu','Jaune','Magenta','Cyan','Blanc','Noir'};

    % =====================================================================
    % CRÉATION DE LA FENÊTRE
    % =====================================================================
    screenSize = get(0, 'ScreenSize');
    figWidth = min(1200, screenSize(3) - 100);
    figHeight = min(900, screenSize(4) - 100);
    fig = uifigure('Name', 'Vérification & Correction (Colorcell)', 'Position', [50 50 figWidth figHeight]);
    
    axWidth = figWidth - 300; 
    axHeight = figHeight - 100;
    ax = uiaxes(fig, 'Position', [20, 50, axWidth, axHeight]);
    
    h_img = imshow(cat(3, aligned_red, aligned_green, aligned_blue), 'Parent', ax);
    title(ax, 'Vérification : Clic = Forcer en Noir / Annuler', 'FontSize', 14, 'FontWeight', 'bold');
    hold(ax, 'on');
    set(ax, 'YDir', 'reverse');
    
    % =====================================================================
    % DESSIN DES ROIs ET DES NUMÉROS (AVEC CLICS)
    % =====================================================================
    h_rois = gobjects(NCell_total, 1);
    h_texts = gobjects(NCell_total, 1);
    
    for n = 1:NCell_total
        if n <= length(outline_gcampx) && n <= length(outline_gcampy) && ~isnan(colorcell(n))
            x = outline_gcampx{n};
            y = outline_gcampy{n};
            
            if length(x) >= 3
                sorted_indices = convhull(x, y);
                x_contour = x(sorted_indices);
                y_contour = y(sorted_indices);
                
                c_idx = round(colorcell(n));
                if c_idx >= 1 && c_idx <= 8
                    couleur_roi = cmap{c_idx};
                else
                    couleur_roi = [1 1 1]; 
                end
                
                h_rois(n) = plot(ax, x_contour, y_contour, 'Color', couleur_roi, 'LineWidth', 1.5, ...
                                 'HitTest', 'on', 'PickableParts', 'all');
                             
                h_rois(n).ButtonDownFcn = @(src, event) toggle_cell_black(n);
                
                cx = mean(x_contour);
                cy = mean(y_contour);
                h_texts(n) = text(ax, cx, cy, num2str(n), ...
                    'Color', 'w', 'FontSize', 8, 'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center', 'Visible', 'off', ...
                    'HitTest', 'off'); 
            end
        end
    end
    
    % =====================================================================
    % PANNEAU 1 : CANAUX & LUMINOSITÉ
    % =====================================================================
    pnlX = figWidth - 250;
    pnl1 = uipanel(fig, 'Position', [pnlX, figHeight-380, 220, 330], 'Title', 'Canaux & Luminosité');
    
    cb_red = uicheckbox(pnl1, 'Position', [10, 270, 100, 22], 'Text', 'Rouge', 'Value', 1);
    sld_red = uislider(pnl1, 'Position', [15, 250, 180, 3], 'Limits', [0 5], 'Value', 1);
    
    cb_green = uicheckbox(pnl1, 'Position', [10, 200, 100, 22], 'Text', 'Vert', 'Value', 1);
    sld_green = uislider(pnl1, 'Position', [15, 180, 180, 3], 'Limits', [0 5], 'Value', 1);
    
    cb_blue = uicheckbox(pnl1, 'Position', [10, 130, 100, 22], 'Text', 'Bleu', 'Value', 1);
    sld_blue = uislider(pnl1, 'Position', [15, 110, 180, 3], 'Limits', [0 5], 'Value', 1);
    
    cb_ca = uicheckbox(pnl1, 'Position', [10, 60, 100, 22], 'Text', 'GCaMP', 'Value', 0);
    sld_ca = uislider(pnl1, 'Position', [15, 40, 180, 3], 'Limits', [0 5], 'Value', 1);

    % =====================================================================
    % PANNEAU 2 : FILTRES ROIs (PAR COULEUR)
    % =====================================================================
    pnl2 = uipanel(fig, 'Position', [pnlX, figHeight-660, 220, 260], 'Title', 'Affichage des ROIs');
    
    % Création des 8 cases à cocher en 2 colonnes
    cb_colors = gobjects(8, 1);
    y_pos = [190, 155, 120, 85];
    x_pos = [10, 110];
    
    for i = 1:8
        col = ceil(i/4);       % Colonne 1 ou 2
        row = mod(i-1, 4) + 1; % Ligne 1 à 4
        
        cb_colors(i) = uicheckbox(pnl2, 'Position', [x_pos(col), y_pos(row), 90, 22], ...
                                  'Text', noms_couleurs{i}, 'Value', 1);
        cb_colors(i).ValueChangedFcn = @(cb,event) update_rois_visibility();
    end
    
    % Checkbox globale pour les numéros
    cb_nums = uicheckbox(pnl2, 'Position', [10, 20, 180, 22], 'Text', 'Afficher les Numéros', ...
                         'Value', 0, 'FontWeight', 'bold');
    cb_nums.ValueChangedFcn = @(cb,event) update_rois_visibility();
    
    % =====================================================================
    % LIAISONS & BOUTON
    % =====================================================================
    cb_red.ValueChangedFcn = @(cb,event) update_image_channels();
    sld_red.ValueChangedFcn = @(sld,event) update_image_channels();
    cb_green.ValueChangedFcn = @(cb,event) update_image_channels();
    sld_green.ValueChangedFcn = @(sld,event) update_image_channels();
    cb_blue.ValueChangedFcn = @(cb,event) update_image_channels();
    sld_blue.ValueChangedFcn = @(sld,event) update_image_channels();
    cb_ca.ValueChangedFcn = @(cb,event) update_image_channels();
    sld_ca.ValueChangedFcn = @(sld,event) update_image_channels();

    uibutton(fig, 'Position', [pnlX, figHeight-750, 220, 60], 'Text', 'VALIDER ET SAUVER', ...
             'FontWeight', 'bold', 'FontSize', 16, 'ButtonPushedFcn', @(btn,event) uiresume(fig));

    % Mise en attente pour récupérer le colorcell_updated final
    uiwait(fig);
    if isgraphics(fig)
        close(fig);
    end

    % --- Fonctions imbriquées ---
    function toggle_cell_black(cell_id)
        if colorcell_updated(cell_id) == 8 && colorcell(cell_id) ~= 8
            colorcell_updated(cell_id) = colorcell(cell_id);
        else
            colorcell_updated(cell_id) = 8;
        end
        
        c_idx = round(colorcell_updated(cell_id));
        if c_idx >= 1 && c_idx <= 8
            h_rois(cell_id).Color = cmap{c_idx};
        end
        % Met à jour la visibilité au cas où on transformerait la ROI 
        % en une couleur actuellement masquée
        update_rois_visibility(); 
    end

    function update_rois_visibility()
        show_nums = cb_nums.Value;
        
        for i = 1:NCell_total
            if isgraphics(h_rois(i))
                % Vérifie la couleur ACTUELLE de la cellule
                c_idx = round(colorcell_updated(i));
                
                % Détermine si cette couleur doit être affichée
                if c_idx >= 1 && c_idx <= 8
                    is_visible = cb_colors(c_idx).Value;
                else
                    is_visible = 1;
                end
                
                % Applique la visibilité
                if is_visible
                    h_rois(i).Visible = 'on';
                    if show_nums
                        h_texts(i).Visible = 'on';
                    else
                        h_texts(i).Visible = 'off';
                    end
                else
                    h_rois(i).Visible = 'off';
                    h_texts(i).Visible = 'off';
                end
            end
        end
    end

    function update_image_channels()
        r_layer = zero_img; g_layer = zero_img; b_layer = zero_img;
        
        if cb_red.Value == 1,   r_layer = aligned_red * sld_red.Value;       end
        if cb_green.Value == 1, g_layer = aligned_green * sld_green.Value;   end
        if cb_blue.Value == 1,  b_layer = aligned_blue * sld_blue.Value;     end
        
        if cb_ca.Value == 1
            ca_boosted = calcium_norm * sld_ca.Value;
            r_layer = r_layer + ca_boosted;
            g_layer = g_layer + ca_boosted;
            b_layer = b_layer + ca_boosted;
        end
        
        r_layer(r_layer > 1) = 1; 
        g_layer(g_layer > 1) = 1;
        b_layer(b_layer > 1) = 1;
        
        h_img.CData = cat(3, r_layer, g_layer, b_layer);
    end
end
