function [roi_red_corrected, roi_green_corrected, roi_blue_corrected, id_low_intensity2] = ...
    interactive_threshold_selection(aligned_red, aligned_green, aligned_blue, calcium_norm, ...
                                    R, G, B, outline_gcampx, outline_gcampy, seuil_initial, id_low_intensity_prev)
    
    % =====================================================================
    % INITIALISATION DES VARIABLES
    % =====================================================================
    I_all = sqrt(R.^2 + G.^2 + B.^2); 
    NCell_total = length(I_all);
    
    manual_force_reject = false(NCell_total, 1);
    manual_force_keep = false(NCell_total, 1);
    
    % --- GESTION DE L'HISTORIQUE ET DU SEUIL INITIAL ---
    if nargin >= 11 && ~isempty(id_low_intensity_prev)
        % Si on charge un historique, on force le slider à 0
        seuil_choisi = 0; 
        idx_low_global = false(NCell_total, 1); % Le seuil global ne rejette rien au départ
        
        % On applique strictement l'historique
        manual_force_reject(id_low_intensity_prev) = true; 
    else
        % Sinon, comportement normal : on utilise le seuil fourni
        seuil_choisi = seuil_initial;
        idx_low_global = I_all < seuil_choisi;
    end
    
    % Pré-calcul des centres pour la sélection par rectangle
    cell_centers_x = zeros(NCell_total, 1);
    cell_centers_y = zeros(NCell_total, 1);
    for n = 1:NCell_total
        if n <= length(outline_gcampx) && n <= length(outline_gcampy) && ~isempty(outline_gcampx{n})
            cell_centers_x(n) = mean(outline_gcampx{n});
            cell_centers_y(n) = mean(outline_gcampy{n});
        end
    end
    
    zero_img = zeros(size(aligned_red));
    
    % =====================================================================
    % FENÊTRE ET AXE (Optimisés pour remplissage total)
    % =====================================================================
    screenSize = get(0, 'ScreenSize');
    figWidth = min(1600, screenSize(3) - 50); 
    figHeight = min(1000, screenSize(4) - 80);
    fig = uifigure('Name', 'Sélection interactive du seuil & Clics', 'Position', [25 40 figWidth figHeight]);
    
    axWidth = figWidth - 260;  
    axHeight = figHeight - 120;
    
    ax = uiaxes(fig, 'Position', [20, 100, axWidth, axHeight]);
    
    % Utilisation de imagesc pour forcer l'étirement
    h_img = imagesc(ax, cat(3, aligned_red, aligned_green, aligned_blue));
    
    title(ax, 'Rouge = Rejetée | Vert = Conservée | Clic = Inverser | Rectangle = Groupe', 'FontSize', 14);
    hold(ax, 'on');
    set(ax, 'YDir', 'reverse');
    
    % Configuration de l'axe pour un zoom parfait sans marge blanche
    axis(ax, 'off');
    axis(ax, 'normal'); 
    set(ax, 'DataAspectRatioMode', 'auto');
    set(ax, 'PlotBoxAspectRatioMode', 'auto');
    
    % =====================================================================
    % DESSIN DES ROIs (Contour lissé)
    % =====================================================================
    h_rois = gobjects(NCell_total, 1);
    for n = 1:NCell_total
        if n <= length(outline_gcampx) && n <= length(outline_gcampy)
            x_brut = outline_gcampx{n};
            y_brut = outline_gcampy{n};
            if length(x_brut) >= 3
                sorted_indices = convhull(x_brut, y_brut);
                x_contour = x_brut(sorted_indices);
                y_contour = y_brut(sorted_indices);
                
                % Ligne simple, HitTest = off
                h_rois(n) = plot(ax, x_contour, y_contour, 'Color', 'g', 'LineWidth', 1, ...
                                 'HitTest', 'off', 'PickableParts', 'none');
            end
        end
    end
    
    % =====================================================================
    % CONTRÔLES (SLIDER & PANNEAU)
    % =====================================================================
    lbl = uilabel(fig, 'Position', [100, 60, 600, 22], 'FontWeight', 'bold', 'FontSize', 14);
    
    % Le slider prend la valeur de seuil_choisi (qui vaut 0 en cas d'historique)
    sld = uislider(fig, 'Position', [100, 40, axWidth-100, 3], 'Limits', [0 1], 'Value', seuil_choisi);
    sld.ValueChangedFcn = @(sld, event) update_rois(sld.Value);
    
    pnlX = figWidth - 250;
    pnl = uipanel(fig, 'Position', [pnlX, figHeight-400, 220, 350], 'Title', 'Canaux & Luminosité');
    
    cb_red = uicheckbox(pnl, 'Position', [10, 300, 100, 22], 'Text', 'Rouge', 'Value', 1);
    sld_red = uislider(pnl, 'Position', [15, 280, 180, 3], 'Limits', [0 5], 'Value', 1);
    cb_green = uicheckbox(pnl, 'Position', [10, 230, 100, 22], 'Text', 'Vert', 'Value', 1);
    sld_green = uislider(pnl, 'Position', [15, 210, 180, 3], 'Limits', [0 5], 'Value', 1);
    cb_blue = uicheckbox(pnl, 'Position', [10, 160, 100, 22], 'Text', 'Bleu', 'Value', 1);
    sld_blue = uislider(pnl, 'Position', [15, 140, 180, 3], 'Limits', [0 5], 'Value', 1);
    cb_ca = uicheckbox(pnl, 'Position', [10, 90, 100, 22], 'Text', 'GCaMP', 'Value', 0);
    sld_ca = uislider(pnl, 'Position', [15, 70, 180, 3], 'Limits', [0 5], 'Value', 1);
    cb_rois = uicheckbox(pnl, 'Position', [10, 20, 180, 22], 'Text', 'Afficher les ROIs', 'Value', 1, 'FontWeight', 'bold');
    
    cb_red.ValueChangedFcn = @(cb,event) update_image_channels(); sld_red.ValueChangedFcn = @(sld,event) update_image_channels();
    cb_green.ValueChangedFcn = @(cb,event) update_image_channels(); sld_green.ValueChangedFcn = @(sld,event) update_image_channels();
    cb_blue.ValueChangedFcn = @(cb,event) update_image_channels(); sld_blue.ValueChangedFcn = @(sld,event) update_image_channels();
    cb_ca.ValueChangedFcn = @(cb,event) update_image_channels(); sld_ca.ValueChangedFcn = @(sld,event) update_image_channels();
    cb_rois.ValueChangedFcn = @(cb,event) toggle_rois_visibility(cb.Value);
    
    % =====================================================================
    % BOUTONS
    % =====================================================================
    uibutton(fig, 'Position', [pnlX, figHeight-470, 220, 50], 'Text', 'Sélection par Rectangle', ...
             'FontWeight', 'bold', 'ButtonPushedFcn', @(btn,event) trigger_rectangle_selection());
             
    uibutton(fig, 'Position', [pnlX, figHeight-550, 220, 60], 'Text', 'VALIDER', ...
             'FontWeight', 'bold', 'FontSize', 16, 'ButtonPushedFcn', @(btn,event) valider_callback());
    
    % =====================================================================
    % GESTION DES CLICS GLOBAUX SUR LA FENÊTRE
    % =====================================================================
    fig.WindowButtonDownFcn = @(src, event) handle_image_click();
    
    % =====================================================================
    % BOUCLE D'ATTENTE SÉCURISÉE
    % =====================================================================
    app_is_running = true;
    update_rois(seuil_choisi);
    
    while app_is_running && isgraphics(fig)
        uiwait(fig);
    end
    
    if isgraphics(fig)
        close(fig);
    end
    
    % =====================================================================
    % SORTIES FINALES
    % =====================================================================
    mask_low_intensity = (idx_low_global & ~manual_force_keep) | manual_force_reject;
    
    % --- FORMAT DE SORTIE ---
    % La sortie sera TOUJOURS une liste d'index de cellules, prête à être 
    % réinjectée en entrée la prochaine fois.
    id_low_intensity2 = find(mask_low_intensity); 
    % ------------------------
    
    roi_red_corrected = R; roi_green_corrected = G; roi_blue_corrected = B;
    
    % On nettoie les mauvaises cellules pour la suite du pipeline
    roi_red_corrected(mask_low_intensity) = NaN;
    roi_green_corrected(mask_low_intensity) = NaN;
    roi_blue_corrected(mask_low_intensity) = NaN;
    
    % =====================================================================
    % FONCTIONS IMBRIQUÉES (Callbacks)
    % =====================================================================
    function valider_callback()
        app_is_running = false; 
        uiresume(fig);
    end
    
    function handle_image_click()
        % On récupère les coordonnées du clic
        cp = ax.CurrentPoint;
        x_click = cp(1,1);
        y_click = cp(1,2);
        
        xlims = ax.XLim;
        ylims = ax.YLim;
        
        if x_click >= xlims(1) && x_click <= xlims(2) && y_click >= ylims(1) && y_click <= ylims(2)
            for i = 1:NCell_total
                if i <= length(outline_gcampx) && i <= length(outline_gcampy)
                    x_brut = outline_gcampx{i};
                    y_brut = outline_gcampy{i};
                    
                    if length(x_brut) >= 3
                        % On lisse le contour avant de chercher le clic dedans
                        sorted_indices = convhull(x_brut, y_brut);
                        x_lisse = x_brut(sorted_indices);
                        y_lisse = y_brut(sorted_indices);
                        
                        if inpolygon(x_click, y_click, x_lisse, y_lisse)
                            toggle_manual_state(i, true);
                            break; 
                        end
                    end
                end
            end
        end
    end
    
    function update_rois(current_threshold)
        seuil_choisi = current_threshold;
        idx_low_global = I_all < current_threshold;
        refresh_colors();
    end
    
    function toggle_manual_state(cell_id, update_vis)
        est_noire = (idx_low_global(cell_id) && ~manual_force_keep(cell_id)) || manual_force_reject(cell_id);
        if est_noire
            manual_force_reject(cell_id) = false; manual_force_keep(cell_id) = true;
        else
            manual_force_keep(cell_id) = false; manual_force_reject(cell_id) = true;
        end
        if update_vis
            refresh_colors();
        end
    end
    
    function trigger_rectangle_selection()
        % Désactive le clic global pour dessiner sans bug
        fig.WindowButtonDownFcn = '';
        
        rect = drawrectangle(ax, 'Color', 'y', 'FaceAlpha', 0.2);
        pos = rect.Position;
        x_min = pos(1); x_max = pos(1) + pos(3);
        y_min = pos(2); y_max = pos(2) + pos(4);
        
        for i = 1:NCell_total
            if cell_centers_x(i) >= x_min && cell_centers_x(i) <= x_max && ...
               cell_centers_y(i) >= y_min && cell_centers_y(i) <= y_max
                
                % On force le passage en "rejetée" (négative) sans bascule
                manual_force_keep(i) = false;
                manual_force_reject(i) = true;
            end
        end
        
        delete(rect);
        refresh_colors();
        
        % Réactive le clic global
        fig.WindowButtonDownFcn = @(src, event) handle_image_click();
    end
    
    function refresh_colors()
        idx_final = (idx_low_global & ~manual_force_keep) | manual_force_reject;
        for i = 1:NCell_total
            if isgraphics(h_rois(i))
                if idx_final(i)
                    h_rois(i).Color = [1 0 0];
                else
                    h_rois(i).Color = [0 1 0];
                end
            end
        end
        lbl.Text = sprintf('Seuil: %.2f | Rejetées: %d (Clics: %d sauvées, %d rejetées)', ...
                           seuil_choisi, sum(idx_final), sum(manual_force_keep), sum(manual_force_reject));
    end
    
    function update_image_channels()
        r_layer = zero_img; g_layer = zero_img; b_layer = zero_img;
        if cb_red.Value == 1,   r_layer = aligned_red * sld_red.Value;       end
        if cb_green.Value == 1, g_layer = aligned_green * sld_green.Value;   end
        if cb_blue.Value == 1,  b_layer = aligned_blue * sld_blue.Value;     end
        if cb_ca.Value == 1
            ca_boosted = calcium_norm * sld_ca.Value;
            r_layer = r_layer + ca_boosted; g_layer = g_layer + ca_boosted; b_layer = b_layer + ca_boosted;
        end
        r_layer(r_layer > 1) = 1; g_layer(g_layer > 1) = 1; b_layer(b_layer > 1) = 1;
        h_img.CData = cat(3, r_layer, g_layer, b_layer);
    end
    
    function toggle_rois_visibility(state)
        for i = 1:NCell_total
            if isgraphics(h_rois(i))
                if state == 1, h_rois(i).Visible = 'on'; else, h_rois(i).Visible = 'off'; end
            end
        end
    end
end