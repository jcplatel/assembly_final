function colorcell_updated = interactive_color_correction(aligned_red, aligned_green, aligned_blue, calcium_norm, colorcell_initial, Coomasqx, Coomasqy)
    
    % =====================================================================
    % INITIALISATION
    % =====================================================================
    NCell_total = length(colorcell_initial);
    colorcell_updated = colorcell_initial(:); 
    
    % Dictionnaire des couleurs
    cmap_classes = [
        1 0 0;      % 1. Rouge
        0 1 0;      % 2. Vert
        0 0 1;      % 3. Bleu
        1 1 0;      % 4. Jaune
        1 0 1;      % 5. Magenta
        0 1 1;      % 6. Cyan
        1 1 1;      % 7. Blanc
        0.5 0.5 0.5 % 8. Noir (Gris pour l'interface)
    ];
    
    noms_classes = {'1. Rouge', '2. Vert', '3. Bleu', '4. Jaune', '5. Magenta', '6. Cyan', '7. Blanc', '8. Noir'};
    couleurs_txt = {'r', [0 0.6 0], 'b', [0.8 0.8 0], 'm', 'c', 'k', [0.4 0.4 0.4]};
    current_selected_class = 7; % Par défaut : Blanc
    zero_img = zeros(size(aligned_red));
    
    % =====================================================================
    % FENÊTRE ET AXE
    % =====================================================================
    screenSize = get(0, 'ScreenSize');
    figWidth = min(1600, screenSize(3) - 50); 
    % On limite la hauteur à 750 max pour les petits écrans
    figHeight = min(750, screenSize(4) - 80);
    fig = uifigure('Name', 'Correction Interactive', 'Position', [25 40 figWidth figHeight]);
    
    axWidth = figWidth - 250;  
    axHeight = figHeight - 40;
    
    ax = uiaxes(fig, 'Position', [20, 20, axWidth, axHeight]);
    h_img = imagesc(ax, cat(3, aligned_red, aligned_green, aligned_blue));
    
    title(ax, 'Cochez/Décochez les couleurs à droite | Cliquez sur une cellule pour la corriger', 'FontSize', 14);
    hold(ax, 'on');
    set(ax, 'YDir', 'reverse');
    axis(ax, 'off');
    axis(ax, 'image'); 
    
    % =====================================================================
    % DESSIN DES ROIs (Basé sur Coomasqx/y avec convhull)
    % =====================================================================
    h_rois = gobjects(NCell_total, 1);
    h_markers = gobjects(NCell_total, 1);
    
    for n = 1:NCell_total
        if n <= length(Coomasqx) && n <= length(Coomasqy)
            x_all = Coomasqx{n}(:); 
            y_all = Coomasqy{n}(:);
            
            if length(x_all) >= 3
                k = convhull(x_all, y_all);
                xs = x_all(k);
                ys = y_all(k);
                
                c_idx = colorcell_updated(n);
                if c_idx < 1 || c_idx > 8 || isnan(c_idx), c_idx = 8; end 
                
                h_rois(n) = patch(ax, 'XData', xs, 'YData', ys, ...
                                 'FaceColor', cmap_classes(c_idx, :), 'FaceAlpha', 0.1, ...
                                 'EdgeColor', cmap_classes(c_idx, :), 'LineWidth', 1, ...
                                 'HitTest', 'off', 'PickableParts', 'none');
            end
        end
    end
    
    % =====================================================================
    % MISE EN PAGE VERTICALE ULTRA COMPACTE (X = panneau droite)
    % =====================================================================
    pnlX = figWidth - 220;
    
    % --- PANNEAU 3 : LUMINOSITÉ (EN HAUT) ---
    pnl_img = uipanel(fig, 'Position', [pnlX, figHeight-170, 200, 160], 'Title', 'Luminosité');
    
    % On réduit l'espacement entre sliders
    y_pos = 115; esp = 35;
    cb_red = uicheckbox(pnl_img, 'Position', [5, y_pos, 60, 20], 'Text', 'R', 'Value', 1);
    sld_red = uislider(pnl_img, 'Position', [65, y_pos+10, 120, 3], 'Limits', [0 5], 'Value', 1);
    
    y_pos = y_pos - esp;
    cb_green = uicheckbox(pnl_img, 'Position', [5, y_pos, 60, 20], 'Text', 'V', 'Value', 1);
    sld_green = uislider(pnl_img, 'Position', [65, y_pos+10, 120, 3], 'Limits', [0 5], 'Value', 1);
    
    y_pos = y_pos - esp;
    cb_blue = uicheckbox(pnl_img, 'Position', [5, y_pos, 60, 20], 'Text', 'B', 'Value', 1);
    sld_blue = uislider(pnl_img, 'Position', [65, y_pos+10, 120, 3], 'Limits', [0 5], 'Value', 1);
    
    y_pos = y_pos - esp;
    cb_ca = uicheckbox(pnl_img, 'Position', [5, y_pos, 60, 20], 'Text', 'Ca', 'Value', 0);
    sld_ca = uislider(pnl_img, 'Position', [65, y_pos+10, 120, 3], 'Limits', [0 5], 'Value', 1);

    cb_red.ValueChangedFcn = @(cb,event) update_image_channels(); sld_red.ValueChangedFcn = @(sld,event) update_image_channels();
    cb_green.ValueChangedFcn = @(cb,event) update_image_channels(); sld_green.ValueChangedFcn = @(sld,event) update_image_channels();
    cb_blue.ValueChangedFcn = @(cb,event) update_image_channels(); sld_blue.ValueChangedFcn = @(sld,event) update_image_channels();
    cb_ca.ValueChangedFcn = @(cb,event) update_image_channels(); sld_ca.ValueChangedFcn = @(sld,event) update_image_channels();

    % --- PANNEAU 1 : CLASSE (MILIEU) ---
    pnl_class = uipanel(fig, 'Position', [pnlX, figHeight-390, 200, 210], 'Title', 'Classe à appliquer');
    bg = uibuttongroup(pnl_class, 'Position', [5 5 190 180], 'BorderType', 'none');
    bg.SelectionChangedFcn = @(bg,event) change_selected_class(event.NewValue.Text);
    
    rb = gobjects(8,1);
    btn_y = 155; spacing = 22;
    for c = 1:8
        rb(c) = uiradiobutton(bg, 'Position', [10 btn_y 150 20], 'Text', noms_classes{c}, ...
                              'FontColor', couleurs_txt{c}, 'FontWeight', 'bold');
        if c == current_selected_class
            rb(c).Value = 1;
        end
        btn_y = btn_y - spacing;
    end

    % --- PANNEAU 2 : VISIBILITÉ (BAS) ---
    pnl_vis = uipanel(fig, 'Position', [pnlX, figHeight-610, 200, 210], 'Title', 'Afficher / Masquer');
    cb_vis = gobjects(8,1);
    btn_y = 165;
    for c = 1:8
        cb_vis(c) = uicheckbox(pnl_vis, 'Position', [10 btn_y 150 20], 'Text', noms_classes{c}, ...
                               'Value', 1, 'FontColor', couleurs_txt{c}, 'FontWeight', 'bold');
        cb_vis(c).ValueChangedFcn = @(cb,event) update_visibility();
        btn_y = btn_y - spacing;
    end

    % --- BOUTON VALIDER ---
    % On le place juste en dessous de la visibilité
    uibutton(fig, 'Position', [pnlX, figHeight-670, 200, 50], 'Text', 'VALIDER', ...
             'FontWeight', 'bold', 'FontSize', 16, 'ButtonPushedFcn', @(btn,event) valider_callback());
             
    fig.WindowButtonDownFcn = @(src, event) handle_image_click();
    
    app_is_running = true;
    while app_is_running && isgraphics(fig)
        uiwait(fig);
    end
    
    if isgraphics(fig)
        close(fig);
    end
    
    % --- Fonctions imbriquées ---
    function change_selected_class(txt)
        current_selected_class = str2double(txt(1));
    end

    function update_visibility()
        vis_states = zeros(8,1);
        for i = 1:8, vis_states(i) = cb_vis(i).Value; end
        
        for i = 1:NCell_total
            if isgraphics(h_rois(i))
                c_idx = colorcell_updated(i);
                if c_idx >= 1 && c_idx <= 8
                    if vis_states(c_idx) == 1
                        h_rois(i).Visible = 'on';
                        if isgraphics(h_markers(i)), h_markers(i).Visible = 'on'; end
                    else
                        h_rois(i).Visible = 'off';
                        if isgraphics(h_markers(i)), h_markers(i).Visible = 'off'; end
                    end
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
        r_layer(r_layer > 1) = 1; g_layer(g_layer > 1) = 1; b_layer(b_layer > 1) = 1;
        h_img.CData = cat(3, r_layer, g_layer, b_layer);
    end

    function handle_image_click()
        cp = ax.CurrentPoint;
        x_click = cp(1,1); y_click = cp(1,2);
        xlims = ax.XLim; ylims = ax.YLim;
        
        if x_click >= xlims(1) && x_click <= xlims(2) && y_click >= ylims(1) && y_click <= ylims(2)
            for i = 1:NCell_total
                if isgraphics(h_rois(i)) && strcmp(h_rois(i).Visible, 'on')
                    x_data = h_rois(i).XData; 
                    y_data = h_rois(i).YData;
                    
                    if inpolygon(x_click, y_click, x_data, y_data)
                        colorcell_updated(i) = current_selected_class;
                        
                        h_rois(i).FaceColor = cmap_classes(current_selected_class, :);
                        h_rois(i).EdgeColor = cmap_classes(current_selected_class, :);
                        
                        if isgraphics(h_markers(i))
                            delete(h_markers(i));
                        end
                        h_markers(i) = plot(ax, mean(x_data), mean(y_data), 'w*', 'MarkerSize', 8, ...
                                            'HitTest', 'off', 'PickableParts', 'none');
                        break; 
                    end
                end
            end
        end
    end

    function valider_callback()
        app_is_running = false; 
        uiresume(fig);
    end
end