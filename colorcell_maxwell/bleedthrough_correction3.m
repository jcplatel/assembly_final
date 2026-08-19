function [roi_g_corr, roi_r_corr, roi_b_corr, img_g_corr, img_r_corr, img_b_corr] = ...
    bleedthrough_correction3(roi_g, roi_r, roi_b, roi_ca, limits_k, img_g, img_r, img_b, img_ca)
% BLEEDTHROUGH_CORRECTION3 Interface interactive pour corriger la fuite spectrale.
% 
% Arguments d'entrée :
%   - roi_g, roi_r, roi_b : Vecteurs des intensités des cellules pour chaque canal
%   - roi_ca              : Vecteur des intensités des cellules pour le canal GCaMP
%   - limits_k            : Vecteur [k_g_max, k_r_max, k_b_max] contenant les limites biologiques
%   - img_g, img_r, img_b : (Optionnel) Images entières des canaux couleurs
%   - img_ca              : (Optionnel) Image entière du canal GCaMP
%
% Arguments de sortie :
%   - roi_X_corr          : Vecteurs des intensités corrigées
%   - img_X_corr          : Images entières corrigées

    % Gestion des arguments optionnels (images)
    if nargin < 6
        img_g = []; img_r = []; img_b = []; img_ca = [];
    end

    k_green_limit = limits_k(1);
    k_red_limit   = limits_k(2);
    k_blue_limit  = limits_k(3);

    % =====================================================================
    % 1. ESTIMATION INITIALE (via fonction locale)
    % =====================================================================
    [k_g, b_g, x_g, y_g, xl_g, yl_g] = estimate_params(roi_g, roi_ca, k_green_limit);
    [k_r, b_r, x_r, y_r, xl_r, yl_r] = estimate_params(roi_r, roi_ca, k_red_limit);
    [k_b, b_b, x_b, y_b, xl_b, yl_b] = estimate_params(roi_b, roi_ca, k_blue_limit);

    % =====================================================================
    % 2. CRÉATION DE L'INTERFACE GRAPHIQUE INTERACTIVE
    % =====================================================================
    f_ui = uifigure('Name', 'Ajustement Interactif du Bleedthrough', 'Position', [100 100 1200 650]);
    grid_main = uigridlayout(f_ui, [3 3], 'RowHeight', {'1x', 50, 60});

    h_lines = struct(); % Structure pour stocker les handles des droites

    % --- Panneau VERT ---
    pnl_g = uipanel(grid_main); pnl_g.Layout.Row = 1; pnl_g.Layout.Column = 1;
    ax_g = uiaxes(uigridlayout(pnl_g, [1 1], 'Padding', [0 0 0 0]));
    scatter(ax_g, x_g, y_g, 30, [0.7 0.7 0.7], 'o'); hold(ax_g, 'on');
    scatter(ax_g, xl_g, yl_g, 50, 'k', 'filled');
    x_line_g = linspace(0, max(x_g), 100);
    h_lines.g = plot(ax_g, x_line_g, k_g * x_line_g + b_g, 'g-', 'LineWidth', 3);
    xlim(ax_g, [0 max(x_g)*1.05]); ylim(ax_g, [0 max(y_g)*1.05]); hold(ax_g, 'off');

    pnl_ctrl_g = uipanel(grid_main); pnl_ctrl_g.Layout.Row = 2; pnl_ctrl_g.Layout.Column = 1;
    uilabel(pnl_ctrl_g, 'Position', [10 15 20 20], 'Text', 'k:');
    ef_k_g = uieditfield(pnl_ctrl_g, 'numeric', 'Position', [30 15 60 22], 'Value', k_g);
    uilabel(pnl_ctrl_g, 'Position', [110 15 20 20], 'Text', 'b:');
    ef_b_g = uieditfield(pnl_ctrl_g, 'numeric', 'Position', [130 15 60 22], 'Value', b_g);

    % --- Panneau ROUGE ---
    pnl_r = uipanel(grid_main); pnl_r.Layout.Row = 1; pnl_r.Layout.Column = 2;
    ax_r = uiaxes(uigridlayout(pnl_r, [1 1], 'Padding', [0 0 0 0]));
    scatter(ax_r, x_r, y_r, 30, [0.7 0.7 0.7], 'o'); hold(ax_r, 'on');
    scatter(ax_r, xl_r, yl_r, 50, 'k', 'filled');
    x_line_r = linspace(0, max(x_r), 100);
    h_lines.r = plot(ax_r, x_line_r, k_r * x_line_r + b_r, 'g-', 'LineWidth', 3);
    xlim(ax_r, [0 max(x_r)*1.05]); ylim(ax_r, [0 max(y_r)*1.05]); hold(ax_r, 'off');

    pnl_ctrl_r = uipanel(grid_main); pnl_ctrl_r.Layout.Row = 2; pnl_ctrl_r.Layout.Column = 2;
    uilabel(pnl_ctrl_r, 'Position', [10 15 20 20], 'Text', 'k:');
    ef_k_r = uieditfield(pnl_ctrl_r, 'numeric', 'Position', [30 15 60 22], 'Value', k_r);
    uilabel(pnl_ctrl_r, 'Position', [110 15 20 20], 'Text', 'b:');
    ef_b_r = uieditfield(pnl_ctrl_r, 'numeric', 'Position', [130 15 60 22], 'Value', b_r);

    % --- Panneau BLEU ---
    pnl_b = uipanel(grid_main); pnl_b.Layout.Row = 1; pnl_b.Layout.Column = 3;
    ax_b = uiaxes(uigridlayout(pnl_b, [1 1], 'Padding', [0 0 0 0]));
    scatter(ax_b, x_b, y_b, 30, [0.7 0.7 0.7], 'o'); hold(ax_b, 'on');
    scatter(ax_b, xl_b, yl_b, 50, 'k', 'filled');
    x_line_b = linspace(0, max(x_b), 100);
    h_lines.b = plot(ax_b, x_line_b, k_b * x_line_b + b_b, 'g-', 'LineWidth', 3);
    xlim(ax_b, [0 max(x_b)*1.05]); ylim(ax_b, [0 max(y_b)*1.05]); hold(ax_b, 'off');

    pnl_ctrl_b = uipanel(grid_main); pnl_ctrl_b.Layout.Row = 2; pnl_ctrl_b.Layout.Column = 3;
    uilabel(pnl_ctrl_b, 'Position', [10 15 20 20], 'Text', 'k:');
    ef_k_b = uieditfield(pnl_ctrl_b, 'numeric', 'Position', [30 15 60 22], 'Value', k_b);
    uilabel(pnl_ctrl_b, 'Position', [110 15 20 20], 'Text', 'b:');
    ef_b_b = uieditfield(pnl_ctrl_b, 'numeric', 'Position', [130 15 60 22], 'Value', b_b);

    % --- Callbacks (mise à jour en direct) ---
    ef_k_g.ValueChangedFcn = @(src, event) update_line(h_lines.g, x_line_g, src.Value, ef_b_g.Value, ax_g, 'Vert');
    ef_b_g.ValueChangedFcn = @(src, event) update_line(h_lines.g, x_line_g, ef_k_g.Value, src.Value, ax_g, 'Vert');
    ef_k_r.ValueChangedFcn = @(src, event) update_line(h_lines.r, x_line_r, src.Value, ef_b_r.Value, ax_r, 'Rouge');
    ef_b_r.ValueChangedFcn = @(src, event) update_line(h_lines.r, x_line_r, ef_k_r.Value, src.Value, ax_r, 'Rouge');
    ef_k_b.ValueChangedFcn = @(src, event) update_line(h_lines.b, x_line_b, src.Value, ef_b_b.Value, ax_b, 'Bleu');
    ef_b_b.ValueChangedFcn = @(src, event) update_line(h_lines.b, x_line_b, ef_k_b.Value, src.Value, ax_b, 'Bleu');

    % Initialisation des titres
    update_line(h_lines.g, x_line_g, k_g, b_g, ax_g, 'Vert');
    update_line(h_lines.r, x_line_r, k_r, b_r, ax_r, 'Rouge');
    update_line(h_lines.b, x_line_b, k_b, b_b, ax_b, 'Bleu');

    % --- Bouton Valider ---
    pnl_boutons = uipanel(grid_main); pnl_boutons.Layout.Row = 3; pnl_boutons.Layout.Column = [1 3]; 
    uibutton(pnl_boutons, 'push', 'Text', '✅ Appliquer ces corrections finales', ...
        'Position', [450 10 300 40], 'FontSize', 14, ...
        'ButtonPushedFcn', @(src, event) uiresume(f_ui));

    % Attente de la validation utilisateur
    uiwait(f_ui);

    if ~isvalid(f_ui)
        error('Validation annulée par l''utilisateur.');
    end

    % Récupération des valeurs finales
    k_g_final = ef_k_g.Value; b_g_final = ef_b_g.Value;
    k_r_final = ef_k_r.Value; b_r_final = ef_b_r.Value;
    k_b_final = ef_k_b.Value; b_b_final = ef_b_b.Value;

    close(f_ui); 

    % =====================================================================
    % 3. CALCUL DES CORRECTIONS FINALES (ROIs & Images)
    % =====================================================================
    roi_g_corr = max(0, roi_g - (k_g_final * roi_ca + b_g_final));
    roi_r_corr = max(0, roi_r - (k_r_final * roi_ca + b_r_final));
    roi_b_corr = max(0, roi_b - (k_b_final * roi_ca + b_b_final));

    if ~isempty(img_g) && ~isempty(img_ca)
        img_g_corr = max(0, img_g - (k_g_final * img_ca + b_g_final));
        img_r_corr = max(0, img_r - (k_r_final * img_ca + b_r_final));
        img_b_corr = max(0, img_b - (k_b_final * img_ca + b_b_final));
    else
        img_g_corr = []; img_r_corr = []; img_b_corr = [];
    end
end

% =========================================================================
% FONCTIONS LOCALES (Invisibles de l'extérieur)
% =========================================================================

function update_line(h_line, x_data, k, b, ax, col_name)
    % Met à jour la droite et le titre en direct
    h_line.YData = k * x_data + b;
    legend(ax, 'Cellules', 'Points Convex Hull', sprintf('Fit: y = %.3fx + %.3f', k, b), 'Location', 'best');
    title(ax, sprintf('%s (k=%.3f, b=%.3f)', col_name, k, b));
end

function [k_est, b_est, x_v, y_v, x_l, y_l] = estimate_params(roi_col, roi_ca, k_limit)
    % Calcule les paramètres k et b initiaux (méthode de l'élastique inférieur)
    idx_v = (roi_ca > 0) & (roi_col > 0) & isfinite(roi_ca) & isfinite(roi_col);
    x_v = roi_ca(idx_v); y_v = roi_col(idx_v);
    
    if length(x_v) < 10
        k_est = 0; b_est = 0; x_l = []; y_l = []; return;
    end
    
    x_aug = [x_v; max(x_v)]; y_aug = [y_v; max(y_v)*10]; 
    try
        k_h = convhull(x_aug, y_aug);
        x_h = x_aug(k_h); y_h = y_aug(k_h);
        [~, i_min_x] = min(x_h); i_max_x = max(x_v); 
        x_l = []; y_l = [];
        for i = 1:length(x_h)
            if x_h(i) >= x_h(i_min_x) && x_h(i) <= i_max_x && y_h(i) <= max(y_v)
                x_l = [x_l; x_h(i)]; y_l = [y_l; y_h(i)];
            end
        end
        if length(x_l) < 2
            x_l = [min(x_v); max(x_v)]; y_l = [min(y_v); min(y_v)];
        end
    catch
        x_l = [min(x_v); max(x_v)]; y_l = [min(y_v); min(y_v)];
    end
    
    if length(x_l) >= 2
        p = robustfit(x_l, y_l); b_est = p(1); k_est = p(2);
    else
        k_est = 0; b_est = min(y_v);
    end
    
    if k_est < 0, k_est = 0; end
    b_fit_orig = b_est; k_brid = false;
    
    if k_est > k_limit, k_est = k_limit; k_brid = true; end
    y_virt = y_v - k_est * x_v; b_est = prctile(y_virt, 5); 
    
    if k_brid
        b_max = min([b_fit_orig, max(y_v)*0.05]); 
        if b_est > b_max, b_est = b_max; end
    end
    
    if b_est < 0, b_est = 0; end
    if b_est > min(y_v), b_est = min(y_v); end
end