function [k_est, b_est] = estimer_fuite3(roi_col, roi_ca, k_limit)
% ESTIMER_FUITE3 (17 nov.) : k et b que bleedthrough_correction3 propose
% (fonction locale estimate_params, copiée à l'identique), sans la fenêtre.
% Méthode de l'élastique inférieur ; k borné par k_limit (limites de
% l'utilisateur : [vert rouge bleu] = [0.4 0.2 0.1]).
idx_v = (roi_ca > 0) & (roi_col > 0) & isfinite(roi_ca) & isfinite(roi_col);
x_v = roi_ca(idx_v); y_v = roi_col(idx_v);
if length(x_v) < 10
    k_est = 0; b_est = 0; return
end
x_aug = [x_v; max(x_v)]; y_aug = [y_v; max(y_v) * 10];
try
    k_h = convhull(x_aug, y_aug);
    x_h = x_aug(k_h); y_h = y_aug(k_h);
    [~, i_min_x] = min(x_h); i_max_x = max(x_v);
    x_l = []; y_l = [];
    for i = 1:length(x_h)
        if x_h(i) >= x_h(i_min_x) && x_h(i) <= i_max_x && y_h(i) <= max(y_v)
            x_l = [x_l; x_h(i)]; y_l = [y_l; y_h(i)]; %#ok<AGROW>
        end
    end
    if length(x_l) < 2
        x_l = [min(x_v); max(x_v)]; y_l = [min(y_v); min(y_v)];
    end
catch
    x_l = [min(x_v); max(x_v)]; y_l = [min(y_v); min(y_v)];
end
if length(x_l) >= 2
    pp = robustfit(x_l, y_l); b_est = pp(1); k_est = pp(2);
else
    k_est = 0; b_est = min(y_v);
end
if k_est < 0, k_est = 0; end
b_fit_orig = b_est; k_brid = false;
if k_est > k_limit, k_est = k_limit; k_brid = true; end
y_virt = y_v - k_est * x_v; b_est = prctile(y_virt, 5);
if k_brid
    b_max = min([b_fit_orig, max(y_v) * 0.05]);
    if b_est > b_max, b_est = b_max; end
end
if b_est < 0, b_est = 0; end
if b_est > min(y_v), b_est = min(y_v); end
end
