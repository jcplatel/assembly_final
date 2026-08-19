function [DFF0, Fzero] = baseline_calculation_glissante2(Tr1b, bad_frames, sampling_rate) 
    [NCell, Nz] = size(Tr1b);
    
    % Paramètres de la fenêtre glissante
    window_size = floor(sampling_rate * 120); 
    half_win = floor(window_size / 2);        
    step_size = floor(sampling_rate * 5);     
    percentile_value = 10;
    
    frames_pour_1sec = max(1, floor(sampling_rate * 1)); 
    
    centers = 1:step_size:Nz;
    if centers(end) ~= Nz
        centers = [centers, Nz]; 
    end
    num_steps = length(centers);
    
    DFF0 = zeros(NCell, Nz);
    Fzero = zeros(NCell, Nz);
    
    % ==============================================================
    % OPTIMISATION : Application globale du masque avant la boucle
    % ==============================================================
    Tr1b_masked = Tr1b; 
    if nargin >= 2 && ~isempty(bad_frames)
        % On remplace par NaN les colonnes (temps) correspondant aux bad_frames
        % pour TOUTES les lignes (cellules) d'un seul coup.
        Tr1b_masked(:, bad_frames) = NaN;
    end
    
    parfor n = 1:NCell
        trace_brute = Tr1b(n, :);       % Pour le calcul final de (F-F0)/F0
        trace_masked = Tr1b_masked(n, :); % Pour l'estimation de F0
        
        anchor_X = zeros(1, num_steps);
        anchor_Y = zeros(1, num_steps);
        
        for i = 1:num_steps
            c = centers(i);
            idx_s = max(1, c - half_win);
            idx_e = min(Nz, c + half_win);
            
            segment = trace_masked(idx_s:idx_e);
            segment_lisse = movmean(segment, frames_pour_1sec, 'omitnan');
            
            anchor_X(i) = c;
            anchor_Y(i) = prctile(segment_lisse, percentile_value);
        end
        
        valid_anchors = ~isnan(anchor_Y);
        aX = anchor_X(valid_anchors);
        aY = anchor_Y(valid_anchors);
        
        if length(aX) > 1
            F0 = interp1(aX, aY, 1:Nz, 'pchip', 'extrap');
        else
            F0 = repmat(nanmean(aY), 1, Nz);
        end
        
        % On utilise la trace brute pour le dF/F final
        DFF0(n, :) = (trace_brute - F0) ./ F0;
        Fzero(n, :) = F0;
    end
end