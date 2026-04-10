function [DFF0,Fzero] = baseline_calculation (Tr1b,bad_frames,sampling_rate) 

[NCell, Nz] = size(Tr1b);
window_size = floor(sampling_rate*120);%2 minutes was 60 sec
percentile_value = 10;
num_blocks = ceil(Nz / window_size);

for n = 1:NCell
    trace = Tr1b(n, :);
    Nz = length(trace);

    % 1. Masquer les bad frames
    trace_masked = trace;
    if exist('bad_frames', 'var') && ~isempty(bad_frames)
        trace_masked(bad_frames) = NaN;
    end

    % On va stocker les positions X (temps) et Y (valeur baseline)
    anchor_X = zeros(1, num_blocks);
    anchor_Y = zeros(1, num_blocks);

    for i = 1:num_blocks
        idx_s = (i-1) * window_size + 1;
        idx_e = min(i * window_size, Nz);

        segment = trace_masked(idx_s:idx_e);

        % Valeur du percentile pour ce bloc
        val = prctile(segment, percentile_value);

        % On définit le "centre" temporel de ce bloc
        anchor_X(i) = (idx_s + idx_e) / 2;
        anchor_Y(i) = val;
    end

    % Nettoyage des ancres (si un bloc entier était NaN)
    valid_anchors = ~isnan(anchor_Y);
    anchor_X = anchor_X(valid_anchors);
    anchor_Y = anchor_Y(valid_anchors);

    % Interpolation pour créer le vecteur F0 complet (1 x Nz)
    % 'pchip' ou 'spline' crée une courbe douce. 'linear' fait des lignes droites.
    % extrapolation ('extrap') est nécessaire pour les bords gauche/droite
    if length(anchor_X) > 1
        F0 = interp1(anchor_X, anchor_Y, 1:Nz, 'pchip', 'extrap');
        traceF0(n,:)= F0;
    else
        % Cas rare : signal trop court ou tout est NaN
        F0 = repmat(nanmean(anchor_Y), 1, Nz);
    end
  % % Calcul dF/F
    
    DFF0(n, :) = (trace - F0) ./ F0;
    Fzero(n,:)=F0;
    % figure(Visible="on");plot(trace);hold on; plot(F0);hold on
    % yyaxis right
    % plot(Fdetrend(n,:));
  % 

end