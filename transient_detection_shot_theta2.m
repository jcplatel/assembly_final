function [Raster, Acttmp2, th_detection, BurstMask, PromALL,WidthALL] = ...
    transient_detection_shot_theta2(Tr1b, MinPeakDistance, threshold_peak, WinRest, DFF0, namefull, Fzero)

[NCell, Nz] = size(Tr1b);
Raster = zeros(NCell, Nz);
Acttmp2 = cell(1, NCell);
ampALL  = cell(1, NCell);
PromALL = cell(1, NCell);
WidthALL = cell(1, NCell); 
th_detection = zeros(1, NCell);
BurstMask = false(NCell, Nz);
RejectedLocs   = cell(1, NCell);
RejectedReason = cell(1, NCell);

% =========================================================================
% PARAMÈTRES DE REJET DES GRANDS PICS / BURSTS
% =========================================================================
use_absolute_cutoff  = true;   % true = seuil fixe en dF/F0, false = seuil relatif au bruit
abs_prom_cutoff       = 0.5;   % <-- votre seuil dérivé du graphe (1.5 ou 2)
max_threshold_factor  = 5.0;   % utilisé seulement si use_absolute_cutoff = false
max_peak_width_frames = 10;    % Rejeter si le pic est trop large (> 1s à 10 Hz) test5 ou 10 frame...
% =========================================================================

total_peaks_rest      = 0;
total_peaks_kept      = 0;
rejected_amplitude    = 0;
rejected_width        = 0;
rejected_baseline_F0  = 0;

for i = 1:NCell     
    noise = 1.4826 * mad(diff(DFF0(i, WinRest)), 1) / sqrt(2);
    th_detection(i) = threshold_peak * noise;
    
    if use_absolute_cutoff
        th_max_detection = abs_prom_cutoff;          % seuil fixe, identique pour toutes les cellules
    else
        th_max_detection = max_threshold_factor * noise; % seuil relatif au bruit local
    end
    
    [amplitude, locs, w, p] = findpeaks(Tr1b(i, :), ...
        'MinPeakProminence', th_detection(i), ...
        'MinPeakDistance',   MinPeakDistance);
    
    peaks_Rest = ismember(locs, WinRest);
    locs_in_rest      = locs(peaks_Rest);
    amplitude_in_rest = amplitude(peaks_Rest);
    prom_in_rest      = p(peaks_Rest);
    width_in_rest     = w(peaks_Rest);
    
    n_rest_cell = length(locs_in_rest);
    total_peaks_rest = total_peaks_rest + n_rest_cell;
    
    if n_rest_cell == 0
        ampALL{i}  = [];
        PromALL{i} = [];
        Acttmp2{i} = [];
        RejectedLocs{i}   = [];  
        RejectedReason{i} = []; 
        continue;
    end
    
    valid_amplitude = prom_in_rest <= th_max_detection;
    valid_width     = width_in_rest <= max_peak_width_frames;
    
    F0_median = median(Fzero(i, :), 'omitnan');
    F0_iqr    = iqr(Fzero(i, :));
    ThBurst_F0 = F0_median + 2.5 * F0_iqr; 
    valid_F0 = Fzero(i, locs_in_rest) < ThBurst_F0;
    
    rejected_amplitude   = rejected_amplitude   + sum(~valid_amplitude);
    rejected_width       = rejected_width       + sum(~valid_width);
    % rejected_baseline_F0 = rejected_baseline_F0 + sum(~valid_F0);
    
    valid_peaks = valid_amplitude & valid_width ;%& valid_F0;
    
    reject_mask = ~valid_peaks;
    RejectedLocs{i} = locs_in_rest(reject_mask);
    reason = zeros(1, sum(reject_mask));
    idx_rej = find(reject_mask);
    for r = 1:length(idx_rej)
        bad_amp = ~valid_amplitude(idx_rej(r));
        bad_w   = ~valid_width(idx_rej(r));
        if bad_amp && bad_w
            reason(r) = 3;
        elseif bad_amp
            reason(r) = 1;
        elseif bad_w
            reason(r) = 2;
        end
    end
    RejectedReason{i} = reason;

    locs_final = locs_in_rest(valid_peaks);
    amp_final  = amplitude_in_rest(valid_peaks);
    prom_final = prom_in_rest(valid_peaks);
    width_final =  width_in_rest (valid_peaks);
    
    total_peaks_kept = total_peaks_kept + length(locs_final);
    
    ampALL{i}  = amp_final;
    PromALL{i} = prom_final;
    WidthALL{i} = width_final; 
    Acttmp2{i} = locs_final;
    Raster(i, locs_final) = 1;
    
    BurstMask(i, :) = Fzero(i, :) >= ThBurst_F0;
end

total_peaks_rejected = total_peaks_rest - total_peaks_kept;

fprintf('\n--- BILAN DU FILTRAGE DES PICS (WinRest) ---\n');
fprintf('  * Mode de coupure haute     : %s\n', ...
    ternary(use_absolute_cutoff, sprintf('Absolu (proéminence > %.2f dF/F0)', abs_prom_cutoff), sprintf('Relatif (> %.1fx bruit)', max_threshold_factor)));
fprintf('  * Total pics bruts détectés : %d\n', total_peaks_rest);
fprintf('  * Total pics conservés      : %d (%.1f%%)\n', total_peaks_kept, (total_peaks_kept/max(1,total_peaks_rest))*100);
fprintf('  * Total pics SUPPRIMÉS      : %d (%.1f%%)\n', total_peaks_rejected, (total_peaks_rejected/max(1,total_peaks_rest))*100);
fprintf('    - Rejets par amplitude/proéminence excessive : %d\n', rejected_amplitude);
fprintf('    - Rejets par largeur excessive (> %d frames)  : %d\n', max_peak_width_frames, rejected_width);
fprintf('    - Rejets par instabilité de ligne de base F0  : %d\n', rejected_baseline_F0);
fprintf('--------------------------------------------\n\n');

stats_filter = struct();
stats_filter.total_peaks_rest     = total_peaks_rest;
stats_filter.total_peaks_kept     = total_peaks_kept;
stats_filter.total_peaks_rejected = total_peaks_rejected;
stats_filter.rejected_amplitude   = rejected_amplitude;
stats_filter.rejected_width       = rejected_width;
stats_filter.rejected_baseline_F0 = rejected_baseline_F0;
%%
% =========================================================================
% PLOT QC : 10 traces aléatoires, CHAQUE TRACE DANS SA PROPRE FIGURE
% avec pics gardés (étoile verte) et pics rejetés car trop grands/trop
% larges (étoile rouge)
% =========================================================================
% Nécessite en workspace :
%   Tr1b            : matrice [NCell x Nz] des traces utilisées pour findpeaks
%   Acttmp2         : cell array des positions des pics CONSERVÉS (sortie fonction)
%   RejectedLocs    : cell array des positions des pics REJETÉS (sortie version QC)
%   RejectedReason  : (optionnel) 1=amplitude, 2=largeur, 3=les deux
%   sampling_rate   : fréquence d'acquisition (Hz), pour l'axe temps en secondes
% =========================================================================

[NCell, Nz] = size(Tr1b);
n_traces_to_plot = min(10, NCell);

% Sélection aléatoire de cellules (on privilégie celles avec au moins 1 pic rejeté)
has_rejected = cellfun(@(x) ~isempty(x), RejectedLocs);
pool_priority = find(has_rejected);
pool_fallback = setdiff(1:NCell, pool_priority);

if length(pool_priority) >= n_traces_to_plot
    selected_cells = pool_priority(randperm(length(pool_priority), n_traces_to_plot));
else
    n_missing = n_traces_to_plot - length(pool_priority);
    extra = pool_fallback(randperm(length(pool_fallback), min(n_missing, length(pool_fallback))));
    selected_cells = [pool_priority, extra];
end
selected_cells = sort(selected_cells);

t = (1:Nz);
if exist('sampling_rate', 'var') && ~isempty(sampling_rate)
    t = t / sampling_rate;
    xlab_str = 'Temps (s)';
else
    xlab_str = 'Frame';
end

% Dossier de sortie
if exist('namefull', 'var') && isfolder(namefull)
    out_dir = fullfile(namefull, 'QC_transients');
else
    out_dir = 'QC_transients';
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end

for k = 1:length(selected_cells)
    c = selected_cells(k);

    fig_qc = figure('Color', 'w', 'Position', [50, 50, 1600, 350], 'Visible', 'off');
    hold on;

    plot(t, Tr1b(c, :), 'k-', 'LineWidth', 1.0);

    % Pics conservés : étoile verte
    locs_kept = Acttmp2{c};
    if ~isempty(locs_kept)
        plot(t(locs_kept), Tr1b(c, locs_kept), 'g*', 'MarkerSize', 2, 'LineWidth', 1.8);
    end

    % Pics rejetés : étoile rouge
    locs_rej = RejectedLocs{c};
    if ~isempty(locs_rej)
        plot(t(locs_rej), Tr1b(c, locs_rej), 'r*', 'MarkerSize', 2, 'LineWidth', 1.8);
    end

    xlabel(xlab_str);
    ylabel('dF/F0 (Tr1b)');
    title(sprintf('Cellule %d — vert = conservé, rouge = rejeté (amplitude/largeur)', c));
    legend({'Trace', 'Conservé', 'Rejeté'}, 'Location', 'northeastoutside');
    box off;
    xlim([t(1), t(end)]);

    fname = fullfile(out_dir, sprintf('QC_transient_cell%03d.png', c));
    exportgraphics(fig_qc, fname, 'Resolution', 300, 'BackgroundColor', 'white');
    close(fig_qc);
end

%%
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end