function [Raster,Acttmp2,th_detection,BurstMask] = transient_detection_shot_theta (Tr1b,MinPeakDistancesce,MinPeakDistance,threshold_peak,WinRest,WinActive,DFF0,namefull,Fzero)

[NCell, Nz] = size(Tr1b);
Raster = zeros(NCell,Nz);
Acttmp2 = cell(1,NCell);
ampALL  = cell(1,NCell);
th_detection = zeros(1, NCell);
fenetre_correction = 400; 

for i=1:NCell     
    noise= 1.4826 * mad(diff(DFF0(i, :)))/ sqrt(2);
    th_detection(i)= threshold_peak * noise ;%default=3.09
    % baseline_locale = movmedian(Tr1b(i, :), fenetre_correction);
    % Tr1b_plate = Tr1b(i, :) - baseline_locale;
    [amplitude, locs] = findpeaks(Tr1b(i, :),...
    'MinPeakHeight', th_detection(i), ...
    'MinPeakProminence', th_detection(i)*0.8, ...
    'MinPeakDistance', MinPeakDistance);

    peaks_Rest = ismember(locs, WinRest);
    locs_in_rest = locs(peaks_Rest);
    amplitude_in_rest = amplitude(peaks_Rest);

    F0_median = median(Fzero(i, :), 'omitnan');
    F0_iqr = iqr(Fzero(i, :));
    
    % % Seuil de Burst : par exemple, Médiane + 2 IQR (à ajuster selon tes données)
    % ThBurst_F0 = F0_median + F0_iqr; 
    % BurstMask(i, :) = Fzero(i, :) >= ThBurst_F0;
    % 
    % % On filtre les locs : on ne garde le pic que si le Fzero à cet instant
    % % est inférieur au seuil de Burst.
    % valid_non_burst = Fzero(i, locs_in_rest) < ThBurst_F0;
    % 
    % % Appliquer le filtre
    % locs_in_rest = locs_in_rest(valid_non_burst);
    % amplitude_in_rest = amplitude_in_rest(valid_non_burst);
    % % ========================================================
    
    % 5. Enregistrement final
    % ampALL{i} = amplitude_in_rest;
    % Acttmp2{i} = locs_in_rest;
    % Raster(i, locs_in_rest) = 1;
    ampALL{i} = amplitude_in_rest;
    Acttmp2{i} = locs_in_rest;
    Raster(i, locs_in_rest) = 1;
end

% for i = 1:NCell
%     if Acttmp2{i}>0
%         Raster(i,Acttmp2{i}) = 1;           %Raster = real raster of cell activity
%     end
% end

% cellules_choisies = randperm(NCell, 10);
% 
% for i=1
%     figure(Visible='on');
%     plot(Tr1b(i, :), 'b-', 'LineWidth', 1.2);
%     hold on;
%     pics_x = Acttmp2{i}; 
%     plot( pics_x, Tr1b(i, pics_x), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
%     yline(th_detection(i),'Color','red')
%     title (['cellule' , num2str(i)])
%     hold off
%     % yyaxis right
%     % % plot (speed)
%     % fenetre=zeros(1,Nz);
%     % fenetre(WinRest)=1;
%     % plot (fenetre)
%     %color coding rest period
%     y_lims = ylim;
%     y_min = y_lims(1);
%     y_max = y_lims(2);
%     is_active=ones(1,Nz);
%     is_active(WinRest)=0;
%     transitions = diff([0, is_active, 0]);
%     debuts = find(transitions == 1);
%     fins = find(transitions == -1) - 1;
%     X_rect = [debuts; debuts; fins; fins];
%     Y_rect = repmat([y_min; y_max; y_max; y_min], 1, length(debuts));
%     hold on;
%     p = patch(X_rect, Y_rect, [1, 0.8, 0.8], 'EdgeColor', 'none');
%     p.FaceAlpha = 0.4;
%     uistack(p, 'bottom');
% 
%     exportgraphics(gcf,strcat(namefull, 'trace_cell_', num2str(i), '.png'),'Resolution',300)
%     close gcf
% end