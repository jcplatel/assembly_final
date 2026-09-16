function [Raster, Acttmp2, th_detection, BurstMask, PromALL,WidthALL] = ...
    transient_detection_shot_theta(Tr1b, MinPeakDistance, threshold_peak, WinRest, DFF0, namefull, Fzero)

[NCell, Nz] = size(Tr1b);
Raster = zeros(NCell, Nz);
Acttmp2 = cell(1, NCell);
ampALL  = cell(1, NCell);
PromALL = cell(1, NCell);  
WidthALL = cell(1, NCell); 
th_detection = zeros(1, NCell);
BurstMask = false(NCell, Nz);
fenetre_correction = 400; 

for i = 1:NCell     
    noise = 1.4826 * mad(diff(DFF0(i, WinRest)), 1) / sqrt(2);
    th_detection(i) = threshold_peak * noise;
    
    [amplitude, locs, w, p] = findpeaks(Tr1b(i, :), ...
        'MinPeakProminence', th_detection(i), ...
        'MinPeakDistance', MinPeakDistance);
    
    peaks_Rest = ismember(locs, WinRest);
    locs_in_rest      = locs(peaks_Rest);
    amplitude_in_rest = amplitude(peaks_Rest);
    prom_in_rest      = p(peaks_Rest);
    width_in_rest     = w(peaks_Rest);
    
    F0_median = median(Fzero(i, :), 'omitnan');
    F0_iqr    = iqr(Fzero(i, :));
    ThBurst_F0 = F0_median + F0_iqr;
    BurstMask(i, :) = Fzero(i, :) >= ThBurst_F0;
    
    % Enregistrement final
    ampALL{i}  = amplitude_in_rest;  % hauteur brute du sommet (si besoin un jour)
    PromALL{i} = prom_in_rest;       % proéminence nette (à utiliser pour l'analyse d'amplitude)
    Acttmp2{i} = locs_in_rest;
    WidthALL{i} = width_in_rest; 
    Raster(i, locs_in_rest) = 1;
end

end

% for i = 1:NCell
%     if Acttmp2{i}>0
%         Raster(i,Acttmp2{i}) = 1;           %Raster = real raster of cell activity
%     end
% end

% cellules_choisies = randperm(NCell, 10);
% 
% for i=cellules_choisies
%     figure(Visible='off');
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
%                 exportgraphics(gcf,strcat(namefull, 'trace_cell_', num2str(i), '.png'),'Resolution',300, 'ContentType', 'image')
%     close gcf
% end