function figure_graphSCE(A, Tr1b, TRace, speed, namefull)
%% = graphSCE : vitesse + raster "HOT" (Tr1b trié) + lignes cyan aux SCE
%  où au moins une assemblée est significativement active.
%  Changement : la figure est fermée (elle s'accumulait d'un K à l'autre).
if A.NCl < 1, return; end
x1 = A.x1;
NCell = size(Tr1b, 1);
speed = double(speed(:)');

fig = figure('Name', 'Raster HOT et SCE', 'Color', 'k', 'Position', [100, 100, 1200, 800], 'Visible', 'off');

% --- 1. SCE "Ok" : au moins une assemblée active ---
is_SCE_Ok = sum(A.PCl_final, 1) > 0;
frames_SCE_Ok = TRace(is_SCE_Ok);

% --- 2. Haut : vitesse ---
ax1 = subplot(4, 1, 1);
set(ax1, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
hold on;
plot(speed, 'Color', [0.8 0.8 0.8], 'LineWidth', 1.5);
ylabel('Vitesse (cm/s)');
title('Vitesse', 'Color', 'w');
grid on; ax1.GridColor = 'w'; ax1.GridAlpha = 0.3;
if ~isempty(frames_SCE_Ok)
    valid_frames = frames_SCE_Ok(frames_SCE_Ok <= length(speed));
    plot(valid_frames, speed(valid_frames), '.', 'Color', 'c', 'MarkerSize', 10);
end
hold off;

% --- 3. Bas : raster HOT ---
ax2 = subplot(4, 1, [2 3 4]);
imagesc(Tr1b(x1, :));
colormap(ax2, hot);
caxis([0 1]);
hold on;
if ~isempty(frames_SCE_Ok)
    try
        xline(frames_SCE_Ok, 'Color', 'c', 'Alpha', 0.4, 'LineWidth', 1);
    catch
        yl = ylim;
        for t = frames_SCE_Ok
            plot([t t], yl, 'Color', 'c', 'LineWidth', 1);
        end
    end
end
set(ax2, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
ylabel('Sorted Cells', 'Color', 'w');
xlabel('Temps (Frames)', 'Color', 'w');
title(sprintf('Raster Plot - %d  SCE ', length(frames_SCE_Ok)), 'Color', 'w');
linkaxes([ax1, ax2], 'x');
xlim([1 length(speed)]);
ylim([0.5 NCell + 0.5]);
hold off;

exportgraphics(fig, strcat(namefull, 'graphSCE.png'), 'Resolution', 150, ...
    'BackgroundColor', 'black', 'ContentType', 'image');
close(fig)
end
