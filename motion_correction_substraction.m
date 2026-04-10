function [Tr1b,bad_frames] = motion_correction_substraction (Tr1b,path,speed) 

load(strcat(path, 'Fall.mat'), 'ops');

corrXY = ops.corrXY;

% Approche robuste : Déviation par rapport à la tendance locale
rolling_median = movmedian(corrXY, 300); 
deviation = corrXY - rolling_median;

% Bad frames = celles qui dévient fortement vers le bas
sigma_dev = std(deviation(deviation < 0));
seuil_bad = -3 * sigma_dev;
bad_frames = deviation < seuil_bad;
bad_frames = conv(double(bad_frames), [1 1 1], 'same') > 0;% to be safe we remove previous and next frames

% fprintf('Bad frames détectées : %d (%.2f%%)\n', ...
%     sum(bad_frames), 100*sum(bad_frames)/length(corrXY));

% Bad frames sans mouvement
bad_frames_no_movement = bad_frames & (speed' < 1);
n_bad_no_move = sum(bad_frames_no_movement);
% fprintf('Bad frames avec speed < 2 cm/s : %d (%.2f%%)\n', ...
%     n_bad_no_move, 100*n_bad_no_move/length(corrXY));
Tr1b_clean = Tr1b;
Tr1b_clean(:, bad_frames) = NaN;
Tr1b_clean = fillmissing(Tr1b_clean, 'linear', 2, 'EndValues', 'nearest'); %interpolation
Tr1b = Tr1b_clean;