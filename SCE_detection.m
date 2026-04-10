function [Race,SumAct,MAct,NRace,RasterRace,TRace] = SCE_detection (Raster,opts,bad_frames) 


c
%% activité moyenne
sce_n_cells_threshold = opts.sce_n_cells_threshold;
MinPeakDistancesce = opts.MinPeakDistancesce;
synchronous_frames = opts.synchronous_frames;
[NCell, Nz] = size(Raster);

SumAct = sum(Raster, 1);
[pks, locs] = findpeaks(SumAct, 'MinPeakHeight', sce_n_cells_threshold, 'MinPeakDistance', MinPeakDistancesce);

% % Filtrer les SCE : enlever celles qui tombent sur des bad frames

if exist('bad_frames','var')
    sces_valid_mask = ~bad_frames(locs)';
    locs = locs(sces_valid_mask);
    pks = pks(sces_valid_mask);
end 

%% Filtrer les SCE aberrantes (outliers)
max_cells_allowed = max(round(NCell * 0.30),100); % Si NCell = 420, le max est 126
TF_too_big = pks > max_cells_allowed;

if any(TF_too_big)
    fprintf('Attention : %d SCE massifs (artefacts potentiels > %d cellules) ignorés.\n', sum(TF_too_big), max_cells_allowed);
    locs(TF_too_big) = [];
    pks(TF_too_big) = [];
end

% On met à jour TRace
TRace = locs;
NRace = length(TRace);


% absolute_threshold = 100;
% TF_relative = isoutlier(pks, "percentiles", [0 99]);
% TF_absolute = pks > absolute_threshold;
% TF_combined = TF_relative & TF_absolute;
% idx_to_remove = find(TF_combined);
% if ~isempty(idx_to_remove)
%     Raster(:, locs(idx_to_remove)) = 0;
% end

% Recalculer après filtrage
% SumAct = sum(Raster, 1);
% [pks, TRace] = findpeaks(SumAct, 'MinPeakHeight', sce_n_cells_threshold, 'MinPeakDistance', MinPeakDistancesce);

% NRace = length(TRace);
% fprintf('nSCE après filtrage bad_frames : %d\n', NRace);

% ===== CRÉER RACE MATRIX (Participation des cellules aux SCE) =====
MAct = zeros(1, Nz - synchronous_frames);
for i = 1:Nz - synchronous_frames
    MAct(i) = sum(max(Raster(:, i:i+synchronous_frames), [], 2));
end

%fprintf('Sum transient: %d\n', sum(MAct));

Race = false(NCell, NRace);
RasterRace = zeros(NCell, Nz);
for i = 1:NRace
    idx_start = max(1, TRace(i) - 1);
    idx_end = min(Nz, TRace(i) + 2);
    % Race(:, i) = max(Raster(:, TRace(i)-1:TRace(i)+2), [], 2);
    Race(:, i) = max(Raster(:, idx_start:idx_end), [], 2);
    RasterRace(Race(:, i) == 1, TRace(i)) = 1;
end

Race = double(Race);
%Race_For_Clustering = Race;