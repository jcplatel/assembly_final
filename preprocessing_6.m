function [Tr1b,speedsm,Raster,SumAct,MAct,Race,Race_For_Clustering,RasterRace,WinRest, WinActive,TRace] = preprocessing_6(F,opts,speed,path)

%% Load settings
MinPeakDistancesce = opts.MinPeakDistancesce;      % 5 default
MinPeakDistance = opts.MinPeakDistance;          % 3 default
threshold_peak = opts.threshold_peak;           % 3 default
synchronous_frames = opts.synchronous_frames;       % 2 default
sce_n_cells_threshold = opts.sce_n_cells_threshold; %10
percentile = opts.percentile;
minithreshold=opts.minithreshold;
SG_window = opts.SG_window ; %default 7
use_PCA=opts.use_PCA; %false
motion_correction=opts.motion_correction;%false
colorsubstraction=opts.colorsubstraction;%false

%% import Fluo
 Tr1b=double(F);

%% correct for motion artifact
if motion_correction==true
    [Tr1b,bad_frames_no_movement] = motion_correction_substraction (Tr1b,path,speed); 
end
%% correct for motion artifact
if colorsubstraction==true
    [Tr1b,colorcell] = color_substraction (Tr1b,path) ;
end
%% filtering, normalising
Tr1b = sgolayfilt(Tr1b', 3, SG_window)';
[NCell, Nz] = size(Tr1b);

speedsm = smoothdata(speed, 'gaussian', 50);

%% baseline calculation

if motion_correction==true
    Tr1b = baseline_calculation (Tr1b,bad_frames);  
else
    Tr1b = baseline_calculation (Tr1b);  
end

%% définition fenetre active rest

WinRest = find(speedsm <= 1) ;
% WinActive = find((speedsm > 2) & ~bad_frames');
WinActive = find(speedsm > 2);

%%  DÉTECTION DES TRANSIENTS

% Raster = transient_detection_mad (Tr1b,MinPeakDistancesce,MinPeakDistance,threshold_peak,WinRest,WinActive);
Raster = transient_detection_std (Tr1b,MinPeakDistancesce,MinPeakDistance,threshold_peak,WinRest,WinActive);


%% 
% ÉLIMINER LES CELLULES HYPERACTIVES 

% nombre_transients_par_cellule = cellfun(@length, Acttmp2);
% seuil_frequence = prctile(nombre_transients_par_cellule, 99);
% cellules_hyperactives_idx = find(nombre_transients_par_cellule > seuil_frequence);
% if ~isempty(cellules_hyperactives_idx)
%     Raster(cellules_hyperactives_idx, :) = 0;
% end

% ===== DÉTECTION DES SCE (Synchronous Events) =====
%%%%shuffling to find threshold for number of cell for sce detection
% MActsh = zeros(1,Nz-synchronous_frames);   
% Rastersh=zeros(NCell,Nz);   
% NShfl=100;
% Sumactsh=zeros(Nz-synchronous_frames,NShfl);   
% for n=1:NShfl
% 
%         for c=1:NCell
%             k = randi(Nz-length(WinActive));
%             Rastersh(c,:)= circshift(Raster(c,:),k,2);
%         end
% 
%         for i=1:Nz-synchronous_frames   %need to use WinRest???
%             MActsh(i) = sum(max(Rastersh(:,i:i+synchronous_frames),[],2));
%         end
% 
%     Sumactsh(:,n)=MActsh;
% end
% 
% percentile = 95; % Calculate the 5% highest point or 99
% sce_n_cells_threshold = prctile(Sumactsh, percentile,"all");
%% activité moyenne

SumAct = sum(Raster, 1);
[pks, locs] = findpeaks(SumAct, 'MinPeakHeight', sce_n_cells_threshold, 'MinPeakDistance', MinPeakDistancesce);

% % Filtrer les SCE : enlever celles qui tombent sur des bad frames
% if exist('bad_frames','var')
%     sces_valid_mask = ~bad_frames(locs)';
%     locs = locs(sces_valid_mask);
%     pks = pks(sces_valid_mask);
% end 

% Filtrer les SCE aberrantes (outliers)
% absolute_threshold = 100;
% TF_relative = isoutlier(pks, "percentiles", [0 99]);
% TF_absolute = pks > absolute_threshold;
% TF_combined = TF_relative & TF_absolute;
% idx_to_remove = find(TF_combined);
% if ~isempty(idx_to_remove)
%     Raster(:, locs(idx_to_remove)) = 0;
% end

% Recalculer après filtrage
SumAct = sum(Raster, 1);
[pks, TRace] = findpeaks(SumAct, 'MinPeakHeight', sce_n_cells_threshold, 'MinPeakDistance', MinPeakDistancesce);

NRace = length(TRace);
fprintf('nSCE après filtrage bad_frames : %d\n', NRace);

% ===== CRÉER RACE MATRIX (Participation des cellules aux SCE) =====
MAct = zeros(1, Nz - synchronous_frames);
for i = 1:Nz - synchronous_frames
    MAct(i) = sum(max(Raster(:, i:i+synchronous_frames), [], 2));
end

fprintf('Sum transient: %d\n', sum(MAct));

Race = false(NCell, NRace);
RasterRace = zeros(NCell, Nz);
for i = 1:NRace
    Race(:, i) = max(Raster(:, TRace(i)-1:TRace(i)+2), [], 2);
    RasterRace(Race(:, i) == 1, TRace(i)) = 1;
end

Race = double(Race);
Race_For_Clustering = Race;
