function [Tr1b,speedsm,Raster,SumAct,MAct,Race,RasterRace,WinRest, WinActive,TRace,Fzero,DFF0,th_detection,bad_frames,max_cells_allowed,opts] = preprocessing_6(F,opts,speed,path,sampling_rate,namefull)

%% Load settings
MinPeakDistancesce = opts.MinPeakDistancesce;      % 5 default
MinPeakDistance = opts.MinPeakDistance;          % 3 default
threshold_peak = opts.threshold_peak;           % 3 default
synchronous_frames = opts.synchronous_frames;       % 2 default
sce_n_cells_threshold = opts.sce_n_cells_threshold; %10
SG_window = opts.SG_window ; %default 7
motion_correction=opts.motion_correction;%false
colorsubstraction=opts.colorsubstraction;%false

%% import Fluo
 Tr1b=double(F);

%% correct for motion artifact
if ~isfile(strcat(path, 'Fall.mat')); motion_correction = false ;end
if motion_correction==true 
    %path="E:\Data\Aurelie\data\nocues\444119\220919_plane0\";
    [Tr1b,bad_frames] = motion_correction_substraction (Tr1b,path,speed);
    save(strcat(path,'badframes.mat'),"bad_frames") 
end
%% correct for motion artifact
if colorsubstraction==true
    [Tr1b,colorcell] = color_substraction (Tr1b,path) ;
end

%% filtering, normalising

[NCell, Nz] = size(Tr1b);
speedsm = smoothdata(speed, 'gaussian', 50);

%% baseline calculation

if motion_correction==true
    [DFF0,Fzero] = baseline_calculation_glissante2 (Tr1b,bad_frames,sampling_rate);  
else
    bad_frames = false(1,Nz);
    [DFF0,Fzero] = baseline_calculation_glissante2 (Tr1b,bad_frames,sampling_rate);  
end

%% filtering
Tr1b = sgolayfilt(DFF0', 3, SG_window)';

%% définition fenetre active rest

% WinRest = find(speedsm <= 1) ;
% WinActive = find((speedsm > 2) & ~bad_frames');
isActive = (speedsm > 1); 
isActive_expanded = movmax(isActive, [3 3]);
if length(isActive_expanded)>Nz 
    isActive_expanded=isActive_expanded(1:Nz);
end
WinActive = find(isActive_expanded);
WinRest = find(~isActive_expanded);

%% Rest period has to be 3sec
WinRest_bin = false(1, Nz);
WinRest_bin(WinRest) = true;
taille_minimum = floor(sampling_rate*3); %3seconds
WinRest_bin = bwareaopen(WinRest_bin, taille_minimum);
WinRest = find(WinRest_bin);

%%  DÉTECTION DES TRANSIENTS

% Raster = transient_detection_mad (Tr1b,MinPeakDistancesce,MinPeakDistance,threshold_peak,WinRest,WinActive);
% [Raster,th] = transient_detection_std (Tr1b,MinPeakDistancesce,MinPeakDistance,threshold_peak,WinRest,WinActive);
[Raster,Acttmp2,th_detection] = transient_detection_shot_theta (Tr1b,MinPeakDistancesce,MinPeakDistance,threshold_peak,WinRest,WinActive,DFF0,namefull,Fzero);

%% 
% ÉLIMINER LES CELLULES HYPERACTIVES 
n_transients_total = sum(cellfun(@length, Acttmp2));
transients_percell = cellfun(@length, Acttmp2);
n_transients_percell = sum(cellfun(@length, Acttmp2))/NCell;
seuil_frequence = prctile(cellfun(@length, Acttmp2), 99);

seuil_frequence = 0.2 ;%1 pic tous les 5 sec
max_peak = max(length(WinRest)/sampling_rate * seuil_frequence,200) ;
cellules_hyperactives_idx = find(transients_percell > max_peak);

% if ~isempty(cellules_hyperactives_idx)
%     Raster(cellules_hyperactives_idx, :) = 0;
% end

% ===== DÉTECTION DES SCE (Synchronous Events) =====

[Race,SumAct,MAct,NRace,RasterRace,TRace, Raster,max_cells_allowed, opts] = SCE_detection3 (Raster,opts,bad_frames);

