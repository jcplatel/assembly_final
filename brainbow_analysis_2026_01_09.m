clear
close all
% delete(gcp('nocreate'))
% parpool ('processes',4)
% profile on; % Activer le profiler

%% load path
if ispc
    % PathSave='E:\Data\Aurelie\analysis\December2025\assembly\nocuesPCA\';%PC
    PathSave='E:\Data\Aurelie\analysis\Jan2026\assembly\';
    % load("E:\Data\Aurelie\data\nocues\nocues.mat");%PC
    load("E:\Data\Aurelie\data\nwb_file\listnwb.mat");
    % PathSave='E:\Data\Aurelie\analysis\November2025\assembly\cuesPCA\';%PC
    % load("E:\Data\Aurelie\data\cues\cues.mat")
elseif ismac
    PathSave='/Users/platel/Desktop/exp/brainbow_hippocampus/analysis/nocues';%mac
    load("/Users/platel/Desktop/exp/brainbow_hippocampus/analysis/final_analysis/bestnocuesold.mat")
end
%%index = contains(matfile,'444152_221130');%227 
% 444119 19 = num 22
% for file_num=1:numel (matfile) %54  444178 at 57  65 crash à cause gros artefact 73 plein d'assemblées ani 198, 1306
% for file_num=[2:5,9,12:14,16,22:27 ,30:33,35,37,46,47,50,53:56,59,61:66,69,70]
    % for file_num=[5,9,12:14,16,22:27 ,30:33,35,37,46,47,50,53:56,59,61:66,69,70]
for file_num=171:184%119
    % try
    clearvars -except file_num matfile PathSave 
    close all
    filename=string(matfile{file_num});
    % filename="C:\Users\jcplatel\Documents\labo\manon\plane1\P22D_230301_230316_13_230323_plane0_2024_05_28.15-20-08.nwb";
    [path,name,ext] = fileparts(filename);
    identifier = name
    file_num
    path =strcat(path ,'\');%pc
    % path=strcat(path ,'/');%mac % filename =
    % '/Users/platel/Desktop/exp/brainbow_hippocampus/data/matfile/411582_230320_plane0.mat';
    % load (filename) % session=name; % str = [name ' / file_num= '
    openingnwb
    if exist ('position',"var"); [nb_lap , PosT] = lap_calculator (position,speed) ;
    else ; nb_lap='NaN';end
    %% options pre processing
    opts = struct();
    opts.MinPeakDistancesce = 5;
    opts.MinPeakDistance    = 3;
    opts.threshold_peak     = 2.576;
    opts.synchronous_frames = 2;
    opts.sce_n_cells_threshold = 10;
    opts.percentile         = NaN;
    opts.minithreshold      = 0.1;
    opts.SG_window          = 9;
    opts.use_PCA=false;
    opts.motion_correction=false;
    opts.colorsubstraction=false;
    %%
    [Tr1b,speedsm,Raster,SumAct,MAct,Race,Race_For_Clustering,RasterRace,WinRest, WinActive,TRace] = preprocessing_6(F,opts,speed,path);
    % preprocessing_4
    [NCell, Nz] = size(Tr1b);
 %cell extraction, denoising, normalisation, baseline substraction, find SCE
    % find_ncluster = true;
    find_ncluster = true;
    
    if find_ncluster == true %%find best K
        for ncluster = 4:20
            ncluster
            %clearvars -except start_PC n_bad_no_move nPC_Final Race_For_Clustering best_SClOK nb_lap Sil DB CH find_ncluster SumAct MAct best_S Nz minithreshold percentile synchronous_frames MinPeakDistance MinPeakDistancesce sce_n_cells_threshold Tr1b WinRest file_num matfile PathSave nanalysis NClini nanalysis path filename F iscell speed name identifier ncluster best_NCl sampling_rate session namefull Race
            % str=['test clustering ' num2str(ncluster) ' clusters'];
            daytime = datestr(now,'yy_mm_dd_HH_MM_SS');
            % fprintf ('%s ; ',str , daytime) 
            namefull = strcat(PathSave ,daytime ,'_',name ,'/');%pc
            % mkdir (namefull) ;   % make folder for saving analysis
            NClini  =ncluster;
            % NClini=5
            kmean_iter = 100;
            kmeans_surrogate = 50;
            kmeans_rnd_iter = 10;
            savefig = 0;
            % clustering_assembly_2
            clustering_PCA1
            % clustering_consensus_1
            % fprintf('silhouette: %f',  mean(sClOK))
            best_NCl(ncluster)=NCl;
            best_S(ncluster)=mean(sCl);
            best_SClOK(ncluster)=mean(sClOK);
        end
    end

    for nanalysis=1%4:20%[12 ,15, 17, 18]%:5

        %clearvars -except NCell n_bad_no_move start_PC nPC_Final Race_For_Clustering best_SClOK Sil CH DB find_ncluster MAct best_S Nz minithreshold percentile synchronous_frames MinPeakDistance MinPeakDistancesce sce_n_cells_threshold Tr1b WinRest file_num matfile PathSave nanalysis NClini path filename F iscell speed name identifier best_NCl sampling_rate session Race
        daytime = datestr(now,'yy_mm_dd_HH_MM_SS');
        % namefull=strcat(PathSave ,daytime ,'_',name ,'_exvivocolorPC2_30/');%pc
        namefull=strcat(PathSave ,daytime ,'_',name ,'exvivocluster_colorcellnew_invivo/');%pc
        mkdir (namefull) ;   % make folder for saving analysis

        if find_ncluster==true
            bestK2 
        else 
            NClini=7;
        end

        % NClini=nanalysis;
        fprintf ('clusters: %d ', NClini)
        kmean_iter=1000; kmeans_surrogate=100; kmeans_rnd_iter=100;savefig=1;
        % clustering_assembly_2 
        clustering_PCA1
        fprintf ('clusters: %d ', NCl)
        
        exportdata
        save(strcat(namefull,'results.mat')) 
        raster_rastermap

        if NCl>1
            rastercolor
            graphSCE
            %brainbowassemblies2025_11_26
            %distance_calculation
            % export_data_brainbow
            % save(strcat(namefull,'brainbow.mat') )
        end
    end

    % catch exception
    %     disp(exception.message);  % Display the error message
    % end
end