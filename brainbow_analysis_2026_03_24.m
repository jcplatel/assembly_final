clear
close all
% delete(gcp('nocreate'))
% parpool ('processes',4)
% profile on; % Activer le profiler

% %%
% % 1. Définir les valeurs à tester pour chaque paramètre
% val_DistSce   = [3,5];            % 2 valeurs
% val_Dist      = [3,5];         % 3 valeurs
% val_Thresh    = [2.33,3.09];     % 2 valeurs
% val_Sync      = [0,1,2, 3];            % 2 valeurs
% val_SceNCells = [5, 10, 15];         % 3 valeurs
% val_SGWin     = [5,7,9];            % 2 valeurs
% 
% % val_DistSce   = [5];            % 2 valeurs
% % val_Dist      = [3];         % 3 valeurs
% % val_Thresh    = [2.33];     % 2 valeurs
% % val_Sync      = [2];            % 2 valeurs
% % val_SceNCells = [5];         % 3 valeurs
% % val_SGWin     = [7];            % 2 valeurs
% 
% [grid_DistSce, grid_Dist, grid_Thresh, grid_Sync, grid_SceNCells, grid_SGWin] = ...
%     ndgrid(val_DistSce, val_Dist, val_Thresh, val_Sync, val_SceNCells, val_SGWin);
% 
% test_DistSce   = grid_DistSce(:);
% test_Dist      = grid_Dist(:);
% test_Thresh    = grid_Thresh(:);
% test_Sync      = grid_Sync(:);
% test_SceNCells = grid_SceNCells(:);
% test_SGWin     = grid_SGWin(:);
% 
% total_tests = length(test_DistSce);

%% load path
if ispc
    PathSave="E:\Data\Aurelie\analysis\March2026\nocues\test_stability\";
    load("E:\Data\Aurelie\data\nocues\nocues.mat");%PC
    % load("E:\Data\Aurelie\data\nwb_file\listnwb.mat");
    % load("E:\Data\Aurelie\data\cues\cues.mat")
elseif ismac
    PathSave='/Users/platel/Desktop/exp/brainbow_hippocampus/analysis/nocues';%mac
    load("/Users/platel/Desktop/exp/brainbow_hippocampus/analysis/final_analysis/bestnocuesold.mat")
end

%% 2. PRÉPARATION DU JOURNAL (LOGS ET UI)
fig_log = uifigure('Name', 'Journal de Traitement', 'Position', [500, 300, 500, 400]);
txt_log = uitextarea(fig_log, 'Position', [20, 20, 460, 360], 'Editable', 'off', 'FontName', 'Consolas'); 

log_file = fullfile(PathSave, sprintf('log_analyse_%s.txt', datestr(now, 'yy_mm_dd_HH_MM')));
historique = {};

% Message de démarrage
historique = ecrire_log(txt_log, log_file, historique, sprintf('=== DÉBUT DE L''ANALYSE : %s ===', datestr(now)));
temps_global = tic;

files_to_process = [1,1,1,22,22, 22];%1:numel (matfile) 

for file_num = files_to_process
    temps_fichier = tic; 
    % try
    close all
    filename=string(matfile{file_num});
    [path,name,ext] = fileparts(filename);
    if strlength(name)> 23; identifier = extractBetween(name,5,24);else identifier='suite2p';end
    id_propre = strrep(identifier, '_', '-'); 
    daytime = datestr(now,'yy_mm_dd_HH_MM_SS');
    path =strcat(path ,'\');%pc
    namefull = strcat(PathSave ,identifier ,'_',daytime  ,'\');%pc
    mkdir (namefull) ;

    % Log de progression
        msg = sprintf('Fichier %d/%d : [%s] en cours...', file_num, numel(matfile), id_propre);
        historique = ecrire_log(txt_log, log_file, historique, msg);
    %%
    openingnwb; clearvars read_nwb
    % open_suite2p

    if exist ('position',"var"); [nb_lap , PosT] = lap_calculator (position,speed) ;
    else ; nb_lap='NaN';end

    %% options pre processing
    % opts.MinPeakDistancesce    = test_DistSce(test);
    % opts.MinPeakDistance       = test_Dist(test);
    % opts.threshold_peak        = test_Thresh(test);
    % opts.synchronous_frames    = test_Sync(test);
    % opts.sce_n_cells_threshold = test_SceNCells(test);
    % opts.SG_window             = test_SGWin(test);

   opts = struct(...
        'MinPeakDistancesce', 5, ...
        'MinPeakDistance', 3, ...
        'threshold_peak', 2.33, ...
        'synchronous_frames', 1, ...
        'sce_n_cells_threshold', 10, ...
        'SG_window', 7, ...
        'motion_correction', true, ...
        'colorsubstraction', false);

    %% Pre-processing: cell extraction, denoising, normalisation, baseline substraction, find SCE
 [Tr1b,speedsm,Raster,SumAct,MAct,Race,RasterRace,WinRest, WinActive,...
        TRace,Fzero,Fdetrend,th_detection,bad_frames,max_cells_allowed,opts] = ...
        preprocessing_6(F,opts,speed,path,sampling_rate,namefull);


    [NCell, Nz] = size(Tr1b);
    find_ncluster = false;

    %%find best K
    if find_ncluster == true 
        for ncluster = 4:20

            NClini  =ncluster;
            kmean_iter = 100;kmeans_surrogate = 50;kmeans_rnd_iter = 10;
            savefig = 0;

            SCE_clustering

            best_NCl(ncluster)=NCl;
            best_S(ncluster)=mean(sCl);
            best_SClOK(ncluster)=mean(sClOK);

        end
    end

    if find_ncluster==true
        bestK2 
        % best_NCl_interp=best_NCl;best_NCl_interp(best_NCl == 0) = NaN;
        % best_S_interp=best_S;best_S_interp(best_S == 0) = NaN;
        % best_SClOK_interp=best_SClOK;best_SClOK_interp(best_SClOK == 0) = NaN;
        % best_NCl_interp   = fillmissing(best_NCl_interp, 'linear');
        % best_S_interp     = fillmissing(best_S_interp, 'linear');
        % best_SClOK_interp = fillmissing(best_SClOK_interp, 'linear');
        % [~, idx_NCl]   = max(best_NCl_interp(4:end));
        % idx_NCl        = idx_NCl + 3; % Corrige l'indice puisqu'on a commencé à 4
        % [~, idx_SClOK] = max(best_SClOK_interp(4:end));
        % idx_SClOK      = idx_SClOK + 3;
        % [~, idx_S]     = max(best_S_interp(4:end));
        % idx_S          = idx_S + 3;

        % NClini = [NClini, idx_NCl, idx_SClOK, idx_S]; NClinit = unique(NClini, 'stable');
        [~, idx_NCl]   = max(best_NCl);
        [~, idx_SClOK] = max(best_SClOK);
        [~, idx_S]     = max(best_S);
        NClini = [NClini, idx_NCl, idx_SClOK, idx_S]; NClinit = unique(NClini, 'stable');
    else 
        NClinit = 4:18;% 5:2:18;[8, 15]
    end

    namefullold = namefull;

    for nanalysis=NClinit%4:20%[12 ,15, 17, 18]%:5

        NClini=nanalysis;

        namefull = strcat (namefullold,'/','k',num2str(nanalysis),'/');
        mkdir (namefull) ;   % make folder for saving analysis

        kmean_iter=1000; kmeans_surrogate=100; kmeans_rnd_iter=100;savefig=1;
        % kmean_iter = 100;kmeans_surrogate = 50;kmeans_rnd_iter = 10;savefig=1;
        SCE_clustering
        clear fig 
        exportdata

        clearvars cMaskCl cRace cCellP cRCl cC0 cNrnd

        % save(strcat(namefull,'results.mat')) 
        % save(strcat(namefull,'results.mat'), '-regexp', '^(?!(fig_log|ax|ax2|ax1|txt_log) $).');
        %raster_rastermap
        % path_colorcell="E:\Data\Aurelie\data\chroms\119\220923\registration\colorcell.mat"; 
        % path_colorcell="E:\Data\Aurelie\data\chroms\119\220919\registration\colorcell.mat"; 
        if NCl>1
            save(strcat(namefull,'results.mat'), '-regexp', '^(?!(fig_log|ax|ax2|ax1|txt_log) $).');
            rastercolor
            graphSCE
        %     % brainbowassemblies2025_11_26
        %     % distance_calculation
        %     % export_data_brainbow
        %     % save(strcat(namefull,'brainbow.mat') )
        end
    end

    % -- G. Fin du fichier & Logs --
    sec_fichier = toc(temps_fichier);
    sec_global = toc(temps_global);
    
    str_temps_f = string(duration(0, 0, sec_fichier, 'Format', 'mm:ss'));
    str_temps_g = string(duration(0, 0, sec_global, 'Format', 'hh:mm:ss'));
    
    msg = sprintf('✅ Fichier %d/%d [%s] terminé (Tps Fichier : %s | Total : %s)', ...
                  file_num, numel(matfile), id_propre, str_temps_f, str_temps_g);
    historique = ecrire_log(txt_log, log_file, historique, msg);

    % catch exception
    %     disp(exception.message);  % Display the error message
    %     nouvelle_ligne = sprintf('❌ ERREUR sur [%s] : %s', id_propre, exception.message);
    % end



end
% end

% =========================================================================
% FONCTION LOCALE : GESTION DES LOGS (UI + Fichier)
% =========================================================================
function historique = ecrire_log(ui_txt, path_file, historique, nouveau_message)
    % 1. Mise à jour de l'UI
    historique{end+1} = nouveau_message;
    ui_txt.Value = historique;
    scroll(ui_txt, 'bottom');
    drawnow;
    
    % 2. Écriture physique dans le fichier texte (.txt)
    fid = fopen(path_file, 'a', 'n', 'UTF-8');
    if fid ~= -1
        fprintf(fid, '%s\r\n', nouveau_message);
        fclose(fid);
    end
end