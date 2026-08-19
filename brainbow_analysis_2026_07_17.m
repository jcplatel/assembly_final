clear;clc;close all
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

PathSave = "E:\Data\Aurelie\analysis\Final\assembly\nocues\";
% load("E:\Data\Aurelie\data\nocues\nocues_2026_06_30.mat");%PC
nocues = matfile("E:\Data\Aurelie\data\nocues_session.mat");
details = whos(nocues);
noms_toutes_variables = {details.name};
noms_sessions = noms_toutes_variables(startsWith(noms_toutes_variables, 'S_') | startsWith(noms_toutes_variables, 'Session_'));
nombre_de_manips = length(noms_sessions);

%% 2. PRÉPARATION DU JOURNAL (LOGS ET UI)
fig_log = uifigure('Name', 'Journal de Traitement', 'Position', [500, 300, 500, 400], 'Visible', 'on');
txt_log = uitextarea(fig_log, 'Position', [20, 20, 460, 360], 'Editable', 'off', 'FontName', 'Consolas'); 

log_file = fullfile(PathSave, sprintf('log_analyse_%s.txt', datestr(now, 'yy_mm_dd_HH_MM')));
historique = {};

% Message de démarrage
historique = ecrire_log(txt_log, log_file, historique, sprintf('=== DÉBUT DE L''ANALYSE : %s ===', datestr(now)));
temps_global = tic;

for file_num = 1:nombre_de_manips
    temps_fichier = tic; 
    % try
    close all

    name = noms_sessions{file_num}; 
    Session_Temp = nocues.(name); 
    identifier = Session_Temp.identifier;
    F = Session_Temp.F;
    sampling_rate = Session_Temp.sampling_rate;
    position = Session_Temp.position;
    path = Session_Temp.path;
    speed = Session_Temp.speed;
    clear Session_Temp; 

    daytime = datestr(now,'yy_mm_dd_HH_MM_SS');

    namefull = strcat(PathSave ,identifier ,'_',daytime  ,'\');%pc
    mkdir (namefull) ;

    % Log de progression
        msg = sprintf('Fichier %d/%d : [%s] en cours...', file_num, name);
        historique = ecrire_log(txt_log, log_file, historique, msg);
    %%


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
        'threshold_peak', 3.09, ...
        'synchronous_frames', 2, ...
        'sce_n_cells_threshold', 10, ...
        'SG_window', 7, ...
        'motion_correction', true, ...
        'colorsubstraction', false);

    %% Pre-processing: cell extraction, denoising, normalisation, baseline substraction, find SCE
        [Tr1b,speedsm,Raster,SumAct,MAct,Race,RasterRace,WinRest, WinActive,...
        TRace,Fzero,DFF0,th_detection,bad_frames,max_cells_allowed,n_transients_total,opts] = ...
        preprocessing_6(F,opts,speed,path,sampling_rate,namefull);
        
        [NCell, Nz] = size(Tr1b);
        find_ncluster = false;

        NClinit = 4:20;% 5:2:18;[8, 15]

    namefullold = namefull;

    for NClini=4:20%4:20%[12 ,15, 17, 18]%:5

        namefull = strcat (namefullold,'k',num2str(NClini),'\');
        mkdir (namefull) ;   % make folder for saving analysis

        % kmean_iter=50; kmeans_surrogate=10; kmeans_rnd_iter=10;savefig=0;
        kmean_iter = 200 ; kmeans_surrogate = 50 ; kmeans_rnd_iter = 20 ; savefig=1 ;
        SCE_clustering
    %     [IDX2, NClOK, assemblyortho, mean_sClOK, mean_Event_Purity, mean_Cellular_Fidelity, AmpTable,inertia] = ...
    % run_kmeans_consensus(Race, TRace, DFF0, Raster, NClini, kmean_iter, kmeans_surrogate, kmeans_rnd_iter, savefig, namefull);
        % results = biclustering_assemblies(Race, DFF0, TRace);
        % plot_bicluster_race(Race, results, namefull);
        % plot_bicluster_race_colored(Race, results, namefull);
        inertia_all(NClini) = mean(inertia);
        clear fig 
        exportdata

        clearvars cMaskCl cRace cCellP cRCl cC0 cNrnd

        % save(strcat(namefull,'results.mat')) 
        % save(strcat(namefull,'results.mat'), '-regexp', '^(?!(fig_log|ax|ax2|ax1|txt_log) $).');
        %raster_rastermap
        % path_colorcell="E:\Data\Aurelie\data\chroms\119\220923\registration\colorcell.mat"; 
        % path_colorcell="E:\Data\Aurelie\data\chroms\119\220919\registration\colorcell.mat"; 
        % path_colorcell="E:\Data\Aurelie\data\nocues\444119\220919_plane0\colorcell_Maxwell.mat";
        path_colorcell=strcat(path,"colorcell_Maxwell.mat");
        % if file_num == 21 || file_num == 26 || file_num == 49
        %     path_colorcell=strcat(path,"colorcell_GT.mat");
        % end% 
        
        if NCl>1
            save(strcat(namefull,'results.mat'), '-regexp', '^(?!(fig_log|ax|ax2|ax1|txt_log) $).');
            rastercolor
            graphSCE
            brainbowassemblies2026_05_07
            distance_calculation
            export_data_brainbow
            save(strcat(namefull,'brainbow.mat'),'-regexp',...
            '^(?!(fig_log|ax|ax2|ax1|txt_log|text_objs|fig|t|titre_obj|p|Table_All_Assemblies|ops|Ancien_Data|Fneu|spks) $).')
        end

    end
    dossier_mere = strcat (namefullold,'\');
    %find best 3 K
    analyse_assembly_results3
    % save them all in same file
    path_mat_file = (strcat(PathSave, 'bestK_all.mat'));

    if ~isfile(path_mat_file)
        % Première souris : création du fichier
        TableFinal_all = TableFinal;
        save(path_mat_file, 'TableFinal_all');
    else
        % Souris suivantes : chargement, concaténation, et écrasement
        Ancien_Data = load(path_mat_file, 'TableFinal_all');
        TableFinal_all = [Ancien_Data.TableFinal_all ; TableFinal];
        save(path_mat_file, 'TableFinal_all');
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