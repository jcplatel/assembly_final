clear;clc;close all
% delete(gcp('nocreate'))
% parpool ('processes',4)
% profile on; % Activer le profiler

%% 1. Définir les paramètres (Grid Search)
% (Laissés commentés ou définis fixes comme dans votre script)
% ... [Votre code pour ndgrid reste ici si vous l'utilisez] ...

%% 2. CHEMINS ET CHARGEMENT DE LA BASE DE DONNÉES
PathSave = "E:\Data\Aurelie\analysis\Final\assembly\cues\";

% --- NOUVEAU : Chargement depuis la base de données globale ---
path_db = "E:\Data\Aurelie\data\Global_Database.mat";
database = matfile(path_db);

% On récupère la liste de toutes les variables
details = whos(database);
noms_toutes_variables = {details.name};

% On ne garde que les noms de variables qui correspondent à des sessions (S_...)
noms_sessions = noms_toutes_variables(startsWith(noms_toutes_variables, 'S_'));
nombre_de_manips_total = length(noms_sessions);

% Définir la condition que l'on veut analyser avec ce script
condition_a_analyser = 'nocues'; 

%% 3. PRÉPARATION DU JOURNAL (LOGS ET UI)
fig_log = uifigure('Name', 'Journal de Traitement', 'Position', [500, 300, 500, 400], 'Visible', 'on');
txt_log = uitextarea(fig_log, 'Position', [20, 20, 460, 360], 'Editable', 'off', 'FontName', 'Consolas'); 
log_file = fullfile(PathSave, sprintf('log_analyse_%s.txt', datestr(now, 'yy_mm_dd_HH_MM')));
historique = {};

% Message de démarrage
historique = ecrire_log(txt_log, log_file, historique, sprintf('=== DÉBUT DE L''ANALYSE : %s ===', datestr(now)));
temps_global = tic;

sessions_traitees = 0; % Compteur des sessions réellement analysées

%% 4. BOUCLE PRINCIPALE SUR LES FICHIERS
for file_num = 1:nombre_de_manips_total
    temps_fichier = tic; 
    close all
    
    % --- CHARGEMENT DYNAMIQUE ---
    name = noms_sessions{file_num}; 
    Session_Temp = database.(name); % Va chercher uniquement cette session sur le disque
    
    % --- FILTRE CONDITION ---
    % Si la session n'est pas un 'nocues' (ou si la condition n'est pas définie), on l'ignore
    if ~isfield(Session_Temp, 'Condition') || ~strcmp(Session_Temp.Condition, condition_a_analyser)
        clear Session_Temp; % On nettoie la RAM
        continue;           % On passe direct au fichier suivant
    end
    
    sessions_traitees = sessions_traitees + 1;
    
    % Extraction des variables nécessaires depuis la structure
    identifier = Session_Temp.identifier;
    F = Session_Temp.F;
    sampling_rate = Session_Temp.sampling_rate;
    
    if isfield(Session_Temp, 'position')
        position = Session_Temp.position;
    else
        position = [];
    end
    
    if isfield(Session_Temp, 'path')
        path = char(Session_Temp.path); % Assure que ce soit une chaine de caractères
    else
        path = '';
    end
    
    if isfield(Session_Temp, 'speed')
        speed = Session_Temp.speed;
    else
        speed = [];
    end
    
    % --- NETTOYAGE RAM TRÈS IMPORTANT ---
    clear Session_Temp; 
    
    % --- PRÉPARATION SAUVEGARDE ---
    daytime = datestr(now,'yy_mm_dd_HH_MM_SS');
    namefull = strcat(PathSave, identifier, '_', daytime, '\'); % pc
    mkdir(namefull);
    
    % Log de progression
    msg = sprintf('Fichier %d/%d : [%s] en cours...', sessions_traitees, nombre_de_manips_total, identifier);
    historique = ecrire_log(txt_log, log_file, historique, msg);
    
    %% CALCULS (Lap, Preprocessing, SCE...)
    if ~isempty(position)
        [nb_lap, PosT] = lap_calculator(position, speed);
    else
        nb_lap = NaN;
    end
    
    % --- Options de pré-processing ---
    opts = struct(...
        'MinPeakDistancesce', 5, ...
        'MinPeakDistance', 3, ...
        'threshold_peak', 3.09, ...
        'synchronous_frames', 2, ...
        'sce_n_cells_threshold', 10, ...
        'SG_window', 7, ...
        'motion_correction', true, ...
        'colorsubstraction', false);
        
    % --- Pre-processing ---
    [Tr1b, speedsm, Raster, SumAct, MAct, Race, RasterRace, WinRest, WinActive,...
     TRace, Fzero, DFF0, th_detection, bad_frames, max_cells_allowed, n_transients_total, opts] = ...
     preprocessing_6(F, opts, speed, path, sampling_rate, namefull);
        
    [NCell, Nz] = size(Tr1b);
    find_ncluster = false;
    NClinit = 4:20; 
    
    namefullold = namefull;
    inertia_all = zeros(1, 20); % Pré-allocation
    
    for NClini = 4:20
        namefull = strcat(namefullold, 'k', num2str(NClini), '\');
        mkdir(namefull); 
        
        kmean_iter = 200; kmeans_surrogate = 50; kmeans_rnd_iter = 20; savefig = 1;
        
        % Lancement du script de clustering (assurez-vous qu'il utilise les bonnes variables)
        SCE_clustering
        
        % Si inertia est défini dans SCE_clustering
        if exist('inertia', 'var')
            inertia_all(NClini) = mean(inertia);
        end
        
        clear fig 
        exportdata
        clearvars cMaskCl cRace cCellP cRCl cC0 cNrnd
        
        path_colorcell = strcat(path, "colorcell_Maxwell.mat");
        
        if NCl > 1
            save(strcat(namefull, 'results.mat'), '-regexp', '^(?!(fig_log|ax|ax2|ax1|txt_log) $).');
            rastercolor
            graphSCE
            brainbowassemblies2026_05_07
            distance_calculation
            export_data_brainbow
            
            save(strcat(namefull, 'brainbow.mat'), '-regexp',...
                '^(?!(fig_log|ax|ax2|ax1|txt_log|text_objs|fig|t|titre_obj|p|Table_All_Assemblies|ops|Ancien_Data|Fneu|spks) $).')
        end
    end % Fin de la boucle K
    
    dossier_mere = strcat(namefullold, '\');
    
    % Find best 3 K
    analyse_assembly_results3
    
    % Save them all in same file
    path_mat_file = strcat(PathSave, 'bestK_all.mat');
    if ~isfile(path_mat_file)
        TableFinal_all = TableFinal;
        save(path_mat_file, 'TableFinal_all');
    else
        Ancien_Data = load(path_mat_file, 'TableFinal_all');
        TableFinal_all = [Ancien_Data.TableFinal_all; TableFinal];
        save(path_mat_file, 'TableFinal_all');
    end
    
    % -- Fin du fichier & Logs --
    sec_fichier = toc(temps_fichier);
    sec_global = toc(temps_global);
    str_temps_f = string(duration(0, 0, sec_fichier, 'Format', 'mm:ss'));
    str_temps_g = string(duration(0, 0, sec_global, 'Format', 'hh:mm:ss'));
    
    msg = sprintf('✅ [%s] terminé (Tps Fichier : %s | Total : %s)', ...
                  identifier, str_temps_f, str_temps_g);
    historique = ecrire_log(txt_log, log_file, historique, msg);
    
end % Fin de la boucle des sessions

historique = ecrire_log(txt_log, log_file, historique, sprintf('=== FIN DE L''ANALYSE === (%d sessions traitées)', sessions_traitees));

% =========================================================================
% FONCTION LOCALE : GESTION DES LOGS (UI + Fichier)
% =========================================================================
function historique = ecrire_log(ui_txt, path_file, historique, nouveau_message)
    historique{end+1} = nouveau_message;
    ui_txt.Value = historique;
    scroll(ui_txt, 'bottom');
    drawnow;
    
    fid = fopen(path_file, 'a', 'n', 'UTF-8');
    if fid ~= -1
        fprintf(fid, '%s\r\n', nouveau_message);
        fclose(fid);
    end
end