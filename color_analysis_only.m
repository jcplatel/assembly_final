clear;clc;close all

% load("E:\Data\Aurelie\data\nocues\nocues_2026_06_30.mat");%PC

dossier_mere = "E:\Data\Aurelie\analysis\Final\assembly\nocues\";
tous_les_fichiers = dir(fullfile(dossier_mere, '**', 'results.mat'));


for file_num = 1:length(tous_les_fichiers)
    % try
    path_results = fullfile(tous_les_fichiers(file_num).folder, tous_les_fichiers(file_num).name);
    parts = strsplit(tous_les_fichiers(file_num).folder, filesep);
    if length(parts) >= 2
        FolderName = parts{end};
        sessionName = parts{end-1};
    else
        continue;
    end
    load(path_results,)
    path_colorcell=strcat(path,"colorcell_Maxwell.mat");
    % if file_num == 21 || file_num == 26 || file_num == 49
    %     path_colorcell=strcat(path,"colorcell_GT.mat");
    % end% 
    
    if NCl>1
        brainbowassemblies2026_05_07
        distance_calculation
        export_data_brainbow
        save(strcat(namefull,'brainbow.mat'),'-regexp',...
        '^(?!(fig_log|ax|ax2|ax1|txt_log|text_objs|fig|t|titre_obj|p|Table_All_Assemblies|ops|Ancien_Data|Fneu|spks) $).')
    end
    
    % catch exception
    %     disp(exception.message);  % Display the error message
    %     nouvelle_ligne = sprintf('❌ ERREUR sur [%s] : %s', id_propre, exception.message);
    % end

end

