%% Traitement automatique par lot de tous les fichiers .nwb
clearvars; clc;

% 1. Sélection du dossier racine
rootFolder = "E:\Data\Aurelie\data\nocues\";

% 2. Recherche récursive de tous les fichiers .nwb dans le dossier et ses sous-dossiers
nwb_files = dir(fullfile(rootFolder, '**', '*.nwb'));

% Exclure les fichiers cachés/temporaires commençant par '._' (courants sur macOS / disques partagés)
valid_idx = ~startsWith({nwb_files.name}, '._');
nwb_files = nwb_files(valid_idx);

total_files = length(nwb_files);
if total_files == 0
    fprintf('Aucun fichier .nwb trouvé dans : %s\n', rootFolder);
    return;
end

fprintf('==> %d fichier(s) .nwb trouvé(s). Début du traitement...\n\n', total_files);

% 3. Boucle de traitement pour chaque enregistrement
success_count = 0;
failed_files = {};

for k = 1:total_files
    current_nwb_name = nwb_files(k).name;
    current_folder   = nwb_files(k).folder;
    selected_fullpath = fullfile(current_folder, current_nwb_name);
    
    fprintf('------------------------------------------------------------\n');
    fprintf('[%d/%d] Traitement de : %s\n', k, total_files, current_nwb_name);
    fprintf('Dossier : %s\n', current_folder);
    
    try
        % A. Variables d'environnement pour le script openingnwb
        dataset_path_folder = current_folder;
        dataset_path = [dataset_path_folder filesep];
        filename = selected_fullpath;
        folder = current_folder;
        
        % B. Recherche locale du fichier Fall.mat correspondant
        Fall_path_candidate = fullfile(current_folder, "Fall.mat");
        if exist(Fall_path_candidate, 'file')
            Fall_path = Fall_path_candidate;
        else
            % Recherche dans les sous-dossiers immédiats du dossier en cours (ex: suite2p/plane0)
            d_fall = dir(fullfile(current_folder, '**', 'Fall.mat'));
            d_fall = d_fall(~startsWith({d_fall.name}, '._'));
            
            if ~isempty(d_fall)
                Fall_path = fullfile(d_fall(1).folder, d_fall(1).name);
            else
                error('Fall.mat introuvable pour ce fichier NWB.');
            end
        end
        fprintf('  - Fall.mat trouvé : %s\n', Fall_path);
        
        % C. Lecture du fichier NWB via votre script existant
        openingnwb;
        
        % D. Chargement et filtrage de l'activité calcique
        fall_data = load(Fall_path, 'F', 'iscell');
        if ~isfield(fall_data, 'F') || ~isfield(fall_data, 'iscell')
            error('Le fichier Fall.mat ne contient pas les variables F ou iscell requises.');
        end
        
        % Filtrage des cellules valides uniquement
        F_filtered = double(fall_data.F(fall_data.iscell(:, 1) > 0, :));
        
        % E. Construction de la structure de sauvegarde
        data_recording = struct();
        data_recording.identifier    = identifier;
        data_recording.F             = F_filtered;
        data_recording.sampling_rate = sampling_rate;
        data_recording.speed         = speed;
        data_recording.position      = position;
        data_recording.path          = current_folder;
        
        % F. Enregistrement dans le dossier de l'enregistrement
        chemin_save_local = fullfile(current_folder, 'data_recording.mat');
        save(chemin_save_local, 'data_recording', '-v7.3');
        
        fprintf('  -> Enregistré avec succès : %s\n', chemin_save_local);
        success_count = success_count + 1;
        
    catch ME
        warning('Erreur lors du traitement de %s :\n%s', current_nwb_name, ME.message);
        failed_files{end+1} = selected_fullpath; %#ok<AGROW>
    end
end

% 4. Rapport de fin de traitement
fprintf('\n============================================================\n');
fprintf('TRAITEMENT TERMINÉ : %d/%d fichiers traités avec succès.\n', success_count, total_files);
if ~isempty(failed_files)
    fprintf('Fichiers en échec (%d) :\n', length(failed_files));
    for f = 1:length(failed_files)
        fprintf('  - %s\n', failed_files{f});
    end
end
fprintf('============================================================\n');