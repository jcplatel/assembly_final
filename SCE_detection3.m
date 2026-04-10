function [Race,SumAct,MAct,NRace,RasterRace,TRace, Raster,max_cells_allowed,opts] = SCE_detection3(Raster,opts,bad_frames) 
    % Note : J'ai ajouté 'Raster' en output, car il sera modifié.

    sce_n_cells_threshold = opts.sce_n_cells_threshold;
    MinPeakDistancesce = opts.MinPeakDistancesce;
    synchronous_frames = opts.synchronous_frames;

   
    % Paramètres optimisés pour 10 Hz
    % synchronous_frames = 1;     % Intégration de 2 frames (200 ms)
    max_jitter = 10;            % Décalage max de +/- 1 seconde (10 frames)
    NShfl = 100;                % Nombre de shuffles pour le seuil
    % 
    [NCell, Nz] = size(Raster);
    max_cells_SCE = floor(NCell*30/100);
    opts.max_cells_SCE = max_cells_SCE ;

    % Optionnel: Sécurité absolue, forcer les bad frames à 0 dans le Raster dès le début
    if exist('bad_frames','var') && ~isempty(bad_frames)
         Raster(:, bad_frames) = 0;
    end
    
    % =====================================================================
    % 1. NETTOYAGE PRÉVENTIF DES ÉVÉNEMENTS ABERRANTS (Artefacts optiques)
    % =====================================================================
    SumAct = sum(Raster, 1); % Activité instantanée (0 frame de décalage)
    
    max_cells_allowed = max(round(NCell * 0.20), max_cells_SCE); 
    toxic_frames = SumAct > max_cells_allowed;
    
    if any(toxic_frames)
        fprintf('Nettoyage Raster : %d frames aberrantes (> %d cellules) mises à zéro.\n', ...
                sum(toxic_frames), max_cells_allowed);
        % Dilatation (padding) pour nettoyer de +/- 1 frame autour de l'artefact
        toxic_frames_padded = conv(double(toxic_frames), [1 1 1], 'same') > 0;
        % On efface ces frames du Raster
        Raster(:, toxic_frames_padded) = 0;
        
        % On recalcule SumAct puisque Raster a été nettoyé
        SumAct = sum(Raster, 1); 
    end
    
    % =====================================================================
    % 2. CRÉATION DE MAct 
    % =====================================================================
    MAct = zeros(1, Nz - synchronous_frames);
    for i = 1:(Nz - synchronous_frames)
        MAct(i) = sum(max(Raster(:, i:i+synchronous_frames), [], 2));
    end
    
    % % =====================================================================
    % % 3. SHUFFLING (Local Jittering) ET CALCUL DU SEUIL STATISTIQUE
    % % =====================================================================
    % 
    % Sumactsh = zeros(Nz - synchronous_frames, NShfl);   
    % disp('Calcul du seuil statistique avec Local Jittering...');
    % 
    % for n = 1:NShfl
    %     Rastersh = zeros(NCell, Nz); 
    % 
    %     % Local Jittering indépendant (Zéro-Padding aux extrémités)
    %     for c = 1:NCell
    %         k = randi([-max_jitter, max_jitter]);
    %         if k > 0
    %             Rastersh(c, k+1:end) = Raster(c, 1:end-k);
    %         elseif k < 0
    %             Rastersh(c, 1:end+k) = Raster(c, 1-k:end);
    %         else
    %             Rastersh(c, :) = Raster(c, :);
    %         end
    %     end
    % 
    %     % Calcul des co-activations du Shuffle avec la MÊME méthode que MAct
    %     MActsh = zeros(1, Nz - synchronous_frames);
    %     for i = 1:(Nz - synchronous_frames)
    %         MActsh(i) = sum(max(Rastersh(:, i:i+synchronous_frames), [], 2));
    %     end
    % 
    %     Sumactsh(:, n) = MActsh;
    % end
    % 
    % % Détermination du seuil rigoureux sur les maximums pour contrer les artefacts
    % % max_per_shuffle = max(Sumactsh, [], 1); 
    % % sce_n_cells_threshold = ceil(prctile(max_per_shuffle, 95));
    % sce_n_cells_threshold = prctile(Sumactsh, 95, 'all');
    % fprintf('--> Seuil final retenu pour un SCE : %d cellules\n', sce_n_cells_threshold);
    
    % =====================================================================
    % 4. DÉTECTION FINALE DES SCE SUR MAct
    % =====================================================================
    % On cherche les pics sur MAct, car le seuil a été calculé avec une fenêtre similaire
    [pks, locs] = findpeaks(MAct, 'MinPeakHeight', sce_n_cells_threshold, 'MinPeakDistance', MinPeakDistancesce);
    
    TRace = locs;
    NRace = length(TRace);
    disp(['Ncells= ' num2str(NCell) ' | Sum transient: ' num2str(sum(SumAct)) ' | nSCE détectés : '  num2str(NRace)]);

    % =====================================================================
    % 5. CRÉATION DE LA MATRICE RACE (Extraction des cellules par SCE)
    % =====================================================================
    Race = zeros(NCell, NRace);       % Race = cellules participant au SCE 
    RasterRace = zeros(NCell, Nz);
    
    % Fenêtre symétrique autour du pic du SCE (pour 10Hz, +/- 1 frame = fenêtre de 300ms globale)
    window_radius = 1; 

    for i = 1:NRace
        % Protection contre les bords de la matrice
        idx_start = max(1, TRace(i) - 1);
        idx_end   = min(Nz, TRace(i) + 2);
        
        % Identification des cellules actives dans la fenêtre du SCE
        Race(:, i) = max(Raster(:, idx_start:idx_end), [], 2);    
        
        % Placement du pic exact dans RasterRace
        for c = 1:NCell
            if Race(c, i) == 1
                % Trouve la frame exacte du pic de cette cellule dans la fenêtre locale
                local_frames = idx_start:idx_end;
                active_frame = local_frames(Raster(c, local_frames) == 1);
                
                % Enregistre le premier pic trouvé dans la fenêtre
                if ~isempty(active_frame)
                    RasterRace(c, active_frame(1)) = 1; 
                end
            end
        end
    end
end