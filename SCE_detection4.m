function [Race,SumAct,MAct,NRace,RasterRace,TRace,Raster,max_cells_allowed,opts] = SCE_detection4(Raster,opts,bad_frames)

    sce_n_cells_threshold = opts.sce_n_cells_threshold;
    MinPeakDistancesce    = opts.MinPeakDistancesce;
    synchronous_frames    = opts.synchronous_frames;
   
    max_jitter = 10;            % Décalage max (+/- 10 frames)
    NShfl      = 100;           % Nombre de shuffles pour le seuil
    
    [NCell, Nz] = size(Raster);
    
    % Plafond maximal absolu de cellules pour une assemblée
    max_cells_SCE     = 150;
    max_cells_allowed = max_cells_SCE;
    opts.max_cells_SCE = max_cells_SCE;

    % 0. Suppression des bad frames connues
    if exist('bad_frames','var') && ~isempty(bad_frames)
         Raster(:, bad_frames) = 0;
    end
    
    % =====================================================================
    % 1. NETTOYAGE PRÉVENTIF DES ARTEFACTS INSTANTANÉS
    % =====================================================================
    SumAct = sum(Raster, 1);
    toxic_frames = SumAct > max_cells_allowed;
    
    if any(toxic_frames)
        fprintf('Nettoyage Raster : %d frames instantanées (> %d cellules) mises à zéro.\n', ...
                sum(toxic_frames), max_cells_allowed);
        toxic_frames_padded = conv(double(toxic_frames), [1 1 1], 'same') > 0;
        Raster(:, toxic_frames_padded) = 0;
        SumAct = sum(Raster, 1); 
    end
    
    % =====================================================================
    % 2. CALCUL DU VECTEUR D'ACTIVITÉ CO-SYNCHRONE (MAct)
    % =====================================================================
    MAct = zeros(1, Nz - synchronous_frames);
    for i = 1:(Nz - synchronous_frames)
        MAct(i) = sum(max(Raster(:, i:i+synchronous_frames), [], 2));
    end
    
    % =====================================================================
    % 3. SHUFFLING (Local Jittering) ET SEUIL STATISTIQUE
    % =====================================================================
    Sumactsh = zeros(Nz - synchronous_frames, NShfl);   
    disp('Calcul du seuil statistique avec Local Jittering...');
    
    parfor n = 1:NShfl
        Rastersh = zeros(NCell, Nz); 
        for c = 1:NCell
            k = randi([-max_jitter, max_jitter]);
            if k > 0
                Rastersh(c, k+1:end) = Raster(c, 1:end-k);
            elseif k < 0
                Rastersh(c, 1:end+k) = Raster(c, 1-k:end);
            else
                Rastersh(c, :) = Raster(c, :);
            end
        end
        
        MActsh = zeros(1, Nz - synchronous_frames);
        for i = 1:(Nz - synchronous_frames)
            MActsh(i) = sum(max(Rastersh(:, i:i+synchronous_frames), [], 2));
        end
        Sumactsh(:, n) = MActsh;
    end

    % Calcul du seuil bas
    Toutes_valeurs_shuffles = Sumactsh(:);
    Valeurs_actives_shuffles = Toutes_valeurs_shuffles(Toutes_valeurs_shuffles > 0);
    Valeurs_triees = sort(Valeurs_actives_shuffles);
    Valeurs_lissees = movmean(Valeurs_triees, 3);
    
    sce_n_cells_threshold = ceil(prctile(Valeurs_lissees, 90));
    
    % Garde-fous seuil bas
    min_cells_absolute = 5;
    max_cells_absolute = 15;
    sce_n_cells_threshold = max(sce_n_cells_threshold, min_cells_absolute);
    sce_n_cells_threshold = min(sce_n_cells_threshold, max_cells_absolute);
    opts.sce_n_cells_threshold = sce_n_cells_threshold;
    
    % =====================================================================
    % 4. DÉTECTION DES PICS AVEC DOUBLE BORNE (MIN ET MAX)
    % =====================================================================
    [pks, locs] = findpeaks(MAct, ...
        'MinPeakHeight', sce_n_cells_threshold, ...
        'MinPeakDistance', MinPeakDistancesce);
    
    % Filtrage du plafond maximum (<= max_cells_allowed)
    valid_peaks = (pks <= max_cells_allowed);
    TRace = locs(valid_peaks);
    NRace = length(TRace);

    % =====================================================================
    % 5. EXTRACTION ET CONTRÔLE DE SÉCURITÉ SUR RACE
    % =====================================================================
 Race = zeros(NCell, NRace);       
    RasterRace = zeros(NCell, Nz);
    
    for i = 1:NRace
        idx_start = TRace(i);
        idx_end   = min(Nz, TRace(i) + synchronous_frames); 
        
        % Extraction de l'identité des cellules actives sur la fenêtre
        Race(:, i) = max(Raster(:, idx_start:idx_end), [], 2);    
        
        % Placement dans RasterRace
        for c = 1:NCell
            if Race(c, i) == 1
                local_frames = idx_start:idx_end;
                active_frame = local_frames(Raster(c, local_frames) == 1);
                
                if ~isempty(active_frame)
                    RasterRace(c, active_frame(1)) = 1; 
                end
            end
        end
    end

    disp(['Ncells= ' num2str(NCell) ' | Sum transient: ' num2str(sum(SumAct)) ' | nSCE détectés : ' num2str(NRace)]);
end