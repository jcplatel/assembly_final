function [Race,SumAct,MAct,NRace,RasterRace,TRace, Raster,max_cells_allowed] = SCE_detection(Raster,opts,bad_frames) 
    % Note : J'ai ajouté 'Raster' en output, car il sera modifié.

    sce_n_cells_threshold = opts.sce_n_cells_threshold;
    MinPeakDistancesce = opts.MinPeakDistancesce;
    synchronous_frames = opts.synchronous_frames;
    [NCell, Nz] = size(Raster);
    
    % Optionnel: Sécurité absolue, forcer les bad frames à 0 dans le Raster dès le début
    if exist('bad_frames','var') && ~isempty(bad_frames)
         Raster(:, bad_frames) = 0;
    end
    
    SumAct = sum(Raster, 1);
    
    %% NOUVEAU : Nettoyage préventif des événements aberrants
    % On cherche directement les moments où Raster est trop plein, avant même findpeaks
    max_cells_allowed = max(round(NCell * 0.20), 100) ; 
    %max_cells_allowed =  150 ; 
    toxic_frames = SumAct > max_cells_allowed;
    
    if any(toxic_frames)
        fprintf('Nettoyage Raster : %d frames aberrantes (> %d cellules) mises à zéro.\n', ...
                sum(toxic_frames), max_cells_allowed);
        % Dilatation (padding) pour nettoyer autour de l'artefact
        toxic_frames_padded = conv(double(toxic_frames), [1 1 1], 'same') > 0;
        % On efface ces frames du Raster !
        Raster(:, toxic_frames_padded) = 0;
        
        % On doit recalculer SumAct puisque Raster a changé
        SumAct = sum(Raster, 1); 
    end
 
    
    %%   
%%%shuffling to find threshold for number of cell for sce detection
% --- Paramètres ---
NShfl = 100;                        
max_jitter = 30;                    % Décalage maximum (+/- 30 frames)
% 'synchronous_frames' (ex: 3 pour 100ms)

% Nz est le nombre total de frames de votre matrice Raster (Rest uniquement)
Sumactsh = zeros(Nz - synchronous_frames, NShfl);   

disp('Calcul du seuil statistique avec Local Jittering sur Rest...');

for n = 1:NShfl
    Rastersh = zeros(NCell, Nz); 
    
    % Local Jittering indépendant (Zéro-Padding aux extrémités)
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
    
    % Calcul des co-activations
    MActsh = zeros(1, Nz - synchronous_frames);
    for i = 1:(Nz - synchronous_frames)
        MActsh(i) = sum(max(Rastersh(:, i:i+synchronous_frames), [], 2));
    end
    
    Sumactsh(:, n) = MActsh;
end

% --- Détermination du seuil de SCE ---
% Utilisation de la méthode rigoureuse des maximums pour contrer les artefacts de suture
max_per_shuffle = max(Sumactsh, [], 1); 
sce_n_cells_threshold = ceil(prctile(max_per_shuffle, 95));

fprintf('--> Seuil final retenu pour un SCE : %d cellules\n', sce_n_cells_threshold);

    %% Détection normale sur un Raster propre
    [pks, locs] = findpeaks(SumAct, 'MinPeakHeight', sce_n_cells_threshold, 'MinPeakDistance', MinPeakDistancesce);
    
    % On met à jour TRace
    TRace = locs;
    NRace = length(TRace);
    
    % ===== CRÉER RACE MATRIX =====
    MAct = zeros(1, Nz - synchronous_frames);
    for i = 1:Nz - synchronous_frames
        MAct(i) = sum(max(Raster(:, i:i+synchronous_frames), [], 2));
    end
    
    Race = false(NCell, NRace);
    RasterRace = zeros(NCell, Nz);
    for i = 1:NRace
        idx_start = max(1, TRace(i) - 1);
        idx_end = min(Nz, TRace(i) + 2);
        Race(:, i) = max(Raster(:, idx_start:idx_end), [], 2);
        RasterRace(Race(:, i) == 1, TRace(i)) = 1;
    end
    Race = double(Race);
end