function results = biclustering_assemblies(Race, DFF0, TRace, options)
% BICLUSTERING_ASSEMBLIES
% Analyse biclustering (co-clustering cellules × SCE) sur la matrice Race
% et extraction de métriques spatio-temporelles + amplitudes des transients.
%
% INPUTS
%   Race    : [NCell × NSCE] matrice binaire (1 = cellule active dans le SCE)
%   DFF0    : [NCell × Nz]   traces de fluorescence (Delta F/F0)
%   TRace   : [1 × NSCE]     indices temporels des SCE (frames)
%   options : struct optionnel :
%       .nCellClusters : nombre de clusters de cellules   (par défaut 5)
%       .nSCEClusters  : nombre de clusters de SCE        (par défaut 5)
%       .nComponents   : nb de composantes SVD (spectral) (par défaut 5)
%       .win_peak      : fenêtre micro pour amplitude DFF0 (par défaut 2)
%
% OUTPUT
%   results : struct avec champs :
%       .cellLabels       : cluster ID pour chaque cellule
%       .sceLabels        : cluster ID pour chaque SCE
%       .Event_Purity     : [1 × nBiclusters]
%       .Cellular_Fidelity: [1 × nBiclusters]
%       .AmpTable         : table des transients (ClusterID, SCE, Cell, Frame, Amp)
%       .options          : les options utilisées
%
% Auteur : toi ;)

%% --------------------- 0. Options par défaut -----------------------------
if nargin < 4
    options = struct();
end

if ~isfield(options, 'nCellClusters'); options.nCellClusters = 5; end
if ~isfield(options, 'nSCEClusters');  options.nSCEClusters  = 5; end
if ~isfield(options, 'nComponents');   options.nComponents   = 5; end
if ~isfield(options, 'win_peak');      options.win_peak      = 2; end

nCellClusters = options.nCellClusters;
nSCEClusters  = options.nSCEClusters;
nComponents   = options.nComponents;
win_peak      = options.win_peak;

%% --------------------- 1. Préparation de la matrice -----------------------
[NCell, NSCE] = size(Race);

% On travaille sur une version centrée / normalisée de Race
Race_double = double(Race);

% Centre par colonne (SCE) pour enlever les différences de "taille" de SCE
Race_col_centered = Race_double - mean(Race_double, 1);

% Optionnel : normalisation des lignes (cellules) par leur norme
row_norm = sqrt(sum(Race_col_centered.^2, 2));
row_norm(row_norm == 0) = 1; % éviter division par zéro
Race_norm = Race_col_centered ./ row_norm;

%% --------------------- 2. Spectral biclustering (SVD) --------------------
% Inspiré des approches de biclustering spectrale : on utilise U (cellules)
% et V (SCE) pour faire un k-means séparé sur lignes et colonnes.[web:401][web:396]

% SVD tronquée (si nComponents < min(size))
nComponents = min([nComponents, NCell, NSCE]);
[U, S, V] = svds(Race_norm, nComponents);  % U : NCell×C, V : NSCE×C

% On peut utiliser U pour clusteriser les cellules, V pour les SCE
cellFeatures = U;          % NCell × nComponents
sceFeatures  = V;          % NSCE × nComponents

% K-means sur cellules
cellLabels = kmeans(cellFeatures, nCellClusters, ...
                    'Replicates', 10, 'MaxIter', 1000);

% K-means sur SCE
sceLabels = kmeans(sceFeatures, nSCEClusters, ...
                   'Replicates', 10, 'MaxIter', 1000);

%% --------------------- 3. Construction des biclusters --------------------
% Un bicluster est défini par (clusterCell = c, clusterSCE = s).
% On parcourt toutes les combinaisons possibles.
%
% On va calculer :
%   - Event_Purity(c,s)
%   - Cellular_Fidelity(c,s)
%
nBiclusters = nCellClusters * nSCEClusters;
Event_Purity      = NaN(1, nBiclusters);
Cellular_Fidelity = NaN(1, nBiclusters);

biclusterIndex = @(c,s) (c-1)*nSCEClusters + s;  % index linéaire

for c = 1:nCellClusters
    cells_c = find(cellLabels == c);  % cellules du cluster c
    
    if isempty(cells_c)
        continue;
    end
    
    for s = 1:nSCEClusters
        idx = biclusterIndex(c,s);
        
        sce_s = find(sceLabels == s); % SCE du cluster s
        
        if isempty(sce_s)
            continue;
        end
        
        % --- CALCUL VERTICAL (Event Purity) ---
        % Tirs totaux de toutes les cellules pendant les SCE (c,s)
        tirs_totaux_colonne = sum(Race(:, sce_s), 'all');
        
        % Tirs des cellules de c pendant ces SCE
        tirs_in_box = sum(Race(cells_c, sce_s), 'all');
        
        if tirs_totaux_colonne > 0
            Event_Purity(idx) = tirs_in_box / tirs_totaux_colonne;
        end
        
        % --- CALCUL HORIZONTAL (Cellular Fidelity) ---
        % Tirs des cellules de c dans TOUS les SCE (tous s)
        sce_all_valid = (1:NSCE); % ici : tous les SCE; tu peux filtrer si besoin
        tirs_totaux_ligne = sum(Race(cells_c, sce_all_valid), 'all');
        
        if tirs_totaux_ligne > 0
            Cellular_Fidelity(idx) = tirs_in_box / tirs_totaux_ligne;
        end
    end
end

%% --------------------- 4. Extraction des amplitudes DFF0 -----------------
% On extrait les amplitudes des transients au sein de chaque bicluster,
% un peu comme ce que tu fais déjà pour les clusters consensus.

Amplitudes_Biclusters = []; % [CellCluster, SCECluster, SCE_Index, CellID, PeakFrame, Amp]

for c = 1:nCellClusters
    cells_c = find(cellLabels == c);
    if isempty(cells_c); continue; end
    
    for s = 1:nSCEClusters
        sce_s = find(sceLabels == s);
        if isempty(sce_s); continue; end
        
        for k = 1:length(sce_s)
            current_sce_idx = sce_s(k);      % index dans Race (colonne)
            sce_frame = TRace(current_sce_idx); % frame du SCE
            
            % Cellules actives dans CE SCE précis
            active_cells = find(Race(:, current_sce_idx) > 0);
            
            % Fenêtre locale autour du SCE (ici [-1,+2] comme chez toi)
            idx_start = max(1, sce_frame - 1);
            idx_end   = min(size(DFF0, 2), sce_frame + 2);
            local_frames = idx_start:idx_end;
            
            for iCell = 1:length(active_cells)
                cell_idx = active_cells(iCell);
                
                % Frame du pic dans Raster (si tu as Raster disponible)
                % Ici on suppose que Race reflète déjà le pic, donc on prend sce_frame
                exact_peak_frame = sce_frame;
                
                % Micro-fenêtre pour le max de DFF0 autour du pic
                f_start = max(1, exact_peak_frame - win_peak);
                f_end   = min(size(DFF0, 2), exact_peak_frame + win_peak);
                
                [max_amp, ~] = max(DFF0(cell_idx, f_start:f_end));
                
                % Stockage
                Amplitudes_Biclusters = [Amplitudes_Biclusters; ...
                    c, s, current_sce_idx, cell_idx, exact_peak_frame, max_amp];
            end
        end
    end
end

if ~isempty(Amplitudes_Biclusters)
    AmpTable = array2table(Amplitudes_Biclusters, ...
        'VariableNames', {'CellCluster', 'SCECluster', 'SCE_Index', ...
                          'CellID', 'PeakFrame', 'Amplitude'});
else
    AmpTable = table();
end

%% --------------------- 5. Rassemblement des résultats --------------------
results = struct();
results.cellLabels        = cellLabels;
results.sceLabels         = sceLabels;
results.Event_Purity      = Event_Purity;
results.Cellular_Fidelity = Cellular_Fidelity;
results.AmpTable          = AmpTable;
results.options           = options;

end