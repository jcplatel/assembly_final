function A = SCE_clustering2(Race, TRace, Raster, DFF0, Acttmp2, PromALL, NClini, config_asm, synchronous_frames, graine)
%% =========================================================================
%  SCE_clustering2 : version fonction de SCE_clustering
%  -------------------------------------------------------------------------
%  Entrée : Race (cellules x SCE) et le nombre de clusters NClini.
%  Sortie : struct A dont les champs portent les noms de tes anciennes
%  variables (IDX2, sCl, NClOK, assemblyortho, PCl_final, mean_sClOK, ...).
%  save(fichier, '-struct', 'A') les enregistre sous ces noms.
%
%  Changements par rapport à SCE_clustering :
%   - seuil des clusters au choix (config_asm.seuil_clusters) :
%       'rang' (défaut) : le i-e meilleur cluster est comparé au i-e meilleur
%                cluster des surrogates (95e percentile), en s'arrêtant au
%                premier échec (comme l'analyse parallèle de Horn)
%       'max'  : tous comparés au MEILLEUR cluster de chaque surrogate (très
%                strict : sur simulation, perd 1 assemblée sur 2 à 3)
%       'tous' : ancien seuil, 95e percentile de tous les clusters mélangés
%     Sur simulation, aucune règle ne crée d'assemblée sans vraie assemblée
%     (RACEORTHO et le modèle curveball suffisent)
%   - surrogates : mêmes options k-means que les données, même règle MinSCE,
%     modèle nul qui garde l'activité de chaque cellule (curveball)
%   - chaque assemblée garde son cluster de SCE d'origine (cluster_origine) :
%     purity, fidelity et AmpTable utilisent le BON cluster
%   - mean_sClOK = silhouette des clusters retenus, même définition que sCl
%   - NClOK = 1 n'est plus perdu
%   - amplitudes : fenêtre de Race ([TRace, TRace + synchronous_frames]),
%     et proéminence ajoutée (même mesure que la détection)
%  =========================================================================

kmean_iter       = config_asm.kmean_iter;
% réplicats des surrogates : par défaut les mêmes que les données
kmean_iter_surrogate = kmean_iter;
if isfield(config_asm, 'kmean_iter_surrogate') && ~isempty(config_asm.kmean_iter_surrogate)
    kmean_iter_surrogate = config_asm.kmean_iter_surrogate;
end
kmeans_surrogate = config_asm.kmeans_surrogate;
MinSCE           = config_asm.MinSCE;
null_kmeans      = config_asm.null_kmeans;
% Phase en ligne du k-means (config_asm.kmeans_online, défaut true), la même
% pour les données et les surrogates. online_surrogate (facultatif : réglage
% 'ancien_rang' et retrouver_ancien) la change pour les surrogates seuls.
kmeans_online = true;
if isfield(config_asm, 'kmeans_online') && ~isempty(config_asm.kmeans_online)
    kmeans_online = logical(config_asm.kmeans_online);
end
online_surrogate = kmeans_online;
if isfield(config_asm, 'online_surrogate'), online_surrogate = config_asm.online_surrogate; end

[NCell, NRace] = size(Race);
A = struct();
A.NClini = NClini;
A.NRace  = NRace;

% Pas assez de SCE pour faire NClini clusters
if NRace <= NClini
    A.statut = sprintf('trop peu de SCE (%d) pour K = %d', NRace, NClini);
    A.IDX2 = zeros(1, NRace);  A.sCl = [];  A.threshold = NaN;
    A.NClOK = 0;  A.NRaceOK = 0;  A.NCl = 0;
    A.assemblyortho = cell(1, 0);  A.cluster_origine = zeros(1, 0);
    A.PCl_final = zeros(0, NRace);  A.CellCl_final = zeros(0, NCell);
    A.sClOK = [];  A.mean_sClOK = NaN;
    return
end

%% 1. k-means sur les SCE
[IDX2, sCl, M, S, ~, inertia, nSCE_cluster] = kmeansoptnew2(Race, kmean_iter, 'var', NClini, [], MinSCE, graine, ...
    kmeans_online);
NCl = NClini;
[~, x2] = sort(IDX2);      % SCE triés par cluster (figure)

%% 2. Clusters statistiquement significatifs
% Silhouettes des clusters de chaque surrogate, triées par ordre décroissant
if isfield(config_asm, 'seuil_clusters'), seuil_clusters = config_asm.seuil_clusters; else, seuil_clusters = 'rang'; end
p_seuil = 95;          % percentile des surrogates qui sert de seuil
if isfield(config_asm, 'percentile_seuil') && ~isempty(config_asm.percentile_seuil)
    p_seuil = config_asm.percentile_seuil;
end
sCl_rnd = zeros(kmeans_surrogate, NCl);
parfor i = 1:kmeans_surrogate
    fixer_flux_aleatoire(graine + 1, i, kmeans_surrogate);
    scores_du_run = kmeansoptrndnew2(Race, kmean_iter_surrogate, NCl, MinSCE, null_kmeans, online_surrogate);
    sCl_rnd(i, :) = sort(scores_du_run(:)', 'descend');
end
max_rnd = sCl_rnd(:, 1)';
switch seuil_clusters
    case 'rang'     % i-e cluster des données contre i-e cluster des surrogates
        seuils_rang = prctile(sCl_rnd, p_seuil, 1);
        ok = sCl(:)' > seuils_rang;
        NClOK = find(~ok, 1) - 1;                 % on s'arrête au premier échec
        if isempty(NClOK), NClOK = NCl; end
        threshold = seuils_rang(max(NClOK, 1));
    case 'max'      % tous contre le meilleur cluster de chaque surrogate
        threshold = prctile(max_rnd, p_seuil);
        NClOK = sum(sCl > threshold);             % valable car sCl est trié
        seuils_rang = repmat(threshold, 1, NCl);
    case 'tous'     % ancien : tous les clusters des surrogates mélangés
        threshold = prctile(sCl_rnd(:), p_seuil);
        NClOK = sum(sCl > threshold);
        seuils_rang = repmat(threshold, 1, NCl);
    otherwise
        error('seuil_clusters doit valoir ''rang'', ''max'' ou ''tous''.');
end
RaceOK  = Race(:, IDX2 <= NClOK);      % SCE des clusters significatifs
NRaceOK = size(RaceOK, 2);

%% 3. Assemblées (RACEORTHO_FDR2)
R = RACEORTHO_FDR2(Race, IDX2, NClOK, config_asm, graine + 2);
assemblyortho   = R.assemblyortho;
cluster_origine = R.cluster_origine;    % assemblyortho{c} <-> IDX2 == cluster_origine(c)
NCl = numel(assemblyortho);             % nombre FINAL d'assemblées (comme avant)

%% 4. Silhouette des clusters finaux (même définition que sCl et que le seuil)
if NCl >= 1
    sClOK = sCl(cluster_origine);
    mean_sClOK = mean(sClOK);
else
    sClOK = [];
    mean_sClOK = NaN;
end

%% 5. Spécificité spatio-temporelle (Purity & Fidelity)
Event_Purity      = NaN(1, NCl);    % vertical
Cellular_Fidelity = NaN(1, NCl);    % horizontal
valid_sce_mask = ismember(IDX2, cluster_origine);   % SCE des clusters qui ont donné une assemblée
for c = 1:NCl
    cells_c   = assemblyortho{c};
    idx_SCE_c = (IDX2 == cluster_origine(c));       % SCE du cluster de CETTE assemblée
    if sum(idx_SCE_c) > 0 && ~isempty(cells_c)
        tirs_totaux_colonne = sum(Race(:, idx_SCE_c), 'all');
        tirs_in_box         = sum(Race(cells_c, idx_SCE_c), 'all');
        if tirs_totaux_colonne > 0
            Event_Purity(c) = tirs_in_box / tirs_totaux_colonne;
        end
        tirs_totaux_ligne = sum(Race(cells_c, valid_sce_mask), 'all');
        if tirs_totaux_ligne > 0
            Cellular_Fidelity(c) = tirs_in_box / tirs_totaux_ligne;
        end
    end
end
mean_Event_Purity      = mean(Event_Purity, 'omitnan');
mean_Cellular_Fidelity = mean(Cellular_Fidelity, 'omitnan');

%% 6. Amplitudes des transients dans les SCE des assemblées
% [ClusterID (= numéro de l'assemblée), SCE_Index, CellID, PeakFrame, Amplitude (max de DFF0), Prominence]
Amplitudes_SCE_Clustered = zeros(0, 6);
for c = 1:NCl
    idx_SCE_c = find(IDX2 == cluster_origine(c));
    for k = 1:numel(idx_SCE_c)
        current_sce_idx = idx_SCE_c(k);
        sce_frame = TRace(current_sce_idx);
        active_cells = find(Race(:, current_sce_idx) > 0);
        % même fenêtre que celle qui définit Race
        local_frames = sce_frame:min(size(DFF0, 2), sce_frame + synchronous_frames);
        for i = 1:numel(active_cells)
            cell_idx = active_cells(i);
            active_frame_idx = find(Raster(cell_idx, local_frames) == 1, 1);
            if ~isempty(active_frame_idx)
                exact_peak_frame = local_frames(active_frame_idx);
                win_micro = 2;   % le max de DFF0 brut peut être décalé d'une trame
                f_start = max(1, exact_peak_frame - win_micro);
                f_end   = min(size(DFF0, 2), exact_peak_frame + win_micro);
                max_amp = max(DFF0(cell_idx, f_start:f_end));
                prom = PromALL{cell_idx}(Acttmp2{cell_idx} == exact_peak_frame);
                if isempty(prom), prom = NaN; end
                Amplitudes_SCE_Clustered(end+1, :) = ...
                    [c, current_sce_idx, cell_idx, exact_peak_frame, double(max_amp), double(prom(1))]; %#ok<AGROW>
            end
        end
    end
end
AmpTable = array2table(Amplitudes_SCE_Clustered, ...
    'VariableNames', {'ClusterID', 'SCE_Index', 'CellID', 'PeakFrame', 'Amplitude', 'Prominence'});

%% 7. Sortie (noms habituels)
A.statut          = R.statut;
A.IDX2            = IDX2;
A.sCl             = sCl;
A.S               = S;
A.M               = single(M);
A.x2              = x2;
A.inertia         = inertia;
A.nSCE_cluster    = nSCE_cluster;
A.max_rnd         = max_rnd;
A.sCl_rnd         = sCl_rnd;
A.seuils_rang     = seuils_rang;
A.seuil_clusters  = seuil_clusters;
A.threshold       = threshold;
A.NClOK           = NClOK;
A.NRaceOK         = NRaceOK;
A.assemblyortho   = assemblyortho;
A.cluster_origine = cluster_origine;
A.NCl             = NCl;
A.PCl_final       = R.PCl_final;
A.CellCl_final    = R.CellCl_final;
A.pValues_cells   = R.pValues_cells;
A.pValues_race    = R.pValues_race;
A.x1              = R.x1;
A.RList           = R.RList;
A.sClOK           = sClOK;
A.mean_sClOK      = mean_sClOK;
A.Event_Purity           = Event_Purity;
A.Cellular_Fidelity      = Cellular_Fidelity;
A.mean_Event_Purity      = mean_Event_Purity;
A.mean_Cellular_Fidelity = mean_Cellular_Fidelity;
A.AmpTable        = AmpTable;
end
