%% ========================================================================
% EXPORT_BRAINBOW_SUMMARY
% Objectif : parcourir un dossier et ses sous-dossiers, trouver tous les
% fichiers brainbow.mat, charger les variables nécessaires et reconstruire :
%   1) le fichier Excel de synthèse par session : analysis_brainbow_ALL.xlsx
%   2) un fichier Excel par assemblée (1 ligne par assemblée, toutes
%      sessions confondues), incluant les p-values hypergéométriques
%      par couleur (p_values_hyperG_nobl) : analysis_assemblies_pvalues.xlsx
%
% Toutes les lignes sont accumulées en mémoire puis écrites une seule fois
% à la fin (plus rapide et plus sûr que réécrire le xlsx à chaque fichier).
% =========================================================================

clear; clc; close all;

%% 1. PARAMETRES
rootFolder    = 'E:\Data\Aurelie\analysis\Final\assembly\nocues\';                                              % <-- à adapter
outputExcel   = fullfile(rootFolder, 'analysis_brainbow_ALL_rebuilt.xlsx');
outputExcelPV = fullfile(rootFolder, 'analysis_assemblies_pvalues_rebuilt.xlsx');
useParfor     = true;                                                % mettre false si pas de Parallel Computing Toolbox

%% 2. VARIABLES A CHARGER DANS CHAQUE brainbow.mat
requiredVars = {'daytime','identifier','NRace','freqsce','NRaceOK','freqsceok', ...
                'successrate','sCl','sClOK','NClOK','mean_sClOK','assemblyortho', ...
                'mean_cell_assembly','NCell','Nz','sampling_rate','WinRest', ...
                'average_speed','max_speed','nb_lap','sce_n_cells_threshold', ...
                'mean_distance_assembly','meandistance_gcamp','kmeans_surrogate', ...
                'file_num','NClini','color_assembly','p_values_chi2','p_values_hyperG_nobl'};

%% 3. EN-TETES DES DEUX FICHIERS EXCEL
header_global = {'date', 'identifier', 'SCE', 'sce freq', 'SCEok', 'SCEok freq', ...
    'Assemblies', ' % colored', 'sil all clusters','SClok','SClokfinal', 'cell/assembly', ...
    'mean cell/assembly', 'Ncell', 'frames','sampling rate', 'nb frame not running', ...
    '% not running', 'average speed', 'max speed', 'nblap', 'sce_n_cells_threshold', ...
    'mean distance assemblies', 'mean distance all cells', 'n_kmean_surrogate', 'file_num','NCl_ini'};
nCols = numel(header_global);
t
header_pvalues = {'identifier', 'NCl_ini', 'assembly_num', 'sClOK_assembly', ...
    'color_assembly', 'chi2_p', 'p_hyperG_R', 'p_hyperG_G', 'p_hyperG_B', ...
    'p_hyperG_RG', 'p_hyperG_RB', 'p_hyperG_GB'};
nColsPV = numel(header_pvalues);

%% 4. RECHERCHE RECURSIVE DE TOUS LES brainbow.mat
fprintf('Recherche des fichiers brainbow.mat dans %s ...\n', rootFolder);
tousLesFichiers = dir(fullfile(rootFolder, '**', 'brainbow.mat'));
if isempty(tousLesFichiers)
    error('Aucun fichier brainbow.mat trouvé dans %s', rootFolder);
end
Nfichiers = length(tousLesFichiers);
fprintf('%d fichiers trouvés.\n', Nfichiers);

%% 5. PRE-ALLOCATION
dataRows    = cell(Nfichiers, nCols);
rowValid    = false(Nfichiers, 1);
pvCellArray = cell(Nfichiers, 1);   % chaque élément = matrice cell (nbassemblies x nColsPV) ou []

%% 6. EXTRACTION (parfor ou for classique)
if useParfor
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool;
    end
    D = parallel.pool.DataQueue;
    hwait = waitbar(0, 'Extraction en cours...', 'Name', 'Export Brainbow');
    afterEach(D, @(~) updateWaitbar(hwait, Nfichiers));
end

if useParfor
    parfor i = 1:Nfichiers
        [rowValid(i), dataRows(i,:), pvCellArray{i}] = extractRow(tousLesFichiers(i), requiredVars, nCols, nColsPV);
        send(D, i);
    end
    close(hwait);
else
    for i = 1:Nfichiers
        [rowValid(i), dataRows(i,:), pvCellArray{i}] = extractRow(tousLesFichiers(i), requiredVars, nCols, nColsPV);
        if mod(i,10)==0
            fprintf('%d / %d fichiers traités\n', i, Nfichiers);
        end
    end
end

%% 7. FILTRAGE ET ECRITURE DU RESUME PAR SESSION
dataRows = dataRows(rowValid, :);
fprintf('%d / %d fichiers exploitables (résumé par session).\n', size(dataRows,1), Nfichiers);

if isempty(dataRows)
    error('Aucune ligne valide à écrire pour le résumé par session.');
end

fullTable = [header_global; dataRows];
writecell(fullTable, outputExcel, 'Range', 'A1');
fprintf('Fichier Excel écrit : %s (%d lignes)\n', outputExcel, size(dataRows,1));

%% 8. CONCATENATION ET ECRITURE DES P-VALUES PAR ASSEMBLEE
pvCellArray = pvCellArray(~cellfun(@isempty, pvCellArray));
if isempty(pvCellArray)
    warning('Aucune p-value hypergéométrique trouvée (p_values_hyperG_nobl absent des fichiers). Fichier non créé.');
else
    allPvRows = vertcat(pvCellArray{:});
    fullTablePV = [header_pvalues; allPvRows];
    writecell(fullTablePV, outputExcelPV, 'Range', 'A1');
    fprintf('Fichier Excel écrit : %s (%d lignes, 1 par assemblée)\n', outputExcelPV, size(allPvRows,1));
end


%% =========================================================================
% FONCTIONS LOCALES (doivent être en fin de fichier)
% =========================================================================

function [ok, row, pvRows] = extractRow(fileEntry, requiredVars, nCols, nColsPV)
    ok = false;
    row = cell(1, nCols);
    pvRows = [];

    matFilePath = fullfile(fileEntry.folder, fileEntry.name);

    try
        data = load(matFilePath, requiredVars{:});
    catch
        return; % fichier corrompu ou illisible
    end

    % --- daytime : fallback sur la date de modification du fichier si absente ---
    if isfield(data, 'daytime') && ~isempty(data.daytime)
        daytime = data.daytime;
    else
        daytime = datestr(fileEntry.datenum, 'yy_mm_dd_HH_MM_SS');
    end

    % --- identifier : obligatoire ---
    if ~isfield(data, 'identifier') || isempty(data.identifier)
        return;
    end
    identifier = data.identifier;

    % --- assemblyortho -> writeassemblyortho et idortho ---
    if isfield(data, 'assemblyortho') && ~isempty(data.assemblyortho)
        writeassemblyortho = length(data.assemblyortho);
        idortho = cellfun(@length, data.assemblyortho);
    else
        writeassemblyortho = 0;
        idortho = [];
    end

    % --- NClini : fallback sur le nom du dossier (ex: ...\K7\) si absent ---
    if isfield(data, 'NClini') && ~isempty(data.NClini)
        NClini = data.NClini;
    else
        parts = strsplit(fileEntry.folder, filesep);
        kFolderName = parts{end};
        knumstr = regexp(kFolderName, '\d+', 'match');
        if isempty(knumstr)
            return;
        end
        NClini = str2double(knumstr{1});
    end

    % --- variables numériques simples avec valeur par défaut NaN si absentes ---
    getval = @(name) getFieldOrNaN(data, name);

    NRace                  = getval('NRace');
    freqsce                = getval('freqsce');
    NRaceOK                = getval('NRaceOK');
    freqsceok              = getval('freqsceok');
    successrate            = getval('successrate');
    NClOK                  = getval('NClOK');
    mean_sClOK             = getval('mean_sClOK');
    mean_cell_assembly     = getval('mean_cell_assembly');
    NCell                  = getval('NCell');
    Nz                     = getval('Nz');
    sampling_rate          = getval('sampling_rate');
    average_speed          = getval('average_speed');
    max_speed              = getval('max_speed');
    nb_lap                 = getval('nb_lap');
    sce_n_cells_threshold  = getval('sce_n_cells_threshold');
    mean_distance_assembly = getval('mean_distance_assembly');
    meandistance_gcamp     = getval('meandistance_gcamp');
    kmeans_surrogate       = getval('kmeans_surrogate');
    file_num               = getval('file_num');

    if isfield(data, 'sCl') && ~isempty(data.sCl)
        sCl = data.sCl;
        sil_all = mean(sCl, 'omitnan');
        if ~isnan(NClOK) && NClOK > 0 && NClOK <= length(sCl)
            sil_ok = mean(sCl(1:NClOK), 'omitnan');
        else
            sil_ok = NaN;
        end
    else
        sil_all = NaN;
        sil_ok = NaN;
    end

    if isfield(data, 'WinRest') && ~isempty(data.WinRest)
        nFrameNotRunning = length(data.WinRest);
        if ~isnan(Nz) && Nz > 0
            pctNotRunning = nFrameNotRunning / Nz * 100;
        else
            pctNotRunning = NaN;
        end
    else
        nFrameNotRunning = NaN;
        pctNotRunning = NaN;
    end

    % --- Ligne du résumé par session (même ordre que header_global) ---
    row = {daytime, identifier, NRace, freqsce, NRaceOK, freqsceok, ...
           writeassemblyortho, successrate, sil_all, sil_ok, mean_sClOK, num2str(idortho), ...
           mean_cell_assembly, NCell, Nz, sampling_rate, nFrameNotRunning, pctNotRunning, ...
           average_speed, max_speed, nb_lap, sce_n_cells_threshold, mean_distance_assembly, ...
           meandistance_gcamp, kmeans_surrogate, file_num, NClini};

    ok = true;

    % --- Table des p-values par assemblée (1 ligne par assemblée) ---
    if isfield(data, 'p_values_hyperG_nobl') && ~isempty(data.p_values_hyperG_nobl)
        pvMatrix = data.p_values_hyperG_nobl;      % nbassemblies x 6
        nbAssemblies = size(pvMatrix, 1);

        if isfield(data, 'sClOK') && ~isempty(data.sClOK)
            sClOKvec = nan(nbAssemblies, 1);
            nb_scl = min(length(data.sClOK), nbAssemblies);
            sClOKvec(1:nb_scl) = data.sClOK(1:nb_scl);
        else
            sClOKvec = nan(nbAssemblies, 1);
        end

        if isfield(data, 'color_assembly') && ~isempty(data.color_assembly)
            colorAssembly = repmat({'Not significant'}, nbAssemblies, 1);
            nb_col = min(length(data.color_assembly), nbAssemblies);
            for c_idx = 1:nb_col
                if ~isempty(data.color_assembly{c_idx})
                    colorAssembly{c_idx} = data.color_assembly{c_idx};
                end
            end
        else
            colorAssembly = repmat({'Not significant'}, nbAssemblies, 1);
        end

        if isfield(data, 'p_values_chi2') && ~isempty(data.p_values_chi2)
            chi2p = nan(nbAssemblies, 1);
            nb_chi2 = min(length(data.p_values_chi2), nbAssemblies);
            chi2p(1:nb_chi2) = data.p_values_chi2(1:nb_chi2);
        else
            chi2p = nan(nbAssemblies, 1);
        end

        pvRows = cell(nbAssemblies, nColsPV);
        for a = 1:nbAssemblies
            pvRows(a, :) = {identifier, NClini, a, sClOKvec(a), colorAssembly{a}, chi2p(a), ...
                            pvMatrix(a,1), pvMatrix(a,2), pvMatrix(a,3), ...
                            pvMatrix(a,4), pvMatrix(a,5), pvMatrix(a,6)};
        end
    end
end

function v = getFieldOrNaN(s, fieldName)
    if isfield(s, fieldName) && ~isempty(s.(fieldName))
        v = s.(fieldName);
    else
        v = NaN;
    end
end

function updateWaitbar(h, total)
    persistent count
    if isempty(count), count = 0; end
    count = count + 1;
    if isgraphics(h)
        waitbar(count/total, h, sprintf('Analyse %d/%d', count, total));
    end
    if count >= total, count = []; end % reset pour un prochain run
end
