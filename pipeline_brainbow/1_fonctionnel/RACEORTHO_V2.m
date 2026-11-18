function R = RACEORTHO_V2(Race, IDX2, NClOK, pValues_cells, config_asm, graine)
% =================================================================
% RACEORTHO_V2 : assemblées qui peuvent partager des cellules
% (11 nov., demandé par l'utilisateur ; version 1 = RACEORTHO_FDR2)
% =================================================================
% Même chose que RACEORTHO_FDR2, à une différence près : une cellule
% significative dans plusieurs clusters de SCE reste membre de TOUS ces
% clusters. RACEORTHO_FDR2 la mettait dans le seul cluster où sa
% participation était la plus forte (règle héritée de Malvache 2016 :
% assemblées disjointes). Sur 444175_221121_plane1, K = 10 : 22 cellules
% sur 126 étaient concernées, 26 appartenances sur 152 retirées, et une
% assemblée (9 membres -> 3) supprimée.
%
% Calculé SANS refaire le k-means : à partir des clusters de SCE (IDX2,
% NClOK) et des p des cellules (pValues_cells) enregistrés par l'étape 1.
% Étapes identiques à RACEORTHO_FDR2 :
%   1. membres : FDR (config_asm.FDR_Q, famille config_asm.fdr_famille) sur
%      pValues_cells, et au moins config_asm.min_participations SCE du cluster ;
%   2. clusters d'au plus MIN_CELL_SIZE membres retirés ; PAS de règle du
%      maximum ;
%   3. activations : même modèle nul (config_asm.null_cellules) ; dans le
%      modèle 'pondere', une cellule tirée compte pour TOUTES ses assemblées ;
%      FDR config_asm.FDR_Q_activations par SCE ; au moins
%      min_activations_asm activations.
% Sorties : mêmes champs que RACEORTHO_FDR2 (assemblyortho peut maintenant
% contenir une même cellule dans plusieurs assemblées ; CellCl_final aussi),
% plus :
%   assemblee_principale : 1 x NCell, l'assemblée où la cellule participe le
%       plus (0 = aucune) ; sert SEULEMENT à trier le raster (x1) ;
%   n_appartenances : 1 x NCell, nombre d'assemblées de chaque cellule ;
%   version : 'v2'.

[NCell, NRace] = size(Race);
NCl_input = NClOK;
FDR_Q = config_asm.FDR_Q;
FDR_Q_act = FDR_Q;
if isfield(config_asm, 'FDR_Q_activations') && ~isempty(config_asm.FDR_Q_activations)
    FDR_Q_act = config_asm.FDR_Q_activations;
end
MIN_CELL_SIZE = config_asm.MIN_CELL_SIZE;
null_cellules = 'pondere';
if isfield(config_asm, 'null_cellules'), null_cellules = config_asm.null_cellules; end
NShuf = config_asm.NShuf;
fdr_famille = 'globale';
if isfield(config_asm, 'fdr_famille'), fdr_famille = config_asm.fdr_famille; end
min_part = 1;
if isfield(config_asm, 'min_participations'), min_part = config_asm.min_participations; end
min_act = 1;
if isfield(config_asm, 'min_activations_asm'), min_act = config_asm.min_activations_asm; end

if NCl_input < 1 || isempty(pValues_cells) || size(pValues_cells, 1) ~= NCell
    R = sortie_vide(NCell, NRace, 'aucun cluster de SCE significatif (ou p des cellules absents)');
    return
end
IDX2 = IDX2(:)';

%% Étape 1 : membres (comme RACEORTHO_FDR2, à partir des p enregistrés)
CellP = zeros(NCell, NCl_input);
CellR = zeros(NCell, NCl_input);
for i = 1:NCl_input
    CellP(:, i) = sum(Race(:, IDX2 == i), 2);
    n_i = sum(IDX2 == i);
    if n_i > 0, CellR(:, i) = CellP(:, i) / n_i; end
end
Nact = sum(Race, 2);
pValues_cells = pValues_cells(:, 1:NCl_input);
pValues_cells(Nact == 0, :) = 1;
famille = Nact > 0;
CellCl_logical = false(NCell, NCl_input);
switch fdr_famille
    case 'globale'
        p_fam = pValues_cells(famille, :);
        sig = bh_fdr_mask(p_fam(:)', FDR_Q);
        CellCl_logical(famille, :) = reshape(sig, size(p_fam));
    case 'par_cellule'
        for j = find(famille)'
            CellCl_logical(j, :) = bh_fdr_mask(pValues_cells(j, :), FDR_Q);
        end
    otherwise
        error('fdr_famille doit valoir ''globale'' ou ''par_cellule''.');
end
CellCl_logical = CellCl_logical & (CellP >= min_part);
CellCl = double(CellCl_logical');            % [NCl_input x NCell]

%% Étape 2 : candidats (taille seulement ; pas de règle du maximum)
taille_cl = sum(CellCl, 2);
cand_cl = find(taille_cl > MIN_CELL_SIZE)';
NCl_cand = numel(cand_cl);
if NCl_cand == 0
    R = sortie_vide(NCell, NRace, 'aucune assemblée candidate');
    R.pValues_cells = pValues_cells;
    return
end
CellCl_cand = CellCl(cand_cl, :);             % [NCl_cand x NCell], chevauchements permis
assembly_candidates = cell(1, NCl_cand);
for i = 1:NCl_cand
    assembly_candidates{i} = find(CellCl_cand(i, :));
end
M = CellCl_cand';                              % [NCell x NCl_cand]

%% Étape 3 : activations
RCl = CellCl_cand * Race;                      % membres actifs dans chaque SCE
Nrnd_SCE = sum(Race, 1);
switch null_cellules
    case 'uniforme'
        taille_cand = sum(CellCl_cand, 2);
        pValues_race = hygecdf(RCl - 1, NCell, repmat(taille_cand, 1, NRace), ...
            repmat(Nrnd_SCE, NCl_cand, 1), 'upper');
    case 'pondere'
        activite = Nact';
        pValues_race = ones(NCl_cand, NRace);
        parfor j = 1:NRace
            colonne = ones(NCl_cand, 1);
            if Nrnd_SCE(j) > 0
                fixer_flux_aleatoire(graine + 1, j, NRace);
                RClr = tirages_ponderes_membres(activite, Nrnd_SCE(j), NShuf, M);
                colonne = (sum(RClr >= RCl(:, j), 2) + 1) / (NShuf + 1);
            end
            pValues_race(:, j) = colonne;
        end
    otherwise
        error('null_cellules doit valoir ''pondere'' ou ''uniforme''.');
end
PCl = false(NCl_cand, NRace);
for j = 1:NRace
    PCl(:, j) = bh_fdr_mask(pValues_race(:, j)', FDR_Q_act)';
end
PCl = double(PCl);
keepAsm = sum(PCl, 2) >= min_act;
assemblyortho   = assembly_candidates(keepAsm);
cluster_origine = cand_cl(keepAsm);
PCl_final       = PCl(keepAsm, :);
CellCl_final    = CellCl_cand(keepAsm, :);
NCl = sum(keepAsm);
if NCl == 0
    R = sortie_vide(NCell, NRace, 'aucune assemblée active dans assez de SCE');
    R.pValues_cells = pValues_cells;
    return
end

%% Assemblée principale de chaque cellule (pour le tri du raster seulement)
principale = zeros(1, NCell);
part = CellR(:, cluster_origine)';             % [NCl x NCell]
part(CellCl_final == 0) = -Inf;
[mx, k] = max(part, [], 1);
principale(isfinite(mx)) = k(isfinite(mx));
P_princ = zeros(NCl, NCell);
for c = 1:NCl, P_princ(c, principale == c) = 1; end

%% Tri des événements et des cellules (comme RACEORTHO_FDR2)
Cl0 = find(sum(PCl_final, 1) == 0);
Cl1 = find(sum(PCl_final, 1) == 1);
Cl2 = find(sum(PCl_final, 1) == 2);
Cl3 = find(sum(PCl_final, 1) >= 3);
Bin = 2.^(0:NCl-1);
[~, x01] = sort(Bin * PCl_final(:, Cl1)); Cl1 = Cl1(x01);
[~, x02] = sort(Bin * PCl_final(:, Cl2)); Cl2 = Cl2(x02);
[~, x03] = sort(Bin * PCl_final(:, Cl3)); Cl3 = Cl3(x03);
RList = [Cl0, Cl1, Cl2, Cl3];
[~, x1] = sort(Bin * P_princ);

R = struct('assemblyortho', {assemblyortho}, 'cluster_origine', cluster_origine, ...
    'PCl_final', PCl_final, 'CellCl_final', CellCl_final, ...
    'pValues_cells', pValues_cells, 'pValues_race', pValues_race, ...
    'x1', x1, 'RList', RList, 'statut', 'ok', ...
    'assemblee_principale', principale, 'n_appartenances', sum(CellCl_final, 1), 'version', 'v2');
end


function comptes = tirages_ponderes_membres(poids, n, NShuf, M)
% comme tirages_ponderes, mais une cellule tirée compte pour toutes les
% assemblées dont elle est membre (M : nCellules x nAssemblées, 0/1)
cles = rand(NShuf, numel(poids)) .^ (1 ./ poids);
[~, idx] = maxk(cles, n, 2);                    % NShuf x n
nG = size(M, 2);
comptes = zeros(nG, NShuf);
for g = 1:nG
    Mg = M(:, g);
    comptes(g, :) = sum(reshape(Mg(idx), size(idx)), 2)';
end
end


function R = sortie_vide(NCell, NRace, raison)
R = struct('assemblyortho', {cell(1, 0)}, 'cluster_origine', zeros(1, 0), ...
    'PCl_final', zeros(0, NRace), 'CellCl_final', zeros(0, NCell), ...
    'pValues_cells', [], 'pValues_race', [], ...
    'x1', 1:NCell, 'RList', 1:NRace, 'statut', raison, ...
    'assemblee_principale', zeros(1, NCell), 'n_appartenances', zeros(1, NCell), 'version', 'v2');
end
