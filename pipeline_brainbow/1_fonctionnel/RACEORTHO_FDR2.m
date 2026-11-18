function R = RACEORTHO_FDR2(Race, IDX2, NClOK, config_asm, graine)
% =================================================================
% RACEORTHO_FDR2 : version fonction de RACEORTHO_FDR
% =================================================================
% Entrées : Race (cellules x SCE), IDX2 (cluster de chaque SCE, triés par
% silhouette), NClOK (nb de clusters significatifs).
% Sortie : struct R (assemblyortho, cluster_origine, PCl_final,
% CellCl_final, pValues_cells, pValues_race, x1, RList, statut).
%
% Changements par rapport à RACEORTHO_FDR :
%  1. chaque assemblée garde le numéro de son cluster de SCE :
%     assemblyortho{c} <-> IDX2 == cluster_origine(c)
%  2. modèle nul 'pondere' : à l'étape 1, un gros SCE a plus de chances de
%     contenir une cellule ; à l'étape 3, une cellule active a plus de
%     chances d'être tirée qu'une cellule silencieuse.
%     'uniforme' = ancien modèle, calculé exactement (loi hypergéométrique)
%  3. FDR global sur toute la matrice cellules x clusters ('globale')
%  4. une cellule doit participer à au moins min_participations SCE du cluster
%  5. filtre de taille AVANT les cellules « promiscuous », qui ne sont
%     réassignées qu'entre clusters où elles sont significatives
%  6. une assemblée doit être active dans au moins min_activations_asm SCE
%  7. NClOK = 1 est traité (avant : aucune assemblée)

%% --- Paramètres ---
[NCell, NRace] = size(Race);
NCl_input     = NClOK;                    % nb de clusters significatifs
NShuf         = config_asm.NShuf;         % nb de permutations (modèle 'pondere')
FDR_Q         = config_asm.FDR_Q;         % taux de fausses découvertes toléré (membres)
% activations : FDR_Q_activations s'il existe, sinon FDR_Q (comme avant)
FDR_Q_act = FDR_Q;
if isfield(config_asm, 'FDR_Q_activations') && ~isempty(config_asm.FDR_Q_activations)
    FDR_Q_act = config_asm.FDR_Q_activations;
end
MIN_CELL_SIZE = config_asm.MIN_CELL_SIZE; % taille > MIN_CELL_SIZE (5 -> au moins 6 cellules)
null_cellules = config_asm.null_cellules;

if NCl_input < 1
    R = sortie_vide(NCell, NRace, 'aucun cluster de SCE significatif');
    return
end

%% Étape 1 : association cellule-cluster
% 1.1 Données observées
CellP   = zeros(NCell, NCl_input);   % nb de SCE du cluster où la cellule est active
CellR   = zeros(NCell, NCl_input);   % taux de participation
nSCE_cl = zeros(1, NCl_input);
for i = 1:NCl_input
    CellP(:, i) = sum(Race(:, IDX2 == i), 2);
    nSCE_cl(i)  = sum(IDX2 == i);
    if nSCE_cl(i) > 0
        CellR(:, i) = CellP(:, i) / nSCE_cl(i);
    end
end
Nact = sum(Race, 2);                 % nb de SCE de chaque cellule

% 1.2 p-values
switch null_cellules
    case 'uniforme'
        % Ancien modèle : les Nact SCE de la cellule tirés au hasard, tous
        % équiprobables. Calcul exact : P(X >= observé), loi hypergéométrique.
        pValues_cells = hygecdf(CellP - 1, NRace, repmat(nSCE_cl, NCell, 1), ...
            repmat(Nact, 1, NCl_input), 'upper');
    case 'pondere'
        % Les Nact SCE de la cellule tirés avec une probabilité
        % proportionnelle à leur taille (un gros SCE contient plus de cellules).
        taille_SCE = sum(Race, 1);
        pValues_cells = ones(NCell, NCl_input);
        parfor j = 1:NCell
            ligne = ones(1, NCl_input);
            if Nact(j) > 0
                fixer_flux_aleatoire(graine, j, NCell);
                RClr = tirages_ponderes(taille_SCE, Nact(j), NShuf, IDX2, NCl_input);
                ligne = ((sum(RClr >= CellP(j, :)', 2) + 1) / (NShuf + 1))';
            end
            pValues_cells(j, :) = ligne;
        end
    otherwise
        error('null_cellules doit valoir ''pondere'' ou ''uniforme''.');
end
pValues_cells(Nact == 0, :) = 1;

% 1.3 FDR de Benjamini-Hochberg (seulement les cellules actives dans au moins un SCE)
famille = Nact > 0;
CellCl_logical = false(NCell, NCl_input);
switch config_asm.fdr_famille
    case 'globale'      % une seule famille : toute la matrice cellules x clusters
        p_fam = pValues_cells(famille, :);
        sig = bh_fdr_mask(p_fam(:)', FDR_Q);
        CellCl_logical(famille, :) = reshape(sig, size(p_fam));
    case 'par_cellule'  % ancien : une famille par cellule
        for j = find(famille)'
            CellCl_logical(j, :) = bh_fdr_mask(pValues_cells(j, :), FDR_Q);
        end
    otherwise
        error('fdr_famille doit valoir ''globale'' ou ''par_cellule''.');
end
CellCl_logical = CellCl_logical & (CellP >= config_asm.min_participations);
CellCl = double(CellCl_logical');    % [NCl_input x NCell]

%% Étape 2 : formation des assemblées candidates
% 2.1 Taille d'abord : seuls les clusters avec assez de cellules significatives restent
taille_cl = sum(CellCl, 2);
CellCl(taille_cl <= MIN_CELL_SIZE, :) = 0;

% 2.2 Cellules « promiscuous » : assignées au cluster où leur taux de
%     participation est maximal, PARMI ceux où elles sont significatives
A2 = find(sum(CellCl, 1) >= 2);
for i = A2
    sig_i = find(CellCl(:, i));
    [~, k] = max(CellR(i, sig_i));
    CellCl(:, i) = 0;
    CellCl(sig_i(k), i) = 1;
end

% 2.3 Taille vérifiée à nouveau ; chaque candidat garde son cluster d'origine
C0 = {};
cand_cl = [];
for i = 1:NCl_input
    cell_indices = find(CellCl(i, :));
    if length(cell_indices) > MIN_CELL_SIZE
        C0{end+1} = cell_indices; %#ok<AGROW>
        cand_cl(end+1) = i; %#ok<AGROW>
    end
end
assembly_candidates = C0;
NCl_cand = numel(assembly_candidates);
if NCl_cand == 0
    R = sortie_vide(NCell, NRace, 'aucune assemblée candidate');
    R.pValues_cells = pValues_cells;
    return
end
CellCl_cand = zeros(NCl_cand, NCell);
lab_cellule = zeros(1, NCell);       % numéro du candidat de chaque cellule (0 = aucun)
for i = 1:NCl_cand
    CellCl_cand(i, assembly_candidates{i}) = 1;
    lab_cellule(assembly_candidates{i}) = i;
end

%% Étape 3 : activation des assemblées dans chaque SCE
RCl = zeros(NCl_cand, NRace);
for i = 1:NCl_cand
    RCl(i, :) = sum(Race(assembly_candidates{i}, :), 1);
end
Nrnd_SCE = sum(Race, 1);             % nb de cellules actives dans chaque SCE

switch null_cellules
    case 'uniforme'
        % Ancien modèle : cellules actives tirées parmi les NCell, calcul exact
        taille_cand = cellfun(@numel, assembly_candidates)';
        pValues_race = hygecdf(RCl - 1, NCell, repmat(taille_cand, 1, NRace), ...
            repmat(Nrnd_SCE, NCl_cand, 1), 'upper');
    case 'pondere'
        % Cellules actives tirées avec une probabilité proportionnelle à leur
        % activité : une cellule silencieuse ne peut pas être tirée
        activite = Nact';
        pValues_race = ones(NCl_cand, NRace);
        parfor j = 1:NRace
            colonne = ones(NCl_cand, 1);
            if Nrnd_SCE(j) > 0
                fixer_flux_aleatoire(graine + 1, j, NRace);
                RClr = tirages_ponderes(activite, Nrnd_SCE(j), NShuf, lab_cellule, NCl_cand);
                colonne = (sum(RClr >= RCl(:, j), 2) + 1) / (NShuf + 1);
            end
            pValues_race(:, j) = colonne;
        end
end

% 3.4 FDR, famille = les candidats testés dans un même SCE (comme avant)
PCl = false(NCl_cand, NRace);
for j = 1:NRace
    if Nrnd_SCE(j) == 0, continue; end
    PCl(:, j) = bh_fdr_mask(pValues_race(:, j)', FDR_Q_act)';
end
PCl = double(PCl);

%% Étape 4 : finalisation et tri (visualisation)
keepAsm = sum(PCl, 2) >= config_asm.min_activations_asm;
assemblyortho   = assembly_candidates(keepAsm);   % LISTE FINALE
cluster_origine = cand_cl(keepAsm);
PCl_final       = PCl(keepAsm, :);                % MATRICE EVENTS FINALE
CellCl_final    = CellCl_cand(keepAsm, :);        % MATRICE CELLULES FINALE
NCl             = sum(keepAsm);
if NCl == 0
    R = sortie_vide(NCell, NRace, 'aucune assemblée active dans assez de SCE');
    R.pValues_cells = pValues_cells;
    R.pValues_race  = pValues_race;
    return
end

% Tri des événements (colonnes)
Cl0 = find(sum(PCl_final, 1) == 0);   % aucune assemblée
Cl1 = find(sum(PCl_final, 1) == 1);   % 1 seule assemblée active
Cl2 = find(sum(PCl_final, 1) == 2);   % 2 assemblées
Cl3 = find(sum(PCl_final, 1) >= 3);   % 3+ assemblées
Bin = 2.^(0:NCl-1);
[~, x01] = sort(Bin * PCl_final(:, Cl1)); Cl1 = Cl1(x01);
[~, x02] = sort(Bin * PCl_final(:, Cl2)); Cl2 = Cl2(x02);
[~, x03] = sort(Bin * PCl_final(:, Cl3)); Cl3 = Cl3(x03);
RList = [Cl0, Cl1, Cl2, Cl3];
% Tri des cellules (lignes)
[~, x1] = sort(Bin * CellCl_final);

R = struct('assemblyortho', {assemblyortho}, 'cluster_origine', cluster_origine, ...
    'PCl_final', PCl_final, 'CellCl_final', CellCl_final, ...
    'pValues_cells', pValues_cells, 'pValues_race', pValues_race, ...
    'x1', x1, 'RList', RList, 'statut', 'ok');
end


function R = sortie_vide(NCell, NRace, raison)
R = struct('assemblyortho', {cell(1, 0)}, 'cluster_origine', zeros(1, 0), ...
    'PCl_final', zeros(0, NRace), 'CellCl_final', zeros(0, NCell), ...
    'pValues_cells', [], 'pValues_race', [], ...
    'x1', 1:NCell, 'RList', 1:NRace, 'statut', raison);
end
