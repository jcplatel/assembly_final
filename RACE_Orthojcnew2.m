% =================================================================
% SCRIPT : RACEORTHOJC_NEW (Version Robustesse / Données Bruitées)
% =================================================================
% Objectif : Maximiser la détection d'assemblées (Recall) en relaxant
% les seuils statistiques, puis filtrer par cohérence biologique.

% --- Paramètres ---
[NCell, NRace] = size(Race);
NCl_input = NClOK;        % Nombre de clusters suggéré par K-means initial
NShuf = 5000;             % Nombre de permutations
P_VAL_THRESHOLD = 0.05;   % Seuil relaxé (pas de correction Bonferroni ici)
MIN_CELL_SIZE = 5;        % Nombre min de cellules pour former une assemblée

% fprintf('--- DÉBUT RACEORTHOJC_NEW ---\n');
% fprintf('Paramètres : %d Shuffles, P-val < %.2f, Taille Min > %d\n', ...
%     NShuf, P_VAL_THRESHOLD, MIN_CELL_SIZE);

%% Étape 1 : Calcul de l'association Cellule-Cluster (Permutations)

% 1.1 Calculs préliminaires (Observed Data)
CellP = zeros(NCell, NCl_input); % Participation count
CellR = zeros(NCell, NCl_input); % Participation rate
for i = 1:NCl_input
    CellP(:, i) = sum(Race(:, IDX2 == i), 2);
    mask_size = sum(IDX2 == i);
    if mask_size > 0
        CellR(:, i) = CellP(:, i) / mask_size;
    end
end

% 1.2 Préparation Parallel Pool
cRace   = parallel.pool.Constant(Race);
maskCl = cell(1, NCl_input);
for i = 1:NCl_input, maskCl{i} = (IDX2 == i); end
cMaskCl = parallel.pool.Constant(maskCl);
cCellP  = parallel.pool.Constant(CellP);

% 1.3 Boucle de permutation (Cellules)
pValues_cells = ones(NCell, NCl_input);

parfor j = 1:NCell
    RaceLoc   = cRace.Value;
    maskClLoc = cMaskCl.Value;
    CellPLoc  = cCellP.Value;
    Nrnd = sum(RaceLoc(j, :)); % Nb d'événements pour cette cellule
    
    if Nrnd == 0
        continue; 
    end
    
    % Shuffles
    RClr = zeros(NCl_input, NShuf, 'uint16');
    for l = 1:NShuf
        randIdx = randperm(size(RaceLoc, 2), Nrnd);
        for i = 1:NCl_input
            RClr(i, l) = sum(maskClLoc{i}(randIdx));
        end
    end
    
    % Calcul P-value (Méthode exacte : (b + 1) / (m + 1))
    observed_counts = CellPLoc(j, :)';
    pValues_cells(j, :) = (sum(RClr >= observed_counts, 2) + 1)' ./ (NShuf + 1);
end

% 1.4 Filtrage Statistique (Relaxé)
CellCl = uint16((pValues_cells <= P_VAL_THRESHOLD)'); % [NCl_input x NCell]

%% Étape 2 : Nettoyage et Formation des Candidats

% 2.1 Gestion des cellules "Promiscuous" (Multi-clusters)
A2 = find(sum(CellCl, 1) >= 2); 
for i = A2
    [~, idx] = max(CellR(i, :)); % On assigne au cluster où le taux est max
    CellCl(:, i) = 0;
    CellCl(idx, i) = 1;
end

% 2.2 Filtrage par taille d'assemblée (> 5 cellules)
C0 = cell(0);
k = 0;
for i = 1:NCl_input
    cell_indices = find(CellCl(i, :));
    if length(cell_indices) > MIN_CELL_SIZE
        k = k + 1;
        C0{k} = cell_indices;
    end
end
assembly_candidates = C0;

% --- CORRECTION CRITIQUE ---
% On recrée une matrice binaire propre alignée sur les candidats retenus.
% Cela évite les erreurs de dimensions plus tard.
NCl_cand = length(assembly_candidates);
CellCl_cand = zeros(NCl_cand, NCell);
for i = 1:NCl_cand
    CellCl_cand(i, assembly_candidates{i}) = 1;
end

% fprintf('Candidats retenus (taille > %d) : %d assemblées.\n', MIN_CELL_SIZE, NCl_cand);

if NCl_cand == 0
    % fprintf('Aucune assemblée candidate ne respecte les critères. Fin.\n');
    assemblyortho = cell(0); NCl = 0; return;
end

%% Étape 3 : Test d'Activation des Assemblées (RACE Events)

% 3.1 Somme observée des spikes par assemblée et par event RACE
RCl = zeros(NCl_cand, NRace);
for i = 1:NCl_cand
    RCl(i, :) = sum(Race(assembly_candidates{i}, :));
end

% 3.2 Préparation Parallel Pool (Candidats)
cC0   = parallel.pool.Constant(assembly_candidates);
cRCl  = parallel.pool.Constant(RCl);
cNrnd = parallel.pool.Constant(sum(Race, 1)); % Activité totale par RACE

% 3.3 Boucle de permutation (Events)
pValues_race = ones(NCl_cand, NRace);

parfor j = 1:NRace
    C0loc   = cC0.Value;
    RCl_loc = cRCl.Value;
    Nrnd    = cNrnd.Value(j);
    NCellLoc = size(cRace.Value, 1);
    
    if Nrnd == 0
        continue;
    end
    
    % Shuffles
    RClr = zeros(NCl_cand, NShuf, 'uint16');
    for l = 1:NShuf
        idx = randperm(NCellLoc, Nrnd); % Indices aléatoires de cellules actives
        for i = 1:NCl_cand
            % Combien tombent dans chaque assemblée ?
            RClr(i, l) = sum(ismember(idx, C0loc{i}));
        end
    end
    
    % Calcul P-value
    observed_counts = RCl_loc(:, j);
    pValues_race(:, j) = (sum(RClr >= observed_counts, 2) + 1) ./ (NShuf + 1);
end

% 3.4 Matrice binaire d'activation significative
PCl = double(pValues_race <= P_VAL_THRESHOLD);

%% Étape 4 : Finalisation et Tri (Visualisation)

% 4.1 Filtrage final : On ne garde que ce qui s'active au moins 1 fois
keepAsm = sum(PCl, 2) > 0;

assemblyortho = assembly_candidates(keepAsm); % LISTE FINALE
PCl_final     = PCl(keepAsm, :);              % MATRICE EVENTS FINALE
CellCl_final  = CellCl_cand(keepAsm, :);      % MATRICE CELLULES FINALE
NCl           = sum(keepAsm);                 % NB FINAL

% fprintf('Assemblées finales actives : %d\n', NCl);
% 
% if NCl == 0
%     fprintf('Aucune assemblée active trouvée après test RACE. Fin.\n');
%     return;
% end

% 4.2 Tri pour Raster Plot (Logique binaire)
% Note : Si NCl > 52, la précision double est atteinte, le tri peut être approximatif.
% Mais pour NCl < 50, c'est parfait.

% Tri des événements (Colonnes - RACEs)
Cl0 = find(sum(PCl_final, 1) == 0); % Bruit / Pas d'assemblée
Cl1 = find(sum(PCl_final, 1) == 1); % 1 seule assemblée active
Cl2 = find(sum(PCl_final, 1) == 2); % 2 assemblées
Cl3 = find(sum(PCl_final, 1) >= 3); % 3+ assemblées (Riche/Complexe)

Bin = 2.^(0:NCl-1); % Vecteur de poids binaires

% Tri interne des groupes
[~, x01] = sort(Bin * PCl_final(:, Cl1)); Cl1 = Cl1(x01);
[~, x02] = sort(Bin * PCl_final(:, Cl2)); Cl2 = Cl2(x02);
[~, x03] = sort(Bin * PCl_final(:, Cl3)); Cl3 = Cl3(x03);

% Liste finale des indices d'événements triés
RList = [Cl0, Cl1, Cl2, Cl3]; 

% Tri des cellules (Lignes) pour l'axe Y du raster
[~, x1] = sort(Bin * CellCl_final);

% fprintf('--- FIN DU TRAITEMENT ---\n');

% Sauvegarde optionnelle
% save('Resultats_RaceOrtho_New.mat', 'assemblyortho', 'PCl_final', 'RList', 'x1', 'NCl');
