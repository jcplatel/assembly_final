function [PCl, pValues_race] = activations_assemblees(Race, assemblees, config_asm, graine)
%ACTIVATIONS_ASSEMBLEES  SCE où chaque assemblée est significativement active.
%   Même calcul que l'étape 3 de RACEORTHO_FDR2 (copié tel quel), pour des
%   assemblées venues d'ailleurs (consensus entre les K) :
%     - nombre de cellules de l'assemblée actives dans chaque SCE
%     - modèle nul config_asm.null_cellules ('pondere' : cellules tirées selon
%       leur activité ; 'uniforme' : loi hypergéométrique)
%     - FDR de Benjamini-Hochberg dans chaque SCE (FDR_Q_activations, sinon FDR_Q)
%   assemblees : cellules de chaque assemblée, sans cellule commune à deux
%   assemblées. PCl : assemblées x SCE (1 = activation significative).

[NCell, NRace] = size(Race);
nA = numel(assemblees);
PCl = zeros(nA, NRace);
pValues_race = ones(nA, NRace);
if nA == 0, return; end
FDR_Q_act = config_asm.FDR_Q;
if isfield(config_asm, 'FDR_Q_activations') && ~isempty(config_asm.FDR_Q_activations)
    FDR_Q_act = config_asm.FDR_Q_activations;
end
lab_cellule = zeros(1, NCell);       % numéro de l'assemblée de chaque cellule (0 = aucune)
for i = 1:nA
    lab_cellule(assemblees{i}) = i;
end
RCl = zeros(nA, NRace);
for i = 1:nA
    RCl(i, :) = sum(Race(assemblees{i}, :), 1);
end
Nrnd_SCE = sum(Race, 1);             % nb de cellules actives dans chaque SCE
Nact = sum(Race, 2);

switch config_asm.null_cellules
    case 'uniforme'
        taille = cellfun(@numel, assemblees(:));
        pValues_race = hygecdf(RCl - 1, NCell, repmat(taille, 1, NRace), ...
            repmat(Nrnd_SCE, nA, 1), 'upper');
    case 'pondere'
        activite = Nact';
        NShuf = config_asm.NShuf;
        parfor j = 1:NRace
            colonne = ones(nA, 1);
            if Nrnd_SCE(j) > 0
                fixer_flux_aleatoire(graine + 1, j, NRace);
                RClr = tirages_ponderes(activite, Nrnd_SCE(j), NShuf, lab_cellule, nA);
                colonne = (sum(RClr >= RCl(:, j), 2) + 1) / (NShuf + 1);
            end
            pValues_race(:, j) = colonne;
        end
    otherwise
        error('null_cellules doit valoir ''pondere'' ou ''uniforme''.');
end

for j = 1:NRace
    if Nrnd_SCE(j) == 0, continue; end
    PCl(:, j) = bh_fdr_mask(pValues_race(:, j)', FDR_Q_act)';
end
end
