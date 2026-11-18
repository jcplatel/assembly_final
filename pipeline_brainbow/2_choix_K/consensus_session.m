function [statut, resume] = consensus_session(session_nom, dossier_etape1, dossier_sortie, liste_K, config, graine)
%CONSENSUS_SESSION  Règle de choix 'consensus' pour UNE session (étape 2).
%   Assemblées de consensus entre les K de liste_K (consensus_assemblees),
%   puis leurs activations dans les SCE, avec le même test que RACEORTHO
%   (activations_assemblees) et les mêmes filtres : plus de MIN_CELL_SIZE
%   cellules, au moins min_activations_asm activations significatives.
%   Écrit <dossier_sortie>\consensus\assemblees.mat, avec les mêmes noms de
%   variables qu'un k<K>\assemblees.mat (assemblyortho, PCl_final, x1,
%   RList, resume_K, config_asm, ...) : l'étape 3 et l'interface le lisent
%   comme un K (convention : K = 0 = consensus).
%   Aucune couleur n'est lue.

pc = config.choixK.consensus;
% Réglages des assemblées réellement utilisés à l'étape 1 (ceux du plus
% petit K calculé) : les activations sont testées comme à l'étape 1
asm = config.asm;
resume_K = struct();
for K = liste_K
    fk = fullfile(dossier_etape1, sprintf('k%d', K), 'assemblees.mat');
    if ~isfile(fk), continue; end
    noms = intersect({'resume_K', 'config_asm'}, who('-file', fk));
    R = load(fk, noms{:});
    if isfield(R, 'config_asm'), asm = R.config_asm; end
    if isfield(R, 'resume_K'), resume_K = R.resume_K; end
    break
end
P = load(fullfile(dossier_etape1, 'pretraitement.mat'), 'Race');
Race = double(P.Race > 0);
[NCell, NRace] = size(Race);
opts = struct('seuil', pc.seuil, 'min_K', pc.min_K, 'normalisation', pc.normalisation, ...
    'taille_min', asm.MIN_CELL_SIZE + 1, 'dossier_analyse', dossier_sortie);   % v2 : assemblées de chaque K en version 2
C = consensus_assemblees(dossier_etape1, liste_K, NCell, opts);

dossier = fullfile(dossier_sortie, 'consensus');
fichier = fullfile(dossier, 'assemblees.mat');
if isfile(fichier), delete(fichier); end
resume = struct('n_asm', 0, 'tailles', '', 'K_presence', '', 'n_K', numel(C.K_utilises));
if numel(C.K_utilises) < 2
    statut = sprintf('exclue : %d K calculés à l''étape 1 (il en faut au moins 2 pour un consensus)', ...
        numel(C.K_utilises));
    return
end

% Activations : même test qu'à l'étape 1 (graine propre à la session)
[PCl, ~] = activations_assemblees(Race, C.assemblees, asm, graine + 7);
garde = sum(PCl, 2)' >= asm.min_activations_asm;
assemblyortho   = C.assemblees(garde);
PCl_final       = PCl(garde, :);
consensus_moyen = C.consensus_moyen(garde);
K_presence      = C.K_presence(garde);
NCl = numel(assemblyortho);
if NCl == 0
    statut = 'exclue : aucune assemblée de consensus';
    return
end

% Mêmes sorties qu'un K (tri des cellules et des SCE pour l'affichage)
CellCl_final = zeros(NCl, NCell);
for i = 1:NCl
    CellCl_final(i, assemblyortho{i}) = 1;
end
Cl0 = find(sum(PCl_final, 1) == 0);
Cl1 = find(sum(PCl_final, 1) == 1);
Cl2 = find(sum(PCl_final, 1) == 2);
Cl3 = find(sum(PCl_final, 1) >= 3);
Bin = 2.^(0:NCl-1);
[~, x01] = sort(Bin * PCl_final(:, Cl1)); Cl1 = Cl1(x01);
[~, x02] = sort(Bin * PCl_final(:, Cl2)); Cl2 = Cl2(x02);
[~, x03] = sort(Bin * PCl_final(:, Cl3)); Cl3 = Cl3(x03);
RList = [Cl0, Cl1, Cl2, Cl3];
[~, x1] = sort(Bin * CellCl_final);
cluster_origine = 1:NCl;
NClOK = NCl;
NRaceOK = sum(any(PCl_final > 0, 1));
K_utilises = C.K_utilises;
opts_consensus = opts;
config_asm = asm;

% Ligne de résumé : celle du plus petit K (champs de la session), adaptée
tailles = cellfun(@numel, assemblyortho);
if ~isempty(fieldnames(resume_K))
    resume_K.NCl_ini = 0;                                  % 0 = consensus
    if isfield(resume_K, 'SCE') && resume_K.SCE > 0
        resume_K.SCEok_freq = resume_K.sce_freq * NRaceOK / resume_K.SCE;
    end
    resume_K.SCEok = NRaceOK;
    resume_K.Assemblies = NCl;
    resume_K.SCE_Purity = NaN;
    resume_K.cell_Purity = NaN;
    resume_K.sil_all_clusters = NaN;
    resume_K.SClok = NaN;
    resume_K.SClokfinal = NaN;
    resume_K.cell_assembly = num2str(tailles);
    resume_K.mean_cell_assembly = mean(tailles);
end

if ~isfolder(dossier), mkdir(dossier); end
save(fichier, 'assemblyortho', 'PCl_final', 'CellCl_final', 'x1', 'RList', 'cluster_origine', ...
    'NClOK', 'NRace', 'NRaceOK', 'consensus_moyen', 'K_presence', 'K_utilises', 'opts_consensus', ...
    'resume_K', 'config_asm');
statut = 'ok';
resume.n_asm = NCl;
resume.tailles = strjoin(arrayfun(@(t) sprintf('%d', t), tailles, 'UniformOutput', false), ' ');
resume.K_presence = strjoin(cellfun(@texte_K, K_presence, 'UniformOutput', false), ' | ');
end


function t = texte_K(Kp)
if isempty(Kp)
    t = '-';
elseif numel(Kp) == 1
    t = sprintf('K %d', Kp);
else
    t = sprintf('K %d-%d (%d)', min(Kp), max(Kp), numel(Kp));
end
end
