function [asm, best_K] = assemblees_retenues(session, NCell)
%ASSEMBLEES_RETENUES  Assemblées de la session au K retenu (étapes 2-3) :
%   asm = cellules de chaque assemblée (dans l'ordre de assemblyortho, le
%   même numéro que dans les tableaux de l'étape 3) ; {} si rien.
%   best_K = 0 : assemblées de consensus.
asm = {};
best_K = NaN;
f = fullfile(session.dossier, 'choix_K.mat');
if ~isfile(f), return; end
choix = load(f);
if ~isfield(choix, 'statut_choix') || ~strcmp(choix.statut_choix, 'ok'), return; end
best_K = choix.best_K;
fa = fichier_assemblees(session.dossier, session.dossier_source, best_K);   % v1 ou v2 selon l'analyse
if ~isfile(fa), return; end
S = load(fa, 'assemblyortho');
asm = cell(1, numel(S.assemblyortho));
for a = 1:numel(S.assemblyortho)
    c = unique(S.assemblyortho{a}(:))';
    asm{a} = c(c >= 1 & c <= NCell);
end
end
