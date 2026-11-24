function [stat, est_cellule, corriges] = charger_masques(dossier)
%% =========================================================================
%  CHARGER_MASQUES (18 nov.)
%  Masques des ROI pour les couleurs : masques_corriges.mat s'il existe dans
%  le dossier (corriger_masques, animaux 446178 et 446198), sinon Fall.mat.
%  est_cellule = iscell de Fall.mat (inchangé). corriges = true si les
%  masques corrigés sont utilisés.
%  =========================================================================
d = regexprep(char(dossier), '[\\/]+$', '');
F = load(fullfile(d, 'Fall.mat'), 'iscell');
est_cellule = F.iscell;
fc = fullfile(d, 'masques_corriges.mat');
corriges = isfile(fc);
if corriges
    M = load(fc, 'stat'); stat = M.stat;
else
    M = load(fullfile(d, 'Fall.mat'), 'stat'); stat = M.stat;
end
end
