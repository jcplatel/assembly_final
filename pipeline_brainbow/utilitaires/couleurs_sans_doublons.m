function colorcell = couleurs_sans_doublons(colorcell, session)
%COULEURS_SANS_DOUBLONS  Les ROI retirées à l'étape 1 sont traitées comme
%   non colorées (8 = none) : doublons (detecter_doublons : même cellule
%   comptée deux fois ; sinon sa couleur compterait deux fois) et ROI trop
%   petites ou trop grandes (tailles_roi). Sans effet si l'étape 1 a été
%   faite sans config.pre.retirer_doublons / retirer_tailles.
%   Aussi (15 nov.) les cellules écartées par config.pre.population
%   ('sequence' / 'hors_sequence') : les couleurs sont alors comptées et
%   mélangées parmi les cellules gardées seulement.
dossier = session.dossier;
if isfield(session, 'dossier_source') && ~isempty(session.dossier_source)
    dossier = session.dossier_source;
end
f = fullfile(dossier, 'pretraitement.mat');
if ~isfile(f) || ~ismember('diagnostic', who('-file', f)), return; end
S = load(f, 'diagnostic');
r = zeros(1, 0);
if isfield(S.diagnostic, 'doublons') && isfield(S.diagnostic.doublons, 'retirees')
    r = [r, S.diagnostic.doublons.retirees(:)'];
end
if isfield(S.diagnostic, 'tailles_roi') && isfield(S.diagnostic.tailles_roi, 'retirees')
    r = [r, S.diagnostic.tailles_roi.retirees(:)'];
end
% population (15 nov.) : cellules écartées par config.pre.population
if isfield(S.diagnostic, 'population') && isfield(S.diagnostic.population, 'retirees')
    r = [r, S.diagnostic.population.retirees(:)'];
end
r = r(r >= 1 & r <= numel(colorcell));
colorcell(r) = 8;
end
