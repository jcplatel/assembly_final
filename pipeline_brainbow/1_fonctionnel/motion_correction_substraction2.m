function [Tr1b, bad_frames] = motion_correction_substraction2(Tr1b, path, speed, regle_ancienne) %#ok<INUSL>
% Version 2 de motion_correction_substraction.
% Les trames où le recalage suite2p est mauvais (chute de ops.corrXY par
% rapport à sa tendance locale) sont repérées puis remplacées par
% interpolation linéaire.
% Changements :
%  - seuil à -3 écarts-types VRAIS. Avant, std() était calculé sur les seules
%    déviations négatives : pour un bruit gaussien, cela vaut 0,6 sigma, donc
%    le seuil était à -1,8 sigma (~3,6 % des trames, ~10 % après
%    élargissement à +/- 1 trame, interpolées et retirées du repos).
%    Ici sigma = racine de la moyenne des carrés des déviations négatives,
%    ce qui vaut bien sigma pour une demi-gaussienne.
%  - vérifie que corrXY a autant de trames que Tr1b (sinon : pas de correction)
%  - les lignes « bad frames sans mouvement » sont retirées : elles ne
%    servaient pas, et speed' pouvait créer une matrice Nz x Nz
% speed n'est plus utilisée (gardée pour ne pas changer l'appel).
% regle_ancienne = true (facultatif) : ancien seuil (std des seules déviations
% négatives, soit ~ -1,8 vrais sigma ; sert à retrouver_ancien).

if nargin < 4, regle_ancienne = false; end

bad_frames = [];
Nz = size(Tr1b, 2);

fallPath = strcat(path, 'Fall.mat');
if ~isfile(fallPath)
    return
end
S = load(fallPath, 'ops');
if ~isfield(S, 'ops') || ~isfield(S.ops, 'corrXY')
    return
end
corrXY = double(S.ops.corrXY(:)');
if numel(corrXY) ~= Nz
    warning('corrXY (%d trames) et F (%d trames) n''ont pas la même longueur : pas de correction de mouvement.', ...
        numel(corrXY), Nz);
    return
end

% Approche robuste : déviation par rapport à la tendance locale
rolling_median = movmedian(corrXY, 300);
deviation = corrXY - rolling_median;

% Bad frames = celles qui dévient fortement vers le bas
if regle_ancienne
    sigma_dev = std(deviation(deviation < 0));
else
    sigma_dev = sqrt(mean(deviation(deviation < 0).^2));
end
seuil_bad = -3 * sigma_dev;
bad_frames = deviation < seuil_bad;
bad_frames = conv(double(bad_frames), [1 1 1], 'same') > 0;   % trames voisines aussi

% Interpolation linéaire à la place des bad frames
Tr1b_clean = Tr1b;
Tr1b_clean(:, bad_frames) = NaN;
Tr1b = fillmissing(Tr1b_clean, 'linear', 2, 'EndValues', 'nearest');
end
