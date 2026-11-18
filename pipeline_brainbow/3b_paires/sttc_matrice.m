function S = sttc_matrice(R, dt)
%STTC_MATRICE  Spike Time Tiling Coefficient de toutes les paires
%   (Cutts & Eglen 2014, J Neurosci).
%   R  : cellules x trames, 1 = événement
%   dt : fenêtre de coïncidence en trames (+/- dt)
%   Pour la paire (A, B) :
%     TA = part du temps à moins de dt d'un événement de A
%     PA = part des événements de A à moins de dt d'un événement de B
%     STTC = 1/2 [ (PA - TB) / (1 - PA TB) + (PB - TA) / (1 - PB TA) ]
%   0 = indépendance, 1 = synchronie parfaite. NaN si une cellule n'a aucun
%   événement.

R = double(R > 0);
N = size(R, 1);
T = double(movmax(R, [dt dt], 2) > 0);     % temps « couvert » par chaque cellule
Tf = mean(T, 2);                            % TA de chaque cellule (colonne)
n = sum(R, 2);
P = (R * T') ./ n;                          % P(a, b) = PA pour la paire (a, b)
Ta = repmat(Tf, 1, N);                      % Ta(a, b) = TA
Tb = Ta';                                   % Tb(a, b) = TB
S = 0.5 * ((P - Tb) ./ (1 - P .* Tb) + (P' - Ta) ./ (1 - P' .* Ta));
S(n == 0, :) = NaN;
S(:, n == 0) = NaN;
S(logical(eye(N))) = NaN;
end
