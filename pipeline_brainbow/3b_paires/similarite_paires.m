function S = similarite_paires(Tr1b, Raster, WinRest, sampling_rate, opts_paires, Race)
%SIMILARITE_PAIRES  Similarité d'activité de toutes les paires de cellules,
%   AU REPOS seulement (trames WinRest : immobilité, bad frames exclues).
%   AUCUNE couleur n'entre ici.
%
%   S.sttc : STTC des événements détectés (Cutts & Eglen 2014), fenêtre de
%            coïncidence +/- dt. Événement = trame du PIC de chaque transient
%            détecté au repos (Raster). Indice de corrélation des trains
%            d'événements qui dépend peu du taux de chaque cellule.
%            Aussi S.sttc_100ms, S.sttc_500ms... (fenêtres de sensibilité).
%   S.dff  : corrélation de Pearson du ΔF/F filtré (Tr1b). Plus sensible :
%            deux cellules plus actives (ou avec un meilleur rapport
%            signal/bruit) sont plus corrélées entre elles, quelle que soit
%            leur lignée. Mesure secondaire.
%   S.sce  : co-participation aux SCE (coparticipation_sce : phi des lignes
%            de Race), seulement si Race est donné (8 nov.).
%   S.n_evenements, S.duree_repos_s : pour le diagnostic du taux par couleur.

repos = WinRest(:)';
R = double(Raster(:, repos) > 0);      % 1 = pic d'un transient détecté au repos
[mesures, dt_s] = liste_mesures_paires(opts_paires);
S.dt_trames = struct();
for k = find(~isnan(dt_s))
    dt = max(1, round(dt_s(k) * sampling_rate));        % fenêtre en trames (+/- dt)
    S.(mesures{k}) = sttc_matrice(R, dt);
    S.dt_trames.(mesures{k}) = dt;
end

X = double(Tr1b(:, repos));
X = X - mean(X, 2);
sd = std(X, 1, 2);
Z = X ./ sd;
C = (Z * Z') / size(X, 2);
C(sd == 0, :) = NaN;
C(:, sd == 0) = NaN;
C(logical(eye(size(C)))) = NaN;
S.dff = C;

if nargin >= 6 && ~isempty(Race)
    S.sce = coparticipation_sce(Race);
end

S.n_evenements  = sum(R, 2);
S.duree_repos_s = numel(repos) / sampling_rate;
end
