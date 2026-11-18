function [Raster, Acttmp2, th_detection, PromALL, WidthALL] = ...
    transient_detection_shot_theta2(Tr1b, MinPeakDistance, threshold_peak, WinRest, DFF0, bruit_ancien)
% Version 2 de transient_detection_shot_theta.
% Changements :
%  - le bruit est estimé UNIQUEMENT à l'intérieur des épisodes de repos.
%    Avant, diff(DFF0(i, WinRest)) mettait bout à bout des épisodes séparés,
%    ce qui créait un saut artificiel à chaque raccord.
%  - BurstMask, ampALL et fenetre_correction, qui ne servaient pas, sont retirés.
% À SAVOIR : le seuil est en unités de bruit de DFF0 (brut), et il est
% appliqué à la proéminence de Tr1b (filtré). threshold_peak = 3.09 n'est
% donc PAS un z-score à p = 0.001 : c'est un seuil relatif au bruit de
% chaque cellule, à régler en regardant des traces, sans les couleurs.
% bruit_ancien = true (facultatif) : bruit calculé comme avant, épisodes de
% repos mis bout à bout (sert à retrouver_ancien).

if nargin < 6, bruit_ancien = false; end

[NCell, Nz] = size(Tr1b);
Raster   = zeros(NCell, Nz);
Acttmp2  = cell(1, NCell);
PromALL  = cell(1, NCell);
WidthALL = cell(1, NCell);
th_detection = nan(1, NCell);

repos = false(1, Nz);
repos(WinRest) = true;
diff_dans_repos = repos(1:end-1) & repos(2:end);   % trames t et t+1 toutes deux au repos

for i = 1:NCell
    if bruit_ancien
        noise = 1.4826 * mad(diff(double(DFF0(i, WinRest))), 1) / sqrt(2);
    else
        d = diff(double(DFF0(i, :)));
        noise = 1.4826 * mad(d(diff_dans_repos), 1) / sqrt(2);
    end
    if ~(noise > 0)        % pas de repos, ou trace plate : cellule ignorée
        continue
    end
    th_detection(i) = threshold_peak * noise;

    [~, locs, w, p] = findpeaks(double(Tr1b(i, :)), ...
        'MinPeakProminence', th_detection(i), ...
        'MinPeakDistance', MinPeakDistance);

    peaks_Rest = repos(locs);          % seulement les pics au repos
    Acttmp2{i}  = locs(peaks_Rest);
    PromALL{i}  = p(peaks_Rest);       % proéminence nette (pour l'analyse d'amplitude)
    WidthALL{i} = w(peaks_Rest);
    Raster(i, Acttmp2{i}) = 1;
end
end
