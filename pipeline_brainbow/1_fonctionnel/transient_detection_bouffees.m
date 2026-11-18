function [Raster, Acttmp2, th_detection, PromALL, WidthALL, info] = ...
    transient_detection_bouffees(Tr1b, WinRest, DFF0, sampling_rate, opts)
%TRANSIENT_DETECTION_BOUFFEES  (18 nov. 2026, demandée par l'utilisateur)
%   Détection des transients avec exclusion des bouffées, d'après le code
%   d'Arnaud Malvache (preprocessing_Arnaud.m). Utilisée seulement si
%   opts.detection vaut 'arnaud' ou 'bouffees' (sinon preprocessing_7 garde
%   transient_detection_shot_theta2, rien ne change).
%
%   Bouffée : instant où la médiane locale de Tr1b (fenêtre +/- opts.bouffee_demi_fenetre_s,
%   défaut 4 s = 40 trames à 10 Hz comme le code d'Arnaud) dépasse
%   median(Tr1b) + iqr(Tr1b)/2 (toute la trace de la cellule). Compatible avec
%   la ligne de base glissante : F0 (10e percentile sur 120 s) ne suit pas
%   une bouffée de quelques secondes, qui reste visible dans DFF0 / Tr1b.
%
%   'arnaud'   : comme preprocessing_Arnaud, au repos (WinRest), hors
%                bouffée : transient à la trame j si Tr1b(j) - médiane locale
%                > 3 x IQR local, et aucun transient dans la seconde
%                précédente (opts.refractaire_s, défaut 1 s = 10 trames).
%   'bouffees' : notre détection (transient_detection_shot_theta2,
%                threshold_peak x bruit) avec MinPeakDistance = 1 s, puis
%                pics situés dans une bouffée retirés.
%   info : nombre de transients retirés par la règle des bouffées, % du repos en bouffée.

mode = lower(char(opts.detection));
if isfield(opts, 'bouffee_demi_fenetre_s'), dfs = opts.bouffee_demi_fenetre_s; else, dfs = 4; end
if isfield(opts, 'refractaire_s'), rs = opts.refractaire_s; else, rs = 1; end
W = max(1, round(dfs * sampling_rate));
refr = max(1, round(rs * sampling_rate));
[NCell, Nz] = size(Tr1b);
repos = false(1, Nz); repos(WinRest) = true;
Tr1b = double(Tr1b);

% bouffées (pour les deux modes)
medloc = movmedian(Tr1b, [W W], 2);
seuil_b = median(Tr1b, 2) + iqr(Tr1b, 2) / 2;
bouffee = medloc >= seuil_b;
info = struct('mode', mode, 'demi_fenetre_trames', W, 'refractaire_trames', refr, ...
    'pct_repos_en_bouffee', 100 * mean(bouffee(:, repos), 2), 'n_retires_bouffee', zeros(NCell, 1));

switch mode
    case 'bouffees'
        [Raster, Acttmp2, th_detection, PromALL, WidthALL] = ...
            transient_detection_shot_theta2(Tr1b, refr, opts.threshold_peak, WinRest, DFF0);
        for i = 1:NCell
            l = Acttmp2{i};
            if isempty(l), continue; end
            dedans = bouffee(i, l);
            info.n_retires_bouffee(i) = sum(dedans);
            Raster(i, l(dedans)) = 0;
            Acttmp2{i} = l(~dedans); PromALL{i} = PromALL{i}(~dedans); WidthALL{i} = WidthALL{i}(~dedans);
        end

    case 'arnaud'
        Raster = zeros(NCell, Nz);
        Acttmp2 = cell(1, NCell); PromALL = cell(1, NCell); WidthALL = cell(1, NCell);
        th_detection = nan(1, NCell);
        n_ret = zeros(NCell, 1);
        parfor i = 1:NCell
            tr = Tr1b(i, :); md = medloc(i, :); bf = bouffee(i, :);
            act = false(1, Nz); amp = nan(1, Nz); seuils = nan(1, Nz); nr = 0;
            dernier = -Inf;
            for j = W + 1:Nz - W
                if ~repos(j) || j - dernier <= refr, continue; end
                fen = sort(tr(j - W:j + W));
                n = numel(fen);
                q = iqr_trie(fen, n);
                if tr(j) - md(j) > 3 * q
                    if bf(j), nr = nr + 1; continue; end
                    act(j) = true; amp(j) = tr(j) - md(j); seuils(j) = 3 * q; dernier = j;
                end
            end
            l = find(act);
            Acttmp2{i} = l; PromALL{i} = amp(l); WidthALL{i} = nan(size(l));
            Raster(i, :) = act;
            th_detection(i) = median(seuils(l), 'omitnan');
            n_ret(i) = nr;
        end
        info.n_retires_bouffee = n_ret;

    otherwise
        error('opts.detection doit valoir ''arnaud'' ou ''bouffees''.');
end
end


function q = iqr_trie(x, n)
% IQR d'un vecteur déjà trié (même définition que iqr de MATLAB : quantiles
% avec interpolation, positions (k - 0,5) / n)
q = quantile_trie(x, n, 0.75) - quantile_trie(x, n, 0.25);
end

function v = quantile_trie(x, n, p)
pos = p * n + 0.5;
if pos <= 1, v = x(1); return; end
if pos >= n, v = x(n); return; end
k = floor(pos); f = pos - k;
v = x(k) + f * (x(k + 1) - x(k));
end
