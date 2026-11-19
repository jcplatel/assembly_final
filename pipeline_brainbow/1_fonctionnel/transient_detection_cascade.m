function [Raster, Acttmp2, th_detection, PromALL, WidthALL, info] = ...
    transient_detection_cascade(DFF0, WinRest, path, sampling_rate, opts)
%TRANSIENT_DETECTION_CASCADE  (18 nov. 2026, demandée par l'utilisateur)
%   Détection des événements à partir de la déconvolution CASCADE
%   (cascade_spikes, modèle GCaMP8 à 10 Hz). Utilisée seulement si
%   opts.detection = 'cascade' (sinon rien ne change).
%
%   spike_prob = nombre de potentiels d'action estimé par image.
%   Événement : suite d'images où spike_prob > opts.cascade_seuil_image
%   (défaut 0,1) ; sa taille = somme de spike_prob sur l'événement (nombre
%   de potentiels estimé), placé à l'image du maximum. Gardé si cette image
%   est au repos (WinRest) et si la taille >= opts.cascade_min_spikes :
%     ~0,5-1 = tous les potentiels isolés, 2 ou 3 = bouffées seulement.
%   PromALL = taille de l'événement, WidthALL = durée (images).
%   Cache : opts.cascade_cache (fichier), sinon <path>\cascade_<modele>.mat.
%   info : spike_prob n'est pas rendu (trop gros) ; nombre d'événements
%   par taille, pour choisir le seuil.

if isfield(opts, 'cascade_seuil_image'), s_img = opts.cascade_seuil_image; else, s_img = 0.1; end
if isfield(opts, 'cascade_min_spikes'), s_min = opts.cascade_min_spikes; else, s_min = 1; end
if isfield(opts, 'cascade_modele'), modele = opts.cascade_modele;
else, modele = 'GC8_EXC_10Hz_smoothing150ms_high_noise'; end
if isfield(opts, 'cascade_cache'), fc = opts.cascade_cache;
else, fc = fullfile(path, ['cascade_' modele '.mat']); end
o = opts; o.cascade_modele = modele;
sp = cascade_spikes(DFF0, o, fc);
sp(isnan(sp)) = 0;

[NCell, Nz] = size(sp);
repos = false(1, Nz); repos(WinRest) = true;
Raster = zeros(NCell, Nz);
Acttmp2 = cell(1, NCell); PromALL = cell(1, NCell); WidthALL = cell(1, NCell);
th_detection = s_min * ones(1, NCell);
bornes = [0 0.5 1 1.5 2 3 5 Inf];
compte = zeros(NCell, numel(bornes) - 1);
for i = 1:NCell
    au_dessus = sp(i, :) > s_img;
    d = diff([false au_dessus false]);
    debut = find(d == 1); fin = find(d == -1) - 1;
    taille = zeros(size(debut)); pic = zeros(size(debut));
    for e = 1:numel(debut)
        seg = sp(i, debut(e):fin(e));
        taille(e) = sum(seg);
        [~, k] = max(seg); pic(e) = debut(e) + k - 1;
    end
    au_repos = repos(pic);
    compte(i, :) = histcounts(taille(au_repos), bornes);
    garde = au_repos & taille >= s_min;
    l = pic(garde);
    Raster(i, l) = 1;
    Acttmp2{i} = l; PromALL{i} = taille(garde); WidthALL{i} = fin(garde) - debut(garde) + 1;
end
info = struct('modele', modele, 'seuil_image', s_img, 'min_spikes', s_min, ...
    'bornes_taille', bornes, 'n_evenements_repos_par_taille', compte, ...
    'spikes_par_min_repos', sum(sp(:, repos), 2) / max(eps, numel(WinRest) / sampling_rate / 60));
end
