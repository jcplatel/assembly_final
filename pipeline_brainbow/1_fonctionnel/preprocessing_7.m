function [Tr1b,speedsm,Raster,SumAct,MAct,Race,RasterRace,WinRest,WinActive,TRace,Fzero,DFF0,th_detection,...
    bad_frames,max_cells_allowed,n_transients_total,Acttmp2,PromALL,WidthALL,opts,diagnostic] = ...
    preprocessing_7(F,opts,speed,path,sampling_rate,namefull,graine)
% Version 7 de preprocessing_6. Mêmes entrées et sorties (+ graine et
% diagnostic). Changements :
%  - motion_correction_substraction2 (seuil à -3 vrais écarts-types)
%  - speed est coupée ou complétée (NaN) pour avoir la longueur de F
%  - bad_frames accepte un masque logique OU une liste de trames
%  - les bad frames sont élargies de (SG_window-1)/2 trames de chaque côté
%    (le filtre Savitzky-Golay étale l'artefact) et retirées du repos
%  - transient_detection_shot_theta2 : bruit estimé à l'intérieur des
%    épisodes de repos seulement
%  - SCE_detection5 : seuil sans plafond, SCE attendus par hasard
%  - Acttmp2 / PromALL / WidthALL mis à jour quand des trames sont mises à 0
%  - opts.retirer_hyperactives (false par défaut = comme avant) ; seuil
%    réglable (21 nov.) : opts.hyperactives_seuil_hz (défaut 0,2) et
%    opts.hyperactives_min (défaut 200 transients)
%  - opts.retirer_tailles (7 nov.) : ROI trop petites ou trop grandes
%    (tailles_roi) : raster mis à zéro
%  - opts.retirer_doublons (7 nov.) : ROI qui sont la même cellule comptée
%    deux fois (detecter_doublons) : raster de la plus petite mis à zéro
%  - opts.population (15 nov., facultatif) : 'sequence' ou 'hors_sequence' :
%    les SCE (et leur seuil) sont cherchés parmi les cellules de la séquence
%    de course seulement, ou parmi les autres (cellules_sequence) ; le
%    raster des autres cellules est mis à zéro. Champ absent = toutes.
% namefull n'est plus utilisé (gardé pour ne pas changer l'appel).

if nargin < 7, graine = 1; end
% opts.pretraitement_ancien = true (facultatif : réglage 'ancien_rang' et
% retrouver_ancien) : bad frames, repos et bruit calculés comme dans
% preprocessing_6
ancien = isfield(opts, 'pretraitement_ancien') && opts.pretraitement_ancien;

%% Load settings
MinPeakDistance   = opts.MinPeakDistance;          % 3 default
threshold_peak    = opts.threshold_peak;           % 3 default
SG_window         = opts.SG_window;                % default 7
motion_correction = opts.motion_correction;        % false
colorsubstraction = opts.colorsubstraction;        % false

%% import Fluo
Tr1b = single(F);
[NCell, Nz] = size(Tr1b);

%% correct for motion artifact
if motion_correction && ~isfile(strcat(path, 'Fall.mat'))
    warning(['Fall.mat introuvable dans %s : pas de correction de mouvement. ' ...
        'Vérifier dossier_suite2p dans la table (l''étape 3 en aura besoin aussi).'], path);
    motion_correction = false;
end
if motion_correction
    [Tr1b, bad_frames] = motion_correction_substraction2(Tr1b, path, speed, ancien);
    save(strcat(path, 'badframes.mat'), "bad_frames")
else
    bad_frames = false(1, Nz);
end
% bad_frames : masque logique de longueur Nz, quel que soit le format reçu
if numel(bad_frames) == Nz
    bad_frames = logical(bad_frames(:)');
else
    masque = false(1, Nz);
    masque(bad_frames) = true;
    bad_frames = masque;
end
% 23 nov. (demandé par l'utilisateur) : images rejetées à la main, ex. un
% artefact du champ (444119 220919 : une zone double de luminosité à 357 s).
% Fichier facultatif images_rejetees_manuelles.mat dans le dossier suite2p :
% intervalles_s (N x 2, en secondes depuis la 1re image) et raison (texte).
% Ces images sont traitées comme les bad frames : hors ligne de base, hors
% repos, raster à 0. La fenêtre de la ligne de base fait 120 s : pour une
% marche de fluorescence, rejeter au moins 60 s après la marche, sinon la
% ligne de base mélange les niveaux d'avant et d'après.
n_rejet_manuel = 0;
f_rejet = fullfile(path, 'images_rejetees_manuelles.mat');
if isfile(f_rejet)
    M = load(f_rejet, 'intervalles_s');
    t_img = (0:Nz-1) / sampling_rate;
    for k = 1:size(M.intervalles_s, 1)
        rej = t_img >= M.intervalles_s(k, 1) & t_img <= M.intervalles_s(k, 2);
        n_rejet_manuel = n_rejet_manuel + sum(rej & ~bad_frames);
        bad_frames = bad_frames | rej;
    end
    fprintf('Images rejetées à la main (%s) : %d\n', f_rejet, n_rejet_manuel);
end

%% correction de la fuite de couleur (optionnelle)
% Elle doit corriger la fluorescence à partir des CANAUX de couleur, jamais
% à partir des classes colorcell (sinon la couleur influencerait la détection).
if colorsubstraction
    [Tr1b, ~] = color_substraction(Tr1b, path);
end

%% vitesse : même longueur que la fluorescence
speed = double(speed(:)');
if numel(speed) >= Nz
    speed = speed(1:Nz);
else
    warning('speed (%d trames) plus court que F (%d) : fin complétée par NaN.', numel(speed), Nz);
    speed(end+1:Nz) = NaN;
end
speedsm = smoothdata(speed, 'gaussian', 50);

%% baseline calculation
[DFF0, Fzero] = baseline_calculation_glissante2(Tr1b, bad_frames, sampling_rate);

%% filtering
Tr1b = sgolayfilt(DFF0', 3, SG_window)';
% Le filtre étale un artefact sur (SG_window-1)/2 trames de chaque côté
demi_fenetre = (SG_window - 1) / 2;
bad_frames_etendues = logical(movmax(double(bad_frames), [demi_fenetre demi_fenetre]));
if ancien      % ancien prétraitement : bad frames ni élargies, ni retirées du repos
    bad_frames_etendues = logical(bad_frames);
end

%% définition fenetre active rest
isActive = (speedsm > 1);
isActive_expanded = logical(movmax(double(isActive), [3 3]));
WinActive = find(isActive_expanded);
% repos : immobile, vitesse pas trop négative, hors bad frames
WinRest_bin = ~isActive_expanded & speedsm > -0.1 & ~bad_frames_etendues;
if ancien
    WinRest_bin = ~isActive_expanded & speedsm > -0.1;
end
%% Rest period has to be 3sec
taille_minimum = floor(sampling_rate * 3);
WinRest_bin = bwareaopen(WinRest_bin, taille_minimum);
WinRest = find(WinRest_bin);

%% DÉTECTION DES TRANSIENTS
info_bouffees = [];
if isfield(opts, 'detection') && strcmpi(opts.detection, 'cascade')
    % 18 nov. (demandé par l'utilisateur) : événements de la déconvolution
    % CASCADE (modèle GCaMP8), sur DFF0 non lissé
    [Raster, Acttmp2, th_detection, PromALL, WidthALL, info_bouffees] = ...
        transient_detection_cascade(DFF0, WinRest, path, sampling_rate, opts);
elseif isfield(opts, 'detection') && strcmpi(opts.detection, 'cascade_discrets')
    % 19 nov. (demandé par l'utilisateur) : potentiels d'action discrets de
    % CASCADE (infer_discrete_spikes) ; 1 par image avec au moins 1 PA, au repos
    [Raster, Acttmp2, th_detection, PromALL, WidthALL, info_bouffees] = ...
        transient_detection_cascade_discrets(DFF0, WinRest, path, sampling_rate, opts);
elseif isfield(opts, 'detection') && ~strcmpi(opts.detection, 'standard')
    % 18 nov. (test demandé par l'utilisateur) : exclusion des bouffées
    % d'après le code d'Arnaud ; sans le champ opts.detection, rien ne change
    [Raster, Acttmp2, th_detection, PromALL, WidthALL, info_bouffees] = ...
        transient_detection_bouffees(Tr1b, WinRest, DFF0, sampling_rate, opts);
else
    [Raster, Acttmp2, th_detection, PromALL, WidthALL] = ...
        transient_detection_shot_theta2(Tr1b, MinPeakDistance, threshold_peak, WinRest, DFF0, ancien);
end

%% CELLULES HYPERACTIVES
transients_percell = cellfun(@length, Acttmp2);
seuil_frequence = 0.2;    % 1 pic toutes les 5 s
min_hyper = 200;          % jamais moins de 200 transients
% 21 nov. (demandé par l'utilisateur) : seuil réglable, pour CASCADE (où
% 0,2 Hz retire jusqu'à 2/3 des cellules) ; défauts = comme avant
if isfield(opts, 'hyperactives_seuil_hz') && ~isempty(opts.hyperactives_seuil_hz)
    seuil_frequence = opts.hyperactives_seuil_hz;
end
if isfield(opts, 'hyperactives_min') && ~isempty(opts.hyperactives_min)
    min_hyper = opts.hyperactives_min;
end
max_peak = max(length(WinRest) / sampling_rate * seuil_frequence, min_hyper);
cellules_hyperactives_idx = find(transients_percell > max_peak);
if isfield(opts, 'retirer_hyperactives') && opts.retirer_hyperactives && ~isempty(cellules_hyperactives_idx)
    Raster(cellules_hyperactives_idx, :) = 0;
end

%% ROI TROP PETITES OU TROP GRANDES (tailles_roi), puis DOUBLONS (même
% cellule comptée deux fois, detecter_doublons)
% Une ROI retirée garde sa trace mais n'a plus de transient : elle n'entre
% plus dans les SCE ni dans les assemblées, et les étapes suivantes la
% traitent comme non colorée (couleurs_sans_doublons).
tailles = struct('npix', [], 'mediane', NaN, 'petites', zeros(1, 0), 'grandes', zeros(1, 0), ...
    'retirees', zeros(1, 0), 'min_px', NaN, 'max_facteur', NaN);
doublons = struct('retirees', zeros(1, 0), 'paires', table(), 'r_min', NaN);
faire_tailles = isfield(opts, 'retirer_tailles') && opts.retirer_tailles;
faire_doublons = isfield(opts, 'retirer_doublons') && opts.retirer_doublons;
if faire_tailles || faire_doublons
    [~, ~, ~, ~, ~, Coomasqx, Coomasqy] = charger_geometrie(path);
    if numel(Coomasqx) ~= NCell
        warning('Fall.mat : %d cellules, activité : %d ; tailles et doublons non vérifiés.', ...
            numel(Coomasqx), NCell);
    else
        if faire_tailles
            min_px = []; max_facteur = [];
            if isfield(opts, 'roi_min_px'), min_px = opts.roi_min_px; end
            if isfield(opts, 'roi_max_facteur'), max_facteur = opts.roi_max_facteur; end
            tailles = tailles_roi(Coomasqx, min_px, max_facteur);
            Raster(tailles.retirees, :) = 0;
        end
        if faire_doublons
            r_min = 0.7;
            if isfield(opts, 'doublons_r_min') && ~isempty(opts.doublons_r_min), r_min = opts.doublons_r_min; end
            Tr_d = Tr1b;
            Tr_d(tailles.retirees, :) = 0;          % une ROI déjà retirée n'est le doublon de personne
            doublons = detecter_doublons(Tr_d, WinRest, Coomasqx, Coomasqy, r_min);
            Raster(doublons.retirees, :) = 0;
        end
    end
end

%% POPULATION (15 nov., facultatif) : SCE et assemblées cherchés dans une
% partie des cellules. 'sequence' : cellules activées en séquence pendant la
% course (même calcul que l'étape 5) ; 'hors_sequence' : toutes les autres.
% Le raster des cellules écartées est mis à zéro : elles n'entrent plus dans
% les SCE, ni dans le calcul du seuil, ni dans les assemblées, et les étapes
% suivantes les traitent comme non colorées (couleurs_sans_doublons) : les
% couleurs sont donc mélangées parmi les cellules gardées seulement.
population = struct('regle', 'toutes', 'retirees', zeros(1, 0), 'dans_sequence', [], 'n_courses', NaN, ...
    'seuil', NaN);
if isfield(opts, 'population') && ~isempty(opts.population) && ~strcmp(opts.population, 'toutes')
    opts_seq = options_sequences(struct());
    if isfield(opts, 'population_seq') && isstruct(opts.population_seq)
        opts_seq = options_sequences(struct('seq', opts.population_seq));
    end
    [dans_sequence, info_seq] = cellules_sequence(F, speed, opts_seq);
    if ~isempty(info_seq.statut)
        error('preprocessing_7:population', ['config.pre.population = ''%s'' : cellules de la séquence ' ...
            'impossibles à trouver (%s). Retirer cette session (inclure = 0) pour cette analyse.'], ...
            opts.population, info_seq.statut);
    end
    switch opts.population
        case 'sequence'
            population.retirees = find(~dans_sequence);
        case 'hors_sequence'
            population.retirees = find(dans_sequence);
        otherwise
            error('preprocessing_7:population', ['config.pre.population = ''%s'' inconnu : ''sequence'' ou ' ...
                '''hors_sequence'' (ou ne pas mettre le champ : toutes les cellules).'], opts.population);
    end
    population.regle = opts.population;
    population.dans_sequence = dans_sequence;
    population.n_courses = info_seq.n_courses;
    population.seuil = info_seq.seuil;
    Raster(population.retirees, :) = 0;
    fprintf('Population ''%s'' : %d cellules gardées sur %d (%d dans la séquence, %d courses)\n', ...
        opts.population, NCell - numel(population.retirees), NCell, sum(dans_sequence), info_seq.n_courses);
end

%% DÉTECTION DES SCE (Synchronous Events)
% Les modèles nuls des SCE décalent l'activité à l'intérieur du repos seulement
opts.sampling_rate = sampling_rate;
[Race, SumAct, MAct, NRace, RasterRace, TRace, Raster, max_cells_allowed, opts, diag_SCE] = ...
    SCE_detection5(Raster, opts, bad_frames_etendues, graine, WinRest_bin);

%% Acttmp2, PromALL, WidthALL : on retire les transients mis à zéro
for i = 1:NCell
    garde = Raster(i, Acttmp2{i}) > 0;
    Acttmp2{i}  = Acttmp2{i}(garde);
    PromALL{i}  = PromALL{i}(garde);
    WidthALL{i} = WidthALL{i}(garde);
end
n_transients_total = sum(cellfun(@length, Acttmp2));

%% Diagnostic (à regarder, jamais avec les couleurs)
diagnostic = diag_SCE;
diagnostic.NRace                     = NRace;
diagnostic.transients_percell        = cellfun(@length, Acttmp2);
diagnostic.cellules_hyperactives_idx = cellules_hyperactives_idx;
diagnostic.doublons                  = doublons;
diagnostic.tailles_roi               = tailles;
diagnostic.population                = population;
diagnostic.n_bad_frames_etendues     = sum(bad_frames_etendues);
diagnostic.n_images_rejetees_manuelles = n_rejet_manuel;   % 23 nov.
diagnostic.duree_repos_s             = numel(WinRest) / sampling_rate;
if ~isempty(info_bouffees), diagnostic.bouffees = info_bouffees; end   % 18 nov.
end
