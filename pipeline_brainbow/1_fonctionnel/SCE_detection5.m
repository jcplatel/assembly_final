function [Race,SumAct,MAct,NRace,RasterRace,TRace,Raster,max_cells_allowed,opts,diag_SCE] = ...
    SCE_detection5(Raster,opts,bad_frames,graine,repos)
% Version 5 de SCE_detection4. Changements :
%  - opts.sce_n_cells_threshold n'est plus écrasé ; il ne sert qu'à la règle
%    'fixe' (8 nov., demandé par l'utilisateur : il n'a jamais été censé être
%    un minimum ; avant, il servait aussi de minimum aux règles 'percentile'
%    et 'fdr', et 14 bloquait le seuil de 15 sessions sur 26 au 95e percentile).
%    Minimum des règles 'percentile' et 'fdr' : opts.sce_minimum (5 par défaut).
%    Le seuil réellement utilisé est rangé dans opts.sce_n_cells_threshold_calcule
%  - deux règles de seuil (opts.sce_regle) :
%      'percentile' : percentile opts.sce_percentile de la co-activité dans des
%                     shuffles locaux (+/- 10 trames), comme avant mais sans
%                     plafond à 15. ATTENTION : ces shuffles gardent en partie
%                     les vrais SCE (étalés sur 21 trames) ; plus une session
%                     a de grands SCE, plus le seuil monte.
%      'fixe'       : seuil = opts.sce_n_cells_threshold, sans calcul
%      'fdr'        : plus petit seuil pour lequel les SCE attendus par hasard
%                     ne dépassent pas opts.sce_fdr (ex. 5 %) des SCE détectés.
%                     Hasard : chaque cellule décalée d'un bloc d'au plus
%                     opts.sce_decalage_fdr_s secondes (Inf = circulaire).
%  - diag_SCE.n_SCE_hasard = SCE attendus par hasard au seuil retenu
%    (décalage circulaire de chaque cellule, à l'intérieur du repos)
%  - MAct calculé avec movmax (même résultat, beaucoup plus rapide)
%  - shuffles reproductibles (graine)
%
%  repos (facultatif) : masque logique des trames où des transients peuvent
%  être détectés (repos, hors bad frames). Les modèles nuls décalent
%  l'activité à l'intérieur de ces trames seulement.

if nargin < 4, graine = 1; end
sce_n_cells_threshold = opts.sce_n_cells_threshold;
MinPeakDistancesce    = opts.MinPeakDistancesce;
synchronous_frames    = opts.synchronous_frames;
if isfield(opts, 'sce_percentile'), sce_percentile = opts.sce_percentile; else, sce_percentile = 99; end
if isfield(opts, 'sce_regle'), sce_regle = opts.sce_regle; else, sce_regle = 'percentile'; end
if isfield(opts, 'sce_minimum') && ~isempty(opts.sce_minimum), sce_minimum = opts.sce_minimum; else, sce_minimum = 5; end

max_jitter = 10;            % Décalage max (+/- 10 frames)
NShfl      = 100;           % Nombre de shuffles pour le seuil
NNull      = 20;            % surrogates pour les courbes « SCE attendus par hasard »

[NCell, Nz] = size(Raster);

% Plafond maximal absolu de cellules pour un SCE
max_cells_SCE     = 150;
max_cells_allowed = max_cells_SCE;
opts.max_cells_SCE = max_cells_SCE;

% Trames utilisables : repos (si fourni), hors bad frames
if nargin < 5 || isempty(repos)
    repos = true(1, Nz);
end
repos = logical(repos(:)');
masque_bad = false(1, Nz);
masque_bad(bad_frames) = true;
repos = repos & ~masque_bad;

% 0. Suppression des bad frames connues (déjà élargies par preprocessing_7)
if ~isempty(bad_frames)
    Raster(:, bad_frames) = 0;
end

% =====================================================================
% 1. NETTOYAGE PRÉVENTIF DES ARTEFACTS INSTANTANÉS
% =====================================================================
SumAct = sum(Raster, 1);
toxic_frames = SumAct > max_cells_allowed;
n_trames_artefact = sum(toxic_frames);   % gardé dans le diagnostic (5 nov.)
if any(toxic_frames)
    fprintf('Nettoyage Raster : %d frames instantanées (> %d cellules) mises à zéro.\n', ...
        sum(toxic_frames), max_cells_allowed);
    toxic_frames_padded = conv(double(toxic_frames), [1 1 1], 'same') > 0;
    Raster(:, toxic_frames_padded) = 0;
    SumAct = sum(Raster, 1);
end

% =====================================================================
% 2. VECTEUR D'ACTIVITÉ CO-SYNCHRONE (MAct)
% =====================================================================
MAct = activite_cosynchrone(Raster, synchronous_frames);

% =====================================================================
% 3. SEUIL STATISTIQUE
% =====================================================================
seuil_shuffles = NaN;
courbe_fdr = [];
switch sce_regle
    case 'percentile'
        % Shuffling local (Local Jittering), comme avant
        Sumactsh = zeros(Nz - synchronous_frames, NShfl);
        parfor n = 1:NShfl
            fixer_flux_aleatoire(graine, n, NShfl);
            Rastersh = zeros(NCell, Nz);
            for c = 1:NCell
                k = randi([-max_jitter, max_jitter]);
                if k > 0
                    Rastersh(c, k+1:end) = Raster(c, 1:end-k);
                elseif k < 0
                    Rastersh(c, 1:end+k) = Raster(c, 1-k:end);
                else
                    Rastersh(c, :) = Raster(c, :);
                end
            end
            Sumactsh(:, n) = activite_cosynchrone(Rastersh, synchronous_frames)';
        end
        % percentile des valeurs non nulles des shuffles, jamais sous opts.sce_minimum
        Toutes_valeurs_shuffles  = Sumactsh(:);
        Valeurs_actives_shuffles = Toutes_valeurs_shuffles(Toutes_valeurs_shuffles > 0);
        seuil_shuffles = ceil(prctile(Valeurs_actives_shuffles, sce_percentile));
        sce_n_cells_threshold = max(seuil_shuffles, sce_minimum);

    case 'fixe'
        % seuil imposé : opts.sce_n_cells_threshold (ancien code : plafond à 15)

    case 'fdr'
        if isfield(opts, 'sampling_rate'), fs = opts.sampling_rate; else, fs = 30; end
        decalage = round(opts.sce_decalage_fdr_s * fs);       % Inf reste Inf
        opts_courbe = opts;
        opts_courbe.sce_n_cells_threshold = sce_minimum;      % courbe à partir du minimum
        [seuils, N_obs, N_hasard] = courbe_SCE_hasard(Raster, repos, opts_courbe, decalage, NNull, graine);
        trop = N_hasard > opts.sce_fdr * N_obs;               % trop de hasard à ce seuil
        dernier = find(trop, 1, 'last');
        if isempty(dernier)
            sce_n_cells_threshold = seuils(1);
        else
            sce_n_cells_threshold = seuils(dernier) + 1;      % au-dessus : toujours <= sce_fdr
        end
        courbe_fdr = struct('seuils', seuils, 'N_obs', N_obs, 'N_hasard', N_hasard, ...
            'decalage_trames', decalage, 'fdr', opts.sce_fdr);

    otherwise
        error('opts.sce_regle doit valoir ''percentile'', ''fixe'' ou ''fdr''.');
end
opts.sce_n_cells_threshold_calcule = sce_n_cells_threshold;

% =====================================================================
% 4. DÉTECTION DES PICS AVEC DOUBLE BORNE (MIN ET MAX)
% =====================================================================
if numel(MAct) > MinPeakDistancesce + 1 && max(MAct) > sce_n_cells_threshold
    [pks, locs] = findpeaks(MAct, ...
        'MinPeakHeight', sce_n_cells_threshold, ...
        'MinPeakDistance', MinPeakDistancesce);
else
    pks = []; locs = [];
end
valid_peaks = (pks <= max_cells_allowed);
TRace = locs(valid_peaks);
NRace = length(TRace);

% Faux SCE attendus au seuil retenu : décalage circulaire COMPLET de chaque
% cellule à l'intérieur du repos (détruit toute co-activation, garde
% l'activité de chaque cellule). Les shuffles locaux ne servent pas à ça :
% ils gardent en partie les vraies co-activations.
opts_c = opts;
opts_c.sce_n_cells_threshold = min(sce_minimum, sce_n_cells_threshold);   % la courbe doit contenir le seuil retenu
[seuils_c, ~, N_circ] = courbe_SCE_hasard(Raster, repos, opts_c, Inf, NNull, graine + 1);
i_seuil = find(seuils_c == sce_n_cells_threshold, 1);
if isempty(i_seuil), n_hasard = 0; else, n_hasard = N_circ(i_seuil); end
diag_SCE.regle_seuil           = sce_regle;
diag_SCE.seuil_SCE             = sce_n_cells_threshold;
diag_SCE.seuil_shuffles        = seuil_shuffles;
diag_SCE.n_SCE_hasard          = n_hasard;
diag_SCE.n_trames_artefact     = n_trames_artefact;   % trames à plus de max_cells_SCE cellules, mises à zéro
diag_SCE.proportion_SCE_hasard = n_hasard / max(NRace, 1);
diag_SCE.courbe_fdr            = courbe_fdr;

% =====================================================================
% 5. EXTRACTION DE RACE (cellules x SCE) ET RASTERRACE (cellules x trames)
% =====================================================================
Race = zeros(NCell, NRace);
RasterRace = zeros(NCell, Nz);
for i = 1:NRace
    idx_start = TRace(i);
    idx_end   = min(Nz, TRace(i) + synchronous_frames);
    Race(:, i) = max(Raster(:, idx_start:idx_end), [], 2);
    for c = 1:NCell
        if Race(c, i) == 1
            local_frames = idx_start:idx_end;
            active_frame = local_frames(Raster(c, local_frames) == 1);
            if ~isempty(active_frame)
                RasterRace(c, active_frame(1)) = 1;
            end
        end
    end
end

fprintf('Ncells = %d | transients : %d | SCE détectés : %d (%.1f attendus par hasard, seuil %d cellules, règle %s)\n', ...
    NCell, sum(SumAct), NRace, diag_SCE.n_SCE_hasard, sce_n_cells_threshold, sce_regle);
end
