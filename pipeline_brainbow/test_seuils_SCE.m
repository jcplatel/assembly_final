function [T_session, T_detail] = test_seuils_SCE(config, percentiles, seuils_fixes, minimum_absolu)
%TEST_SEUILS_SCE  Seuil des SCE (nombre de cellules) et nombre de SCE obtenus
%   avec les 3 règles, pour toutes les sessions déjà prétraitées (étape 1).
%
%   test_seuils_SCE                                  % config_brainbow, 90 95 99
%   test_seuils_SCE(config_brainbow('ancien_rang'))  % autre dossier de l'étape 1
%   test_seuils_SCE([], [90 95 99], [10 14 20], 5)
%   config.seulement = {'444119'}; test_seuils_SCE(config)   % quelques sessions
%
%   Ne refait PAS le prétraitement : relit pretraitement.mat (Raster, repos,
%   bad frames) et refait seulement le seuil et la détection des SCE, avec
%   le même calcul que SCE_detection5 (mêmes shuffles, même graine). Le
%   réglage utilisé à l'étape 1 retombe donc sur le même seuil (colonnes
%   etape1_seuil et etape1_nSCE pour vérifier).
%   Aucune couleur n'est lue. Rien n'est écrit dans les dossiers des sessions.
%
%   Règles testées :
%     percentile p : p-ième percentile de la co-activité dans 100 shuffles
%                    locaux (+/- 10 trames)
%     fdr q        : plus petit seuil où les SCE attendus par hasard
%                    (décalage d'au plus config.pre.sce_decalage_fdr_s) restent
%                    sous q des SCE détectés ; q = 1 - p/100 :
%                    90 -> 10 %, 95 -> 5 %, 99 -> 1 %
%     fixe s       : seuil imposé, pour chaque valeur de seuils_fixes
%                    (défaut : config.pre.sce_n_cells_threshold)
%   percentile et fdr : jamais sous minimum_absolu cellules (défaut 5, comme
%   l'interface).
%
%   Pour chaque session et chaque réglage :
%     seuil         nombre de cellules co-actives au-dessus duquel un pic est un SCE
%     n_SCE         SCE détectés (comme à l'étape 1)
%     SCE_par_min   par minute de repos
%     n_SCE_hasard  SCE attendus par hasard au même seuil (chaque cellule
%                   décalée au hasard dans le repos, comme diagnostic.n_SCE_hasard)
%     pct_hasard    n_SCE_hasard / n_SCE
%
%   Sorties : T_session (une ligne par session) et T_detail (une ligne par
%   session et par réglage), aussi écrites dans
%   <dossier de l'étape 1>\test_seuils_SCE.xlsx (feuilles 'par_session' et 'detail').
%   Durée : environ une minute par session (100 shuffles en parallèle).

if nargin < 1 || isempty(config), config = config_brainbow(); end
if nargin < 2 || isempty(percentiles), percentiles = [90 95 99]; end
if nargin < 3 || isempty(seuils_fixes), seuils_fixes = config.pre.sce_n_cells_threshold; end
if nargin < 4 || isempty(minimum_absolu), minimum_absolu = 5; end
dossier_pipeline = fileparts(mfilename('fullpath'));
addpath(genpath(dossier_pipeline));
etat_rng = rng;
remettre_rng = onCleanup(@() rng(etat_rng)); %#ok<NASGU>

% Mêmes constantes que SCE_detection5
NShfl      = 100;      % shuffles de la règle percentile
max_jitter = 10;       % décalage max des shuffles locaux (trames)
NNull      = 20;       % surrogates des courbes « SCE attendus par hasard »

% Réglages testés : {règle, paramètre, étiquette}
reglages = cell(0, 3);
for p = percentiles(:)'
    reglages(end+1, :) = {'percentile', p, sprintf('perc%g', p)}; %#ok<AGROW>
end
for p = percentiles(:)'
    q = (100 - p) / 100;
    reglages(end+1, :) = {'fdr', q, sprintf('fdr%g', 100 - p)}; %#ok<AGROW>
end
for sf = seuils_fixes(:)'
    reglages(end+1, :) = {'fixe', sf, sprintf('fixe%g', sf)}; %#ok<AGROW>
end
reglages(:, 3) = matlab.lang.makeValidName(reglages(:, 3));
nR = size(reglages, 1);

% Sessions (étape 1 lue dans le dossier de l'étape 1)
config = dossiers_analyse(config);
sessions = lire_table_sessions(config);
if isfield(config, 'seulement') && ~isempty(config.seulement)
    seulement = cellstr(config.seulement);
    sessions = sessions(ismember(sessions.nom, seulement) | ismember(sessions.id_animal, seulement), :);
end
nS = height(sessions);
if isempty(gcp('nocreate')), parpool; end
fprintf('=== Seuil des SCE : %d sessions, étape 1 lue dans %s ===\n', nS, config.PathEtape1);
fprintf('percentile %s | fdr %s %% | fixe %s | minimum %d cellules\n\n', mat2str(percentiles), ...
    mat2str(100 - percentiles), mat2str(seuils_fixes), minimum_absolu);

detail = cell(0, 12);
par_session = cell(nS, 1);
sautees = {};
for s = 1:nS
    nom = sessions.nom{s};
    f = fullfile(sessions.dossier_source{s}, 'pretraitement.mat');
    if ~isfile(f)
        sautees{end+1} = nom; %#ok<AGROW>
        fprintf('%d/%d %s : pas de prétraitement, sautée\n', s, nS, nom);
        continue
    end
    t0 = tic;
    P = load(f, 'Raster', 'opts', 'WinRest', 'bad_frames', 'sampling_rate', 'TRace');
    Raster = double(P.Raster);
    [NCell, Nz] = size(Raster);
    opts = P.opts;
    fs = P.sampling_rate;
    sync = opts.synchronous_frames;
    mpd = opts.MinPeakDistancesce;
    max_cells = 150;
    if isfield(opts, 'max_cells_SCE'), max_cells = opts.max_cells_SCE; end
    graine = sessions.graine(s);

    % repos (comme dans SCE_detection5) : épisodes de repos, hors bad frames
    repos = false(1, Nz);
    repos(P.WinRest) = true;
    if isfield(P, 'bad_frames') && ~isempty(P.bad_frames)
        if islogical(P.bad_frames) && numel(P.bad_frames) == Nz
            repos = repos & ~P.bad_frames(:)';
        else
            repos(P.bad_frames) = false;
        end
    end
    repos_min = sum(repos) / fs / 60;

    % Co-activité des données (MAct), comme à l'étape 1
    MAct = activite_cosynchrone(Raster, sync);

    % Règle percentile : les MÊMES 100 shuffles locaux pour tous les percentiles
    Sumactsh = zeros(Nz - sync, NShfl);
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
        Sumactsh(:, n) = activite_cosynchrone(Rastersh, sync)';
    end
    valeurs_shuffles = Sumactsh(:);
    valeurs_shuffles = valeurs_shuffles(valeurs_shuffles > 0);

    % Règle fdr : une seule courbe (données et hasard) pour tous les q
    opts_c = opts;
    opts_c.max_cells_SCE = max_cells;
    opts_c.sce_n_cells_threshold = minimum_absolu;
    decalage = round(config.pre.sce_decalage_fdr_s * fs);
    [seuils_f, N_obs_f, N_hasard_f] = courbe_SCE_hasard(Raster, repos, opts_c, decalage, NNull, graine);

    % SCE attendus par hasard (décalage circulaire), même graine que l'étape 1
    opts_c.sce_n_cells_threshold = min([minimum_absolu, seuils_fixes(:)']);
    [seuils_c, ~, N_circ] = courbe_SCE_hasard(Raster, repos, opts_c, Inf, NNull, graine + 1);

    ligne = cell(1, 2 * nR);
    resume = {};
    for r = 1:nR
        regle = reglages{r, 1};
        param = reglages{r, 2};
        switch regle
            case 'percentile'
                seuil = max(ceil(prctile(valeurs_shuffles, param)), minimum_absolu);
            case 'fdr'
                trop = N_hasard_f > param * N_obs_f;
                dernier = find(trop, 1, 'last');
                if isempty(dernier), seuil = seuils_f(1); else, seuil = seuils_f(dernier) + 1; end
            case 'fixe'
                seuil = param;
        end
        n_SCE = nombre_SCE(MAct, seuil, mpd, max_cells);
        i_c = find(seuils_c == seuil, 1);
        if isempty(i_c), n_hasard = 0; else, n_hasard = N_circ(i_c); end
        if strcmp(regle, 'fdr'), param_affiche = 100 * param; else, param_affiche = param; end
        detail(end+1, :) = {nom, sessions.id_animal{s}, NCell, round(repos_min, 1), regle, ...
            param_affiche, seuil, n_SCE, round(n_SCE / max(repos_min, eps), 2), ...
            round(n_hasard, 1), round(100 * n_hasard / max(n_SCE, 1), 1), reglages{r, 3}}; %#ok<AGROW>
        ligne(2 * r - 1 : 2 * r) = {seuil, n_SCE};
        resume{end+1} = sprintf('%s %d (%d)', reglages{r, 3}, seuil, n_SCE); %#ok<AGROW>
    end

    % Ce que l'étape 1 a réellement utilisé (pour vérifier)
    seuil_e1 = NaN;
    if isfield(opts, 'sce_n_cells_threshold_calcule'), seuil_e1 = opts.sce_n_cells_threshold_calcule; end
    regle_e1 = 'percentile';
    if isfield(opts, 'sce_regle'), regle_e1 = opts.sce_regle; end
    nSCE_e1 = NaN;
    if isfield(P, 'TRace'), nSCE_e1 = numel(P.TRace); end
    par_session{s} = [{nom, sessions.id_animal{s}, NCell, round(repos_min, 1)}, ligne, ...
        {regle_e1, seuil_e1, nSCE_e1}];

    fprintf('%d/%d %s : %d cellules, %.1f min de repos (%.0f s)\n   seuil (nb SCE) : %s\n', ...
        s, nS, nom, NCell, repos_min, toc(t0), strjoin(resume, ' | '));
end

% Tableaux
noms_detail = {'session', 'animal', 'n_cellules', 'repos_min', 'regle', 'parametre', ...
    'seuil', 'n_SCE', 'SCE_par_min', 'n_SCE_hasard', 'pct_hasard', 'reglage'};
T_detail = cell2table(detail, 'VariableNames', noms_detail);
colonnes = cell(1, 2 * nR);
for r = 1:nR
    colonnes(2 * r - 1 : 2 * r) = {[reglages{r, 3} '_seuil'], [reglages{r, 3} '_nSCE']};
end
noms_session = [{'session', 'animal', 'n_cellules', 'repos_min'}, colonnes, ...
    {'etape1_regle', 'etape1_seuil', 'etape1_nSCE'}];
par_session = par_session(~cellfun(@isempty, par_session));
T_session = cell2table(vertcat(par_session{:}), 'VariableNames', noms_session);
if isempty(T_session)
    fprintf('\nAucune session prétraitée.\n');
    return
end

% Résumé sur toutes les sessions : médiane [min - max]
fprintf('\n=== Toutes les sessions (%d) : médiane [min - max] ===\n', height(T_session));
fprintf('%-12s %-22s %-22s %-20s\n', 'réglage', 'seuil (cellules)', 'nombre de SCE', 'part du hasard (%)');
for r = 1:nR
    m = strcmp(T_detail.reglage, reglages{r, 3});
    fprintf('%-12s %-22s %-22s %-20s\n', reglages{r, 3}, med_min_max(T_detail.seuil(m)), ...
        med_min_max(T_detail.n_SCE(m)), med_min_max(T_detail.pct_hasard(m)));
end
if ~isempty(sautees)
    fprintf('\nSessions sans prétraitement (sautées) : %s\n', strjoin(sautees, ', '));
end

% Excel dans le dossier de l'étape 1
fichier = fullfile(config.PathEtape1, 'test_seuils_SCE.xlsx');
try
    if isfile(fichier), delete(fichier); end
    writetable(T_session, fichier, 'Sheet', 'par_session');
    writetable(T_detail, fichier, 'Sheet', 'detail');
    fprintf('\nTableaux : %s\n', fichier);
catch err
    warning('test_seuils_SCE.xlsx non écrit (ouvert dans Excel ?) : %s', err.message);
end
end


function n = nombre_SCE(MAct, seuil, mpd, max_cells)
% Nombre de SCE au seuil donné : même détection que SCE_detection5 (étape 4)
n = 0;
if numel(MAct) > mpd + 1 && max(MAct) > seuil
    pks = findpeaks(MAct, 'MinPeakHeight', seuil, 'MinPeakDistance', mpd);
    n = sum(pks <= max_cells);
end
end


function t = med_min_max(x)
x = x(~isnan(x));
if isempty(x), t = '-'; return; end
t = sprintf('%g [%g - %g]', round(median(x), 1), min(x), max(x));
end
