function [Raster, Acttmp2, th_detection, PromALL, WidthALL, info] = ...
    transient_detection_cascade_discrets(DFF0, WinRest, path, sampling_rate, opts)
%% =========================================================================
%  TRANSIENT_DETECTION_CASCADE_DISCRETS (19 nov. 2026, demandée par l'utilisateur)
%  Raster = potentiels d'action (PA) discrets de CASCADE :
%    1. spike_prob = CASCADE sur DFF0 NON lissé (cascade_spikes, gardé dans
%       <path>\cascade_<modele>.mat)
%    2. PA entiers placés par infer_discrete_spikes de CASCADE
%       (cascade_discrets.py ; aucun réglage : noyau du modèle ; graine
%       fixée), gardés dans <path>\cascade_discrets_<modele>.mat
%    3. Raster(i, t) = 1 si au moins 1 PA à l'image t, au repos seulement
%  Choix de l'utilisateur (19 nov.) : 1 par image avec >= 1 PA ; pas de
%  regroupement en événements ; repos seulement.
%  Sorties comme transient_detection_shot_theta2 : Acttmp2{i} = images avec
%  PA ; PromALL{i} = nombre de PA à ces images ; WidthALL{i} = 1 ;
%  th_detection = NaN (pas de seuil).
%  opts : cascade_modele, cascade_processus (défaut 24), cascade_graine (1),
%         cascade_reechantillonnage_Hz (défaut [] : aucun),
%         cascade_min_pa (défaut 1 : image comptée avec >= 1 PA ; 2 : >= 2 PA).
%         cascade_max_pa_evenement, cascade_max_duree_evenement_s (22 nov. ;
%         défaut Inf) : gros événements retirés (voir plus bas).
%  20 nov. (demandé par l'utilisateur) : cascade_reechantillonnage_Hz = 15
%  pour un modèle entraîné à 15 Hz (ex. GC8m_EXC_15Hz_smoothing50ms_high_noise) :
%  DFF0 interpolé (linéaire) sur une grille à 15 Hz, CASCADE et PA discrets à
%  15 Hz, puis chaque PA ramené à l'image d'origine la plus proche dans le
%  temps (somme des PA par image). Caches suffixés _<Hz>Hz (PA à 15 Hz).
%  L'interpolation n'ajoute pas d'information temporelle et lisse le bruit.
%  =========================================================================
if isfield(opts, 'cascade_modele'), modele = opts.cascade_modele;
else, modele = 'GC8_EXC_10Hz_smoothing150ms_high_noise'; end
if isfield(opts, 'cascade_processus'), n_proc = opts.cascade_processus; else, n_proc = 24; end
if isfield(opts, 'cascade_graine'), graine = opts.cascade_graine; else, graine = 1; end
hz = [];
if isfield(opts, 'cascade_reechantillonnage_Hz'), hz = opts.cascade_reechantillonnage_Hz; end
[NCell0, Nz0] = size(DFF0);
suffixe = '';
if ~isempty(hz)
    suffixe = sprintf('_%gHz', hz);
    t0 = (0:Nz0 - 1) / sampling_rate;
    t_r = 0:1 / hz:t0(end);
    DFF0 = interp1(t0, double(DFF0)', t_r, 'linear')';
    if NCell0 == 1, DFF0 = DFF0(:)'; end
    image_origine = min(Nz0, round(t_r * sampling_rate) + 1);   % image d'origine la plus proche
end
o = opts; o.cascade_modele = modele;
sp = cascade_spikes(DFF0, o, fullfile(path, ['cascade_' modele suffixe '.mat']));

% PA discrets (gardés : relus si même spike_prob, même modèle, même graine)
fc = fullfile(path, ['cascade_discrets_' modele suffixe '.mat']);
controle = [sum(sp(:), 'omitnan'), sum(abs(sp(:)), 'omitnan')];
n_pa = [];
if isfile(fc)
    C = load(fc);
    if strcmp(C.modele, modele) && isequal(size(C.n_pa), size(sp)) && C.graine == graine && ...
            all(abs(C.controle - controle) <= 1e-6 * max(1, abs(controle)))
        n_pa = double(C.n_pa);
    end
end
if isempty(n_pa)
    script = fullfile(fileparts(mfilename('fullpath')), 'cascade_discrets.py');
    python = 'C:\Users\jcplatel\Desktop\matlab\Cascade\venv_cascade\Scripts\python.exe';
    dossier = 'C:\Users\jcplatel\Desktop\matlab\Cascade';
    if isfield(opts, 'cascade_python'), python = opts.cascade_python; end
    if isfield(opts, 'cascade_dossier'), dossier = opts.cascade_dossier; end
    base = tempname; entree = [base '_e.mat']; sortie = [base '_s.mat'];
    spike_prob = sp; %#ok<NASGU>
    save(entree, 'spike_prob', 'modele', 'graine', '-v7');
    env_conda = 'C:\Users\jcplatel\.conda\envs\suite2p';
    ancien = getenv('PATH');
    setenv('PATH', [env_conda '\Library\bin;' env_conda '\DLLs;' ancien]);
    nettoyer = onCleanup(@() setenv('PATH', ancien));
    [st, txt] = system(sprintf('"%s" "%s" "%s" "%s" "%s" %d', python, script, entree, sortie, dossier, n_proc));
    delete(entree);
    if st ~= 0 || ~isfile(sortie), error('PA discrets CASCADE : échec :\n%s', txt); end
    S = load(sortie, 'n_pa'); delete(sortie);
    n_pa = double(S.n_pa);
    n_pa_s = single(n_pa); %#ok<NASGU>
    sauver(fc, n_pa_s, modele, graine, controle);
end
if ~isempty(hz)
    % PA à 15 Hz -> images d'origine (somme des PA tombés dans chaque image)
    M = sparse(1:numel(image_origine), image_origine, 1, numel(image_origine), Nz0);
    n_pa = full(double(n_pa) * M);
end

[NCell, Nz] = size(n_pa);
repos = false(1, Nz); repos(WinRest) = true;
% 21 nov. (demandé par l'utilisateur) : opts.cascade_min_pa = nombre minimal
% de PA dans l'image pour compter un événement (défaut 1 = comme avant)
min_pa = 1;
if isfield(opts, 'cascade_min_pa') && ~isempty(opts.cascade_min_pa), min_pa = opts.cascade_min_pa; end
% 22 nov. (demandé par l'utilisateur) : retirer les GROS événements (une
% cellule tire 0 à 2 PA par ripple ; une bouffée longue ou riche en PA n'est
% probablement pas liée à un ripple). Événement = suite d'images consécutives
% avec >= 1 PA (session entière). Il est retiré en entier si :
%   opts.cascade_max_pa_evenement      : nombre total de PA > cette valeur
%   opts.cascade_max_duree_evenement_s : durée > cette valeur (en images :
%                                        round(durée x fréquence), ex. 0,3 s -> 3 images)
% Absents ou Inf = rien n'est retiré (comme avant). Champs facultatifs, pas
% dans config_brainbow.m : les étapes 1 déjà faites restent valables.
max_pa = Inf; max_images = Inf;
if isfield(opts, 'cascade_max_pa_evenement') && ~isempty(opts.cascade_max_pa_evenement)
    max_pa = opts.cascade_max_pa_evenement;
end
if isfield(opts, 'cascade_max_duree_evenement_s') && ~isempty(opts.cascade_max_duree_evenement_s)
    max_images = round(opts.cascade_max_duree_evenement_s * sampling_rate);
end
info_gros = struct('n_evenements', 0, 'n_retires', 0, 'PA_retires', 0, 'max_pa', max_pa, 'max_images', max_images);
if isfinite(max_pa) || isfinite(max_images)
    d = diff([false(NCell, 1), n_pa > 0, false(NCell, 1)], 1, 2);
    [debut, ci] = find(d' == 1);          % ordre : cellule puis temps (débuts et fins)
    [fin, ~] = find(d' == -1);
    fin = fin - 1;
    cs = cumsum([zeros(NCell, 1), n_pa], 2);
    taille = cs(sub2ind(size(cs), ci, fin + 1)) - cs(sub2ind(size(cs), ci, debut));
    duree = fin - debut + 1;
    gros = find(taille > max_pa | duree > max_images);
    info_gros.n_evenements = numel(debut);
    info_gros.n_retires = numel(gros);
    info_gros.PA_retires = sum(taille(gros));
    for e = gros(:)'
        n_pa(ci(e), debut(e):fin(e)) = 0;
    end
    fprintf('Gros événements retirés : %d / %d (%.1f %%), %.1f %% des PA (max %g PA, %g images)\n', ...
        numel(gros), numel(debut), 100 * numel(gros) / max(1, numel(debut)), ...
        100 * info_gros.PA_retires / max(1, sum(taille)), max_pa, max_images);
end
Raster = double(n_pa >= min_pa & n_pa > 0) .* repos;
Acttmp2 = cell(1, NCell); PromALL = cell(1, NCell); WidthALL = cell(1, NCell);
for i = 1:NCell
    l = find(Raster(i, :));
    Acttmp2{i} = l; PromALL{i} = n_pa(i, l); WidthALL{i} = ones(size(l));
end
th_detection = nan(1, NCell);
min_repos = numel(WinRest) / sampling_rate / 60;
info = struct('modele', modele, 'graine', graine, 'reechantillonnage_Hz', hz, ...
    'PA_par_min_repos', sum(n_pa(:, repos), 2) / max(eps, min_repos), ...
    'images_PA_par_min_repos', sum(Raster, 2) / max(eps, min_repos), 'gros_evenements', info_gros);
end


function sauver(f, n_pa, modele, graine, controle) %#ok<INUSD>
save(f, 'n_pa', 'modele', 'graine', 'controle', '-v7.3');
end
