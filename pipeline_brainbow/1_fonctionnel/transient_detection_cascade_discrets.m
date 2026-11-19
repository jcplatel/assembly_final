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
%  opts : cascade_modele, cascade_processus (défaut 24), cascade_graine (1).
%  =========================================================================
if isfield(opts, 'cascade_modele'), modele = opts.cascade_modele;
else, modele = 'GC8_EXC_10Hz_smoothing150ms_high_noise'; end
if isfield(opts, 'cascade_processus'), n_proc = opts.cascade_processus; else, n_proc = 24; end
if isfield(opts, 'cascade_graine'), graine = opts.cascade_graine; else, graine = 1; end
o = opts; o.cascade_modele = modele;
sp = cascade_spikes(DFF0, o, fullfile(path, ['cascade_' modele '.mat']));

% PA discrets (gardés : relus si même spike_prob, même modèle, même graine)
fc = fullfile(path, ['cascade_discrets_' modele '.mat']);
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

[NCell, Nz] = size(n_pa);
repos = false(1, Nz); repos(WinRest) = true;
Raster = double(n_pa > 0) .* repos;
Acttmp2 = cell(1, NCell); PromALL = cell(1, NCell); WidthALL = cell(1, NCell);
for i = 1:NCell
    l = find(Raster(i, :));
    Acttmp2{i} = l; PromALL{i} = n_pa(i, l); WidthALL{i} = ones(size(l));
end
th_detection = nan(1, NCell);
min_repos = numel(WinRest) / sampling_rate / 60;
info = struct('modele', modele, 'graine', graine, ...
    'PA_par_min_repos', sum(n_pa(:, repos), 2) / max(eps, min_repos), ...
    'images_PA_par_min_repos', sum(Raster, 2) / max(eps, min_repos));
end


function sauver(f, n_pa, modele, graine, controle) %#ok<INUSD>
save(f, 'n_pa', 'modele', 'graine', 'controle', '-v7.3');
end
