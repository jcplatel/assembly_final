function etape1_fonctionnel(session, config)
%% =========================================================================
%  ÉTAPE 1 : prétraitement puis assemblées pour chaque K — AUCUNE COULEUR
%  -------------------------------------------------------------------------
%  Équivalent de la boucle de brainbow_analysis_2026_09_03 pour UN fichier.
%  Fichiers écrits dans le dossier de la session :
%     pretraitement.mat      sorties de preprocessing_7 (mêmes noms qu'avant :
%                            Tr1b, Raster, Race, TRace, DFF0, ...)
%     k9\assemblees.mat ...  sorties de SCE_clustering2 pour chaque K (IDX2,
%                            sCl, NClOK, assemblyortho, PCl_final, ...) et
%                            resume_K (ancien exportdata)
%     k9\analysis.xlsx, clusters.png, rasterall*.png, graphSCE.png
%  Un fichier déjà calculé avec les mêmes paramètres n'est pas refait.
%  Si le prétraitement est refait, tous les K sont refaits.
%  =========================================================================

opts = config.pre;
namefullold = avec_separateur(session.dossier);
if ~isfolder(namefullold), mkdir(namefullold); end
fichier_pre = fullfile(session.dossier, 'pretraitement.mat');

%% 1. Prétraitement (refait seulement s'il manque ou si opts a changé)
refaire_pre = config.forcer_pretraitement || ~isfile(fichier_pre) ...
    || ~memes_parametres(fichier_pre, 'opts_utilisees', opts);

if refaire_pre
    S = load(session.fichier_data, 'data_recording');
    data_recording = S.data_recording;
    clear S
    F             = data_recording.F;
    sampling_rate = data_recording.sampling_rate;
    speed         = data_recording.speed;
    path          = avec_separateur(session.dossier_suite2p);   % Fall.mat, calcium.tif, ...
    namefull      = namefullold;

    % P contient les sorties sous leurs noms habituels : P.Tr1b, P.Race, ...
    P = struct();
    [P.Tr1b, P.speedsm, P.Raster, P.SumAct, P.MAct, P.Race, P.RasterRace, P.WinRest, P.WinActive, ...
        P.TRace, P.Fzero, P.DFF0, P.th_detection, P.bad_frames, P.max_cells_allowed, ...
        P.n_transients_total, P.Acttmp2, P.PromALL, P.WidthALL, P.opts, P.diagnostic] = ...
        preprocessing_7(F, opts, speed, path, sampling_rate, namefull, session.graine);

    if isfield(data_recording, 'position') && ~isempty(data_recording.position)
        P.position = data_recording.position;
        [P.nb_lap, P.PosT] = lap_calculator(data_recording.position, speed);
    else
        P.position = [];
        P.nb_lap   = NaN;
        P.PosT     = [];
    end
    if isfield(data_recording, 'identifier')
        P.identifier = data_recording.identifier;
    else
        P.identifier = session.nom;
    end
    P.sampling_rate  = sampling_rate;
    P.speed          = speed;
    P.NCell          = size(P.Tr1b, 1);
    P.opts_utilisees = opts;     % pour savoir plus tard si ce prétraitement est à jour

    % Grosses matrices en single : fichier deux fois plus petit
    P.Tr1b  = single(P.Tr1b);
    P.DFF0  = single(P.DFF0);
    P.Fzero = single(P.Fzero);
    if isfile(fichier_pre), delete(fichier_pre); end
    save(fichier_pre, '-struct', 'P', '-v7.3');
    clear data_recording F
else
    P = load(fichier_pre);
end

%% 2. Anciens scripts lancés une fois par session (ex. raster_rastermap)
if config.executer_scripts_anciens && ~isempty(config.scripts_anciens_session)
    V = P;
    V.namefull = namefullold;
    V.path     = avec_separateur(session.dossier_suite2p);
    V.savefig  = double(config.faire_figures);
    [V.NCell, V.Nz] = size(P.Tr1b);
    executer_scripts_anciens(config.scripts_anciens_session, V);
end

%% 3. Assemblées pour chaque K
liste_K = config.liste_K;
for iK = 1:numel(liste_K)
    nanalysis = liste_K(iK);
    if isfield(config, 'progression_K') && ~isempty(config.progression_K)
        config.progression_K(iK - 1, numel(liste_K), nanalysis);   % avancement (interface)
    end
    NClini    = nanalysis;
    namefull  = avec_separateur(fullfile(session.dossier, sprintf('k%d', nanalysis)));
    fichier_K = fullfile(namefull, 'assemblees.mat');

    a_jour = ~refaire_pre && ~config.forcer_assemblees && isfile(fichier_K) ...
        && memes_parametres(fichier_K, 'config_asm', config.asm);
    if a_jour, continue; end
    if ~isfolder(namefull), mkdir(namefull); end

    A = SCE_clustering2(double(P.Race), P.TRace, P.Raster, double(P.DFF0), P.Acttmp2, P.PromALL, ...
        NClini, config.asm, opts.synchronous_frames, session.graine + nanalysis);
    A.config_asm = config.asm;
    A.resume_K = exportdata2(session, P, A, config, namefull);   % ligne de analysis.xlsx / analysisALL.xlsx
    if isfile(fichier_K), delete(fichier_K); end
    save(fichier_K, '-struct', 'A');

    if config.faire_figures
        figure_clusters(A, P.Race, namefull);                          % ancien SCE_clustering
        figure_raster_assemblees(A, P.Race, namefull);                 % ancien rastercolor
        figure_graphSCE(A, P.Tr1b, P.TRace, P.speed, namefull);        % ancien graphSCE
    end

    % Anciens scripts lancés pour chaque K (non relus)
    if config.executer_scripts_anciens && ~isempty(config.scripts_anciens_K)
        V = fusionner(P, A);
        V.namefull = namefull;
        V.assembly = A.assemblyortho;
        V.path     = avec_separateur(session.dossier_suite2p);
        V.savefig  = double(config.faire_figures);
        [V.NCell, V.Nz] = size(P.Tr1b);
        executer_scripts_anciens(config.scripts_anciens_K, V);
    end
end
end


function V = fusionner(V, A)
% Ajoute les champs de A à V (A est prioritaire).
noms = fieldnames(A);
for k = 1:numel(noms)
    V.(noms{k}) = A.(noms{k});
end
end
