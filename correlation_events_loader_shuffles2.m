clear; clc; close all

%% =========================================================
% CORRELATIONS / COACTIVATIONS PAR COULEUR — AVEC SHUFFLES
%
% Base : correlation_events_loader.m
% Cette version corrige l'extraction de l'animal : Animal_ID n'est
% JAMAIS renseigné dans data_recording sur cette base. L'animal est
% donc extrait directement de "identifier", qui suit la convention :
%   <préfixe>_<6 chiffres animal>_<date>_<planeX>[...]
%   ex : "P2M_411582_230320_plane0_K4"      -> animal = "411582"
%        "P2M_446178_10062025_plane1"       -> animal = "446178"
% (la date peut faire 6 OU 8 chiffres selon les sessions ; ça n'a
% aucune importance puisqu'on prend directement le 2e token, pas un
% bloc de 6 chiffres cherché dans toute la chaîne).
%% =========================================================

%% 1. CHEMINS ET CHARGEMENT DE LA LISTE
if ispc
    PathSave = "E:\Data\Aurelie\analysis\Final\correlation\";
    database_file = "E:\Data\Aurelie\data\nocues_data_recording_all.mat";
elseif ismac
    PathSave = "/Users/platel/Desktop/exp/brainbow_hippocampus/analysis/nocues";
    database_file = "/Users/platel/Desktop/exp/brainbow_hippocampus/analysis/final_analysis/bestnocuesold.mat";
else
    error('Système non configuré : définir PathSave et database_file.');
end
if ~exist(PathSave,'dir'), mkdir(PathSave); end

figures_shuffle_dir = fullfile(PathSave,'figures_shuffle_sessions');
if ~exist(figures_shuffle_dir,'dir'), mkdir(figures_shuffle_dir); end

DB = load(database_file);
if isfield(DB,'matfile')
    recording_files = DB.matfile;
elseif isfield(DB,'files')
    recording_files = DB.files;
elseif isfield(DB,'filelist')
    recording_files = DB.filelist;
else
    error('La base %s ne contient pas la variable matfile/files/filelist.',database_file);
end
if isstring(recording_files), recording_files = cellstr(recording_files); end
if ischar(recording_files), recording_files = cellstr(recording_files); end

%% 2. PARAMETRES
analysis = struct();
analysis.pixel_size_um = 1.37;
analysis.radii_um = [100 200 300];
analysis.min_same_pairs = 5;
analysis.min_diff_pairs = 5;
analysis.event_source = 'Raster';      % 'Raster' ou 'Race'
analysis.event_width_frames = 2;       % 1 = peak strict ; 2 = tolérance 1 frame
analysis.min_events_per_cell = 1;
analysis.max_exact_animals = 20;
analysis.n_signflip_mc = 1000;

% --- Shuffles ---
analysis.n_perm_color_session = 1000;      % shuffle couleurs, par session (rapide)
analysis.n_perm_temporal_session = 500;    % shuffle temporel, par session (coûteux)
analysis.n_perm_color_hierarchical = 1000; % shuffle couleurs, niveau population (test principal)
analysis.temporal_min_shift_sec = 5;       % décalage minimal (évite un shift quasi nul)
analysis.session_figure_scope = 'global';  % portée illustrée dans les histogrammes par session
analysis.save_session_shuffle_figures = true;
analysis.show_session_shuffle_figures = false;

analysis.random_seed = 20260916;
analysis.stop_on_error = false;
rng(analysis.random_seed,'twister');

metrics = {'phi','spearman','pearson','jaccard'};
scopes = [{'global'}, arrayfun(@(x)sprintf('R%d',x),analysis.radii_um,'UniformOutput',false)];
keep_null_scopes = {analysis.session_figure_scope};

files_to_process = 1:numel(recording_files); % ex. [1:3 23:30]

%% 3. JOURNAL
fig_log = uifigure('Name','Journal de Traitement', ...
    'Position',[500 300 500 400],'Visible','off','Tag','CorrelationProcessingLog');
txt_log = uitextarea(fig_log,'Position',[20 20 460 360], ...
    'Editable','off','FontName','Consolas');
log_file = fullfile(PathSave,sprintf('log_analyse_%s.txt',datestr(now,'yy_mm_dd_HH_MM')));
historique = {};
historique = ecrire_log(txt_log,log_file,historique, ...
    sprintf('=== DÉBUT DE L''ANALYSE : %s ===',datestr(now)));
temps_global = tic;

results_all = struct([]);
summary_sessions = table();

%% 4. BOUCLE SUR LES FICHIERS
for file_num = files_to_process
    temps_fichier = tic;
    filename = string(recording_files{file_num});

    figs = findall(groot,'Type','figure');
    if ~isempty(figs), close(figs); end

    try
        historique = ecrire_log(txt_log,log_file,historique, ...
            sprintf('Traitement %d/%d : %s',file_num,numel(recording_files),filename));

        %% 4.1 Charge l'entrée data_recording et détermine l'animal
        [Session,recording_file] = load_recording(filename);
        identifier = string(get_required_field(Session,'identifier'));
        pathload = get_recording_path(Session,recording_file);

        [animal_id_extracted,animal_id_ok] = extract_animal_id(identifier);

        if isfield(Session,'Animal_ID') && ~isempty(Session.Animal_ID) && string(Session.Animal_ID) ~= ""
            % Champ Animal_ID réellement renseigné : on lui fait confiance,
            % mais on signale toute incohérence avec l'identifiant.
            animal_id = string(Session.Animal_ID);
            if animal_id_ok && animal_id ~= animal_id_extracted
                historique = ecrire_log(txt_log,log_file,historique, ...
                    sprintf(['ATTENTION : Animal_ID champ (''%s'') != animal ' ...
                    'extrait de l''identifiant (''%s'') pour %s.'], ...
                    animal_id,animal_id_extracted,identifier));
            end
        elseif animal_id_ok
            % Cas normal sur cette base : Animal_ID absent, on extrait
            % le numéro à 6 chiffres depuis l'identifiant (2e token).
            animal_id = animal_id_extracted;
        else
            % Échec total d'extraction : ne JAMAIS retomber silencieusement
            % sur l'identifiant complet (ça casserait le regroupement par
            % animal). On le signale bruyamment à la place.
            animal_id = "UNKNOWN_" + identifier;
            historique = ecrire_log(txt_log,log_file,historique, ...
                sprintf(['ERREUR EXTRACTION ANIMAL : impossible de déduire ' ...
                'un numéro à 6 chiffres depuis ''%s''. Session marquée %s, ' ...
                'à corriger manuellement avant l''analyse finale.'], ...
                identifier,animal_id));
        end

        %% 4.2 Charge Fall.mat : stat + iscell, et F en secours
        fall_file = fullfile(pathload,'Fall.mat');
        if ~isfile(fall_file)
            error('Fall.mat introuvable : %s',fall_file);
        end
        Fall = load(fall_file,'stat','iscell','F');
        if ~isfield(Fall,'stat') || ~isfield(Fall,'iscell')
            error('Fall.mat ne contient pas stat et iscell : %s',fall_file);
        end

        stat = Fall.stat;
        iscell = Fall.iscell;
        suite2p_mask = iscell(:,1)==1;
        n_rois_total = numel(stat);
        n_cells_suite2p = nnz(suite2p_mask);

        if size(iscell,1)~=n_rois_total
            error('Fall.mat incohérent : %d stat mais %d lignes iscell.', ...
                n_rois_total,size(iscell,1));
        end

        %% 4.3 F : applique iscell une seule fois
        if isfield(Session,'F') && ~isempty(Session.F)
            F = Session.F;
        elseif isfield(Fall,'F') && ~isempty(Fall.F)
            F = Fall.F;
        else
            error('F absent de data_recording et de Fall.mat.');
        end
        F = filter_suite2p_rows(F,suite2p_mask,n_cells_suite2p,'F');

        %% 4.4 Centroïdes des seules ROI iscell==1
        gcamp_centroid = centroids_from_suite2p(stat,suite2p_mask);
        if size(gcamp_centroid,1)~=size(F,1)
            error('Centroïdes (%d) et F filtré (%d) ne correspondent pas.', ...
                size(gcamp_centroid,1),size(F,1));
        end

        %% 4.5 Variables d'enregistrement
        speed = get_required_field(Session,'speed');
        sampling_rate = get_required_field(Session,'sampling_rate');

        namefull = fullfile(PathSave,char(identifier + "_" + datestr(now,'yy_mm_dd_HH_MM_SS')));
        if ~exist(namefull,'dir'), mkdir(namefull); end

        opts = struct('MinPeakDistancesce',5,'MinPeakDistance',3, ...
            'threshold_peak',3.09,'synchronous_frames',2, ...
            'sce_n_cells_threshold',10,'SG_window',7, ...
            'motion_correction',true,'colorsubstraction',false);

        [Tr1b,speedsm,Raster,SumAct,MAct,Race,RasterRace,WinRest,WinActive, ...
            TRace,Fzero,DFF0,th_detection,bad_frames,max_cells_allowed, ...
            n_transients_total,opts] = ...
            preprocessing_6(F,opts,speed,pathload,sampling_rate,namefull); %#ok

        %% 4.6 Couleurs, alignées sur le même masque iscell
        if isfield(Session,'Colorcell_Version') && ...
                ~isempty(Session.Colorcell_Version) && ...
                ~strcmp(Session.Colorcell_Version,'Non trouvé')
            color_file = fullfile(pathload,Session.Colorcell_Version);
        else
            color_file = fullfile(pathload,'colorcell_Maxwell.mat');
        end
        if ~isfile(color_file), error('Fichier couleur introuvable : %s',color_file); end
        C = load(color_file,'colorcell');
        if ~isfield(C,'colorcell'), error('colorcell absent de %s.',color_file); end
        colorcell = filter_suite2p_rows(C.colorcell,suite2p_mask,n_cells_suite2p,'colorcell');
        colorcell = double(colorcell(:));

        %% 4.7 Repos, détection et binarisation
        tCorr = tic;

        Fluo = double(Tr1b(:,WinRest));
        switch lower(analysis.event_source)
            case 'raster', EventRaw = Raster;
            case 'race',   EventRaw = Race;
            otherwise, error('event_source doit être Raster ou Race.');
        end
        EventRaw = select_time_window(EventRaw,WinRest,size(Tr1b,2),'EventRaw');

        EventPeaks = detected_runs_to_peaks(EventRaw>0,Fluo);
        EventRaster = widen_binary_events(EventPeaks,analysis.event_width_frames);
        n_events = sum(EventPeaks,2);

        %% 4.8 Cellules valides finales
        valid = isfinite(colorcell) & colorcell>0 & ...
            n_events>=analysis.min_events_per_cell & ...
            std(EventRaster,0,2)>0 & ...
            all(isfinite(gcamp_centroid),2);
        if nnz(valid)<3
            error('Moins de trois cellules valides après filtrage.');
        end

        FluoV = Fluo(valid,:);
        EventV = EventRaster(valid,:);
        PeaksV = EventPeaks(valid,:);
        colorsV = colorcell(valid);
        eventsV = n_events(valid);
        coordsV_um = gcamp_centroid(valid,:)*analysis.pixel_size_um;
        DistMat = squareform(pdist(coordsV_um));
        Nvalid = size(FluoV,1);

        %% 4.9 Similarités observées
        M = compute_similarity_matrices(FluoV,EventV);

        out = struct();
        for m = 1:numel(metrics)
            met = metrics{m};
            out.(met) = summarize_same_diff_radii(M.(met),colorsV,DistMat, ...
                analysis.radii_um,analysis.min_same_pairs,analysis.min_diff_pairs);
        end
        runtime_corr_s = toc(tCorr);

        %% 4.10 Masques de portée (paires candidates, indépendants des couleurs)
        poolMask = build_pool_masks(Nvalid,DistMat,analysis.radii_um);

        %% 4.11 Shuffle des COULEURS (test principal, par session)
        tColor = tic;
        colorRes = color_shuffle_all_metrics_scopes(M,colorsV,poolMask,metrics, ...
            analysis.n_perm_color_session,analysis.min_same_pairs, ...
            analysis.min_diff_pairs,keep_null_scopes);
        runtime_color_s = toc(tColor);

        %% 4.12 Shuffle TEMPOREL (contrôle secondaire, par session)
        tTime = tic;
        min_shift_frames = max(1,round(analysis.temporal_min_shift_sec*sampling_rate));
        temporalRes = temporal_shuffle_all_metrics_scopes(FluoV,PeaksV,colorsV, ...
            poolMask,M,metrics,analysis.n_perm_temporal_session, ...
            analysis.min_same_pairs,analysis.min_diff_pairs, ...
            analysis.event_width_frames,min_shift_frames,keep_null_scopes);
        runtime_time_s = toc(tTime);

        %% 4.13 Figure par session : histogrammes des deux shuffles
        tFig = tic;
        if analysis.save_session_shuffle_figures
            scope_fig = analysis.session_figure_scope;
            colorResScope = colorRes.(scope_fig);
            temporalResScope = temporalRes.(scope_fig);
            fig_path = fullfile(figures_shuffle_dir, ...
                sprintf('Shuffle_%s_%s.png',identifier,scope_fig));
            plot_session_shuffle_histograms(colorResScope,temporalResScope, ...
                metrics,scope_fig,char(identifier),fig_path, ...
                analysis.show_session_shuffle_figures);
        end
        runtime_fig_s = toc(tFig);

        %% 4.14 Assemblage de la sortie session
        out.event_peaks = PeaksV;
        out.event_raster = EventV;
        out.event_count = eventsV;
        out.event_rate_hz = eventsV/(size(FluoV,2)/sampling_rate);
        out.colors = colorsV;
        out.centroids_pixels = gcamp_centroid(valid,:);
        out.centroids_um = coordsV_um;
        out.distance_matrix_um = DistMat;
        out.options = analysis;
        out.source_recording = recording_file;
        out.source_fall = fall_file;
        out.animal_id_ok = animal_id_ok;

        out.M = M;
        out.pool = poolMask;
        out.color_shuffle = colorRes;
        out.temporal_shuffle = temporalRes;

        out.runtime.corr_s = runtime_corr_s;
        out.runtime.color_shuffle_s = runtime_color_s;
        out.runtime.temporal_shuffle_s = runtime_time_s;
        out.runtime.figures_s = runtime_fig_s;

        %% 4.15 Résumé session (table)
        row = struct('Animal_ID',animal_id,'identifier',identifier, ...
            'file_num',file_num,'n_cells',Nvalid, ...
            'n_timepoints_rest',size(FluoV,2), ...
            'median_events_per_cell',median(eventsV,'omitnan'), ...
            'median_event_rate_hz',median(out.event_rate_hz,'omitnan'), ...
            'animal_id_ok',animal_id_ok, ...
            'runtime_corr_s',runtime_corr_s, ...
            'runtime_color_shuffle_s',runtime_color_s, ...
            'runtime_temporal_shuffle_s',runtime_time_s, ...
            'runtime_figures_s',runtime_fig_s);

        for m = 1:numel(metrics)
            met = metrics{m};
            row.(['mean_same_' met '_global']) = out.(met).global.mean_same;
            row.(['mean_diff_' met '_global']) = out.(met).global.mean_diff;
            row.(['delta_' met '_global']) = out.(met).global.delta;
            row.(['n_same_' met '_global']) = out.(met).global.n_same;
            row.(['n_diff_' met '_global']) = out.(met).global.n_diff;
            row.(['p_color_' met '_global']) = colorRes.global.(met).p_value;
            row.(['p_temporal_' met '_global']) = temporalRes.global.(met).p_value;

            for r = 1:numel(analysis.radii_um)
                tag = sprintf('R%d',analysis.radii_um(r));
                row.(['mean_same_' met '_' tag]) = out.(met).radius.mean_same(r);
                row.(['mean_diff_' met '_' tag]) = out.(met).radius.mean_diff(r);
                row.(['delta_' met '_' tag]) = out.(met).radius.delta(r);
                row.(['n_same_' met '_' tag]) = out.(met).radius.n_same(r);
                row.(['n_diff_' met '_' tag]) = out.(met).radius.n_diff(r);
                row.(['n_cells_same_' met '_' tag]) = out.(met).radius.n_cells_same(r);
                row.(['n_cells_diff_' met '_' tag]) = out.(met).radius.n_cells_diff(r);
                row.(['enough_' met '_' tag]) = out.(met).radius.enough(r);
                row.(['p_color_' met '_' tag]) = colorRes.(tag).(met).p_value;
                row.(['p_temporal_' met '_' tag]) = temporalRes.(tag).(met).p_value;
            end
        end

        summary_sessions = [summary_sessions;struct2table(row)]; %#ok
        kresult = numel(results_all)+1;
        results_all(kresult).identifier = identifier;
        results_all(kresult).filename = filename;
        results_all(kresult).out = out;
        results_all(kresult).summary = row;

        save(fullfile(PathSave,'results_all_backup.mat'), ...
            'results_all','summary_sessions','analysis','-v7.3');

        historique = ecrire_log(txt_log,log_file,historique, ...
            sprintf(['OK %s (animal %s) : corr=%.1fs, shuffle couleur=%.1fs, ' ...
            'shuffle temporel=%.1fs, figures=%.1fs, total=%.2f min.'], ...
            identifier,animal_id,runtime_corr_s,runtime_color_s,runtime_time_s, ...
            runtime_fig_s,toc(temps_fichier)/60));

    catch ME
        historique = ecrire_log(txt_log,log_file,historique, ...
            sprintf('ERREUR fichier %d (%s) : %s',file_num,filename,ME.message));
        if analysis.stop_on_error, rethrow(ME); end
    end
end

if isempty(summary_sessions), error('Aucune session valide n''a été analysée.'); end

%% Vérification post-boucle : le regroupement a-t-il fonctionné ?
n_sessions_total = height(summary_sessions);
n_animals_final = numel(unique(summary_sessions.Animal_ID));
n_unknown = sum(startsWith(summary_sessions.Animal_ID,"UNKNOWN_"));

fprintf('\n=== VÉRIFICATION DU REGROUPEMENT ANIMAL ===\n');
fprintf('Sessions analysées         : %d\n',n_sessions_total);
fprintf('Animaux identifiés         : %d\n',n_animals_final);
fprintf('Sessions en échec d''extraction (UNKNOWN_*) : %d\n',n_unknown);

if n_animals_final >= n_sessions_total*0.8
    fprintf(2,['\nALERTE PERSISTANTE : le regroupement par animal ne ' ...
        'fonctionne toujours pas (%d animaux pour %d sessions). ' ...
        'Vérifier extract_animal_id() et le format des identifiants.\n'], ...
        n_animals_final,n_sessions_total);
end

if n_unknown > 0
    fprintf(2,'\nSessions à corriger manuellement :\n');
    disp(summary_sessions(startsWith(summary_sessions.Animal_ID,"UNKNOWN_"), ...
        {'file_num','identifier','Animal_ID'}))
end

%% 5. AGRÉGATION NON PONDEREE PAR ANIMAL (deltas observés)
animal_ids = unique(summary_sessions.Animal_ID,'stable');
summary_animals_true = table();

for a = 1:numel(animal_ids)
    id = animal_ids(a);
    S = summary_sessions(summary_sessions.Animal_ID==id,:);
    row = struct('Animal_ID',id,'Nb_Plans_Averaged',height(S));
    for m = 1:numel(metrics)
        met = metrics{m};
        for q = 1:numel(scopes)
            scope = scopes{q};
            row.(['mean_same_' met '_' scope]) = mean(S.(['mean_same_' met '_' scope]),'omitnan');
            row.(['mean_diff_' met '_' scope]) = mean(S.(['mean_diff_' met '_' scope]),'omitnan');
            row.(['delta_' met '_' scope]) = mean(S.(['delta_' met '_' scope]),'omitnan');
            row.(['N_Plans_' met '_' scope]) = sum(isfinite(S.(['delta_' met '_' scope])));
        end
    end
    summary_animals_true = [summary_animals_true;struct2table(row)]; %#ok
end

disp(summary_animals_true)

%% 6. TESTS AU NIVEAU ANIMAL
final_stats = table();
temps_hier = tic;

for m = 1:numel(metrics)
    met = metrics{m};
    for q = 1:numel(scopes)
        scope = scopes{q};
        D = summary_animals_true.(['delta_' met '_' scope]);

        st = signflip_test(D,'right',analysis.max_exact_animals,analysis.n_signflip_mc);
        try
            p_wilcoxon = signrank(D,0,'tail','right');
        catch
            p_wilcoxon = NaN;
        end

        [Tobs,~,p_hier] = hierarchical_color_shuffle(results_all,met,scope, ...
            animal_ids,analysis.n_perm_color_hierarchical, ...
            analysis.min_same_pairs,analysis.min_diff_pairs);

        newrow = table(string(met),string(scope),st.N,st.mean_delta, ...
            st.median_delta,st.p_value,p_wilcoxon,string(st.method), ...
            Tobs,p_hier, ...
            'VariableNames',{'Metric','Scope','N_animals','Mean_delta', ...
            'Median_delta','P_signflip','P_wilcoxon','Signflip_method', ...
            'T_hierarchical','P_hierarchical_color'});
        final_stats = [final_stats;newrow]; %#ok
    end
end
runtime_hier_s = toc(temps_hier);
historique = ecrire_log(txt_log,log_file,historique, ...
    sprintf('Test hiérarchique de shuffle couleur terminé en %.2f min.',runtime_hier_s/60));

final_stats.P_Holm_signflip = nan(height(final_stats),1);
final_stats.P_Holm_hierarchical = nan(height(final_stats),1);
for m = 1:numel(metrics)
    idx = final_stats.Metric==string(metrics{m});
    final_stats.P_Holm_signflip(idx) = holm_adjust(final_stats.P_signflip(idx));
    final_stats.P_Holm_hierarchical(idx) = holm_adjust(final_stats.P_hierarchical_color(idx));
end

disp(final_stats)

%% 7. FIGURES ANIMAUX
for m = 1:numel(metrics)
    met = metrics{m};
    idx = final_stats.Metric==string(met) & final_stats.Scope=="global";
    plot_paired_animals(summary_animals_true.(['mean_same_' met '_global']), ...
        summary_animals_true.(['mean_diff_' met '_global']),metric_label(met), ...
        sprintf('%s global : p_{hier-couleur}=%.4g, p_{sign-flip}=%.4g, N=%d', ...
        upper(met),final_stats.P_hierarchical_color(idx), ...
        final_stats.P_signflip(idx),final_stats.N_animals(idx)), ...
        fullfile(PathSave,sprintf('Figure_%s_global.png',met)));

    D = nan(height(summary_animals_true),numel(analysis.radii_um));
    for r = 1:numel(analysis.radii_um)
        D(:,r) = summary_animals_true.(['delta_' met sprintf('_R%d',analysis.radii_um(r))]);
    end
    plot_radius_deltas(D,analysis.radii_um,metric_label(met), ...
        sprintf('%s : paires dans le rayon indiqué',upper(met)), ...
        fullfile(PathSave,sprintf('Figure_%s_radii.png',met)));
end

%% 8. SAUVEGARDE FINALE
save(fullfile(PathSave,'results_all_final_events_shuffles.mat'), ...
    'results_all','summary_sessions','summary_animals_true','final_stats','analysis','-v7.3');
writetable(summary_sessions,fullfile(PathSave,'summary_sessions_events_shuffles.csv'));
writetable(summary_animals_true,fullfile(PathSave,'summary_animals_events_shuffles.csv'));
writetable(final_stats,fullfile(PathSave,'final_stats_events_shuffles.csv'));
historique = ecrire_log(txt_log,log_file,historique, ...
    sprintf('=== FIN : durée totale %.2f min ===',toc(temps_global)/60));

%% =========================================================
% FONCTIONS LOCALES — CHARGEMENT / PRÉTRAITEMENT
%% =========================================================
function [Session,recording_file] = load_recording(entry)
entry = string(entry);
if isfolder(entry)
    candidates = [fullfile(entry,'data_recording.mat'),fullfile(entry,'brainbow.mat')];
    recording_file = "";
    for k=1:numel(candidates)
        if isfile(candidates(k)), recording_file=candidates(k); break; end
    end
    if recording_file=="", error('Aucun data_recording.mat/brainbow.mat dans %s.',entry); end
elseif isfile(entry)
    recording_file = entry;
else
    error('Entrée introuvable : %s',entry);
end

S = load(recording_file);
if isfield(S,'data_recording') && isstruct(S.data_recording)
    Session = S.data_recording;
elseif isfield(S,'Session') && isstruct(S.Session)
    Session = S.Session;
else
    Session = S;
end
end

function pathload = get_recording_path(Session,recording_file)
if isfield(Session,'path') && ~isempty(Session.path)
    pathload = char(string(Session.path));
elseif isfield(Session,'DossierSource') && ~isempty(Session.DossierSource)
    pathload = char(string(Session.DossierSource));
else
    pathload = char(fileparts(recording_file));
end
if ~isfolder(pathload), error('Dossier path introuvable : %s',pathload); end
end

function value = get_required_field(S,name)
if ~isfield(S,name) || isempty(S.(name)), error('Champ requis absent : %s',name); end
value = S.(name);
end

function [animal_id,ok] = extract_animal_id(identifier)
% "P2M_411582_230320_plane0_K4" -> "411582"
% "P2M_446178_10062025_plane1"  -> "446178"  (date à 8 chiffres, sans effet)
% Règle : le 2e token séparé par "_" est le numéro d'animal à 6 chiffres.
ok = false;
animal_id = "";
identifier = string(identifier);
if identifier=="", return; end

parts = split(identifier,'_');

if numel(parts) >= 2 && ~isempty(regexp(parts(2),'^\d{6}$','once'))
    animal_id = parts(2);
    ok = true;
    return
end

% Repli : n'importe quel token de exactement 6 chiffres, où qu'il soit
for i = 1:numel(parts)
    if ~isempty(regexp(parts(i),'^\d{6}$','once'))
        animal_id = parts(i);
        ok = true;
        return
    end
end
end

function X = filter_suite2p_rows(X,mask,n_cells,name)
if isvector(X) && ~ischar(X) && ~isstring(X)
    was_row = isrow(X);
    X = X(:);
    if numel(X)==numel(mask)
        X = X(mask);
    elseif numel(X)~=n_cells
        error('%s : %d éléments; attendu %d ROI totales ou %d cellules iscell.', ...
            name,numel(X),numel(mask),n_cells);
    end
    if was_row, X=X'; end
elseif size(X,1)==numel(mask)
    X = X(mask,:);
elseif size(X,1)~=n_cells
    error('%s : %d lignes; attendu %d ROI totales ou %d cellules iscell.', ...
        name,size(X,1),numel(mask),n_cells);
end
end

function centroids = centroids_from_suite2p(stat,mask)
idx=find(mask); centroids=nan(numel(idx),2);
for k=1:numel(idx)
    if iscell(stat), s=stat{idx(k)}; else, s=stat(idx(k)); end
    if isfield(s,'xpix') && isfield(s,'ypix') && ~isempty(s.xpix)
        x=double(s.xpix(:))+1; y=double(s.ypix(:))+1;
    elseif isfield(s,'xext') && isfield(s,'yext') && ~isempty(s.xext)
        x=double(s.xext(:))+1; y=double(s.yext(:))+1;
    else
        error('ROI Suite2p %d sans xpix/ypix ni xext/yext.',idx(k));
    end
    centroids(k,:)=[mean(x,'omitnan') mean(y,'omitnan')];
end
end

function X = select_time_window(X,WinRest,n_total,name)
if size(X,2)==n_total, X=X(:,WinRest);
elseif size(X,2)~=numel(WinRest)
    error('%s : nombre de frames incompatible.',name);
end
end

function Peaks = detected_runs_to_peaks(B,Signal)
B=logical(B); Peaks=false(size(B));
if ~isequal(size(B),size(Signal)), error('Raster et Fluo doivent avoir la même taille.'); end
for i=1:size(B,1)
    starts=find(B(i,:) & [true ~B(i,1:end-1)]);
    stops=find(B(i,:) & [~B(i,2:end) true]);
    for e=1:numel(starts)
        frames=starts(e):stops(e);
        [~,loc]=max(Signal(i,frames));
        Peaks(i,frames(loc))=true;
    end
end
end

function B = widen_binary_events(Peaks,width_frames)
width_frames=max(1,round(width_frames));
B=filter(ones(1,width_frames),1,double(Peaks),[],2)>0;
end

function M = compute_similarity_matrices(Fluo,Events)
M.pearson=corr(Fluo','Type','Pearson','Rows','pairwise');
M.spearman=corr(Fluo','Type','Spearman','Rows','pairwise');
M.phi=corr(double(Events'),'Type','Pearson','Rows','pairwise');
B=double(Events); intersection=B*B'; activity=sum(B,2);
union=activity+activity'-intersection;
M.jaccard=intersection./union;
M.jaccard(union==0)=NaN;
M.jaccard(1:size(B,1)+1:end)=1;
end

function out = summarize_same_diff_radii(Sim,labels,Dist,radii,min_same,min_diff)
N=size(Sim,1); tri=triu(true(N),1); same=labels(:)==labels(:)';
same0=tri & same & isfinite(Sim); diff0=tri & ~same & isfinite(Sim);
vs=Sim(same0); vd=Sim(diff0);
out.global.mean_same=mean(vs,'omitnan');
out.global.mean_diff=mean(vd,'omitnan');
out.global.delta=out.global.mean_same-out.global.mean_diff;
out.global.n_same=numel(vs); out.global.n_diff=numel(vd);

nr=numel(radii);
out.radius.values_um=radii;
out.radius.mean_same=nan(1,nr); out.radius.mean_diff=nan(1,nr);
out.radius.delta=nan(1,nr); out.radius.n_same=zeros(1,nr);
out.radius.n_diff=zeros(1,nr); out.radius.n_cells_same=zeros(1,nr);
out.radius.n_cells_diff=zeros(1,nr); out.radius.enough=false(1,nr);

for r=1:nr
    within=Dist<=radii(r);
    ms=same0 & within; md=diff0 & within;
    vs=Sim(ms); vd=Sim(md);
    out.radius.n_same(r)=numel(vs); out.radius.n_diff(r)=numel(vd);
    [is,js]=find(ms); [id,jd]=find(md);
    out.radius.n_cells_same(r)=numel(unique([is;js]));
    out.radius.n_cells_diff(r)=numel(unique([id;jd]));
    out.radius.enough(r)=numel(vs)>=min_same && numel(vd)>=min_diff;
    if out.radius.enough(r)
        out.radius.mean_same(r)=mean(vs,'omitnan');
        out.radius.mean_diff(r)=mean(vd,'omitnan');
        out.radius.delta(r)=out.radius.mean_same(r)-out.radius.mean_diff(r);
    end
end
end

function st = signflip_test(D,tail,max_exact,nmc)
D=D(isfinite(D)); D=D(:); st.N=numel(D);
st.mean_delta=mean(D,'omitnan'); st.median_delta=median(D,'omitnan');
if isempty(D), st.p_value=NaN; st.method='none'; return; end
if numel(D)<=max_exact
    S=dec2bin(0:(2^numel(D)-1))-'0'; S(S==0)=-1;
    T=mean(S.*D',2); st.method='exact';
    mc=false;
else
    S=2*(rand(nmc,numel(D))>=.5)-1;
    T=mean(S.*D',2); st.method='Monte-Carlo';
    mc=true;
end
switch lower(tail)
    case 'right', ne=sum(T>=st.mean_delta);
    case 'left', ne=sum(T<=st.mean_delta);
    otherwise, ne=sum(abs(T)>=abs(st.mean_delta));
end
if mc, st.p_value=(ne+1)/(numel(T)+1); else, st.p_value=ne/numel(T); end
end

function padj = holm_adjust(p)
padj=nan(size(p)); ok=isfinite(p); pv=p(ok); [ps,ord]=sort(pv); m=numel(ps);
adj=zeros(m,1); running=0;
for i=1:m
    running=max(running,(m-i+1)*ps(i)); adj(i)=min(1,running);
end
tmp=nan(m,1); tmp(ord)=adj; padj(ok)=tmp;
end

function txt = metric_label(metric)
switch lower(metric)
    case 'phi', txt='Coefficient phi (événements binaires)';
    case 'spearman', txt='Corrélation de Spearman (traces continues)';
    case 'pearson', txt='Corrélation de Pearson (traces continues)';
    case 'jaccard', txt='Indice de Jaccard (événements binaires)';
end
end

function plot_paired_animals(a,b,ylab,tit,pathout)
ok=isfinite(a)&isfinite(b); a=a(ok); b=b(ok);
f=figure('Position',[100 100 520 500],'Color','w','Visible','off'); hold on
for i=1:numel(a)
    plot([1 2],[a(i) b(i)],'-','Color',[.72 .72 .72]);
    plot(1,a(i),'o','MarkerFaceColor',[.15 .45 .75],'Color',[.15 .45 .75]);
    plot(2,b(i),'o','MarkerFaceColor',[.85 .35 .20],'Color',[.85 .35 .20]);
end
mu=[mean(a,'omitnan') mean(b,'omitnan')];
sem=[std(a,0,'omitnan') std(b,0,'omitnan')]/sqrt(numel(a));
errorbar([1 2],mu,sem,'k-','LineWidth',2,'CapSize',8);
xlim([.5 2.5]); xticks([1 2]); xticklabels({'Même couleur','Différente'});
ylabel(ylab); title(tit); set(gca,'Box','off','TickDir','out','FontName','Arial');
exportgraphics(f,pathout,'Resolution',300);
close(f);
end

function plot_radius_deltas(D,radii,ylab,tit,pathout)
f=figure('Position',[100 100 620 500],'Color','w','Visible','off'); hold on
for i=1:size(D,1), plot(radii,D(i,:),'-o','Color',[.75 .75 .75]); end
n=sum(isfinite(D),1); mu=mean(D,1,'omitnan');
sem=std(D,0,1,'omitnan')./sqrt(n);
errorbar(radii,mu,sem,'k-o','LineWidth',2,'MarkerFaceColor','k','CapSize',8);
yline(0,'k--'); xticks(radii); xlabel('Distance maximale (µm)');
ylabel(['Delta ' ylab ' (même - différente)']); title(tit);
set(gca,'Box','off','TickDir','out','FontName','Arial');
exportgraphics(f,pathout,'Resolution',300);
close(f);
end

function historique = ecrire_log(txt_log,log_file,historique,message)
timestamped=sprintf('[%s] %s',datestr(now,'HH:MM:SS'),message);
historique{end+1}=timestamped;
if isvalid(txt_log)
    txt_log.Value=historique(max(1,end-200):end);
    drawnow limitrate
end
fid=fopen(log_file,'a');
if fid~=-1, fprintf(fid,'%s\n',timestamped); fclose(fid); end
fprintf('%s\n',timestamped);
end

%% =========================================================
% FONCTIONS LOCALES — SHUFFLES
%% =========================================================
function poolMask = build_pool_masks(N,Dist,radii)
tri = triu(true(N),1);
poolMask.global = tri;
for r = 1:numel(radii)
    tag = sprintf('R%d',radii(r));
    poolMask.(tag) = tri & (Dist<=radii(r));
end
end

function res = color_shuffle_all_metrics_scopes(M,labels,poolMask,metrics, ...
    Nperm,min_same,min_diff,keep_null_scopes)
scopeNames = fieldnames(poolMask);
labels = labels(:);
N = numel(labels);

finiteMask = struct();
for m = 1:numel(metrics)
    finiteMask.(metrics{m}) = isfinite(M.(metrics{m}));
end

res = struct();
sameLbl = labels==labels';
for s = 1:numel(scopeNames)
    scope = scopeNames{s};
    sameMask = poolMask.(scope) & sameLbl;
    diffMask = poolMask.(scope) & ~sameLbl;
    for m = 1:numel(metrics)
        met = metrics{m};
        fm = finiteMask.(met);
        vs = M.(met)(sameMask & fm);
        vd = M.(met)(diffMask & fm);
        r = struct('n_same',numel(vs),'n_diff',numel(vd));
        r.enough = r.n_same>=min_same && r.n_diff>=min_diff;
        if r.enough
            r.mean_same = mean(vs,'omitnan');
            r.mean_diff = mean(vd,'omitnan');
            r.delta_obs = r.mean_same-r.mean_diff;
        else
            r.mean_same = NaN; r.mean_diff = NaN; r.delta_obs = NaN;
        end
        res.(scope).(met) = r;
    end
end

deltaNull = struct();
for s = 1:numel(scopeNames)
    for m = 1:numel(metrics)
        deltaNull.(scopeNames{s}).(metrics{m}) = nan(Nperm,1);
    end
end

for b = 1:Nperm
    perm = labels(randperm(N));
    sameLblB = perm==perm';
    for s = 1:numel(scopeNames)
        scope = scopeNames{s};
        sameMaskB = poolMask.(scope) & sameLblB;
        diffMaskB = poolMask.(scope) & ~sameLblB;
        for m = 1:numel(metrics)
            met = metrics{m};
            if ~res.(scope).(met).enough, continue; end
            fm = finiteMask.(met);
            vs = M.(met)(sameMaskB & fm);
            vd = M.(met)(diffMaskB & fm);
            if isempty(vs) || isempty(vd), continue; end
            deltaNull.(scope).(met)(b) = mean(vs,'omitnan')-mean(vd,'omitnan');
        end
    end
end

for s = 1:numel(scopeNames)
    scope = scopeNames{s};
    for m = 1:numel(metrics)
        met = metrics{m};
        dn = deltaNull.(scope).(met);
        ok = isfinite(dn);
        if res.(scope).(met).enough && any(ok)
            res.(scope).(met).p_value = (1+sum(dn(ok)>=res.(scope).(met).delta_obs))/(1+nnz(ok));
            res.(scope).(met).n_valid_perm = nnz(ok);
            if ismember(scope,keep_null_scopes)
                res.(scope).(met).delta_null = dn(ok);
            else
                res.(scope).(met).delta_null = [];
            end
        else
            res.(scope).(met).p_value = NaN;
            res.(scope).(met).n_valid_perm = 0;
            res.(scope).(met).delta_null = [];
        end
    end
end
end

function Xs = circshift_rows_independent(X,shifts)
[N,~] = size(X);
Xs = zeros(size(X),'like',X);
for i = 1:N
    Xs(i,:) = circshift(X(i,:),shifts(i),2);
end
end

function res = temporal_shuffle_all_metrics_scopes(FluoV,PeaksV,labels,poolMask, ...
    M0,metrics,Nperm,min_same,min_diff,width_frames,min_shift_frames,keep_null_scopes)
scopeNames = fieldnames(poolMask);
labels = labels(:);
[N,T] = size(FluoV);
sameLbl = labels==labels';

res = struct();
for s = 1:numel(scopeNames)
    scope = scopeNames{s};
    sameMask = poolMask.(scope) & sameLbl;
    diffMask = poolMask.(scope) & ~sameLbl;
    for m = 1:numel(metrics)
        met = metrics{m};
        fm = isfinite(M0.(met));
        vs = M0.(met)(sameMask & fm);
        vd = M0.(met)(diffMask & fm);
        r = struct('n_same',numel(vs),'n_diff',numel(vd));
        r.enough = r.n_same>=min_same && r.n_diff>=min_diff;
        if r.enough
            r.mean_same = mean(vs,'omitnan');
            r.mean_diff = mean(vd,'omitnan');
            r.delta_obs = r.mean_same-r.mean_diff;
        else
            r.mean_same = NaN; r.mean_diff = NaN; r.delta_obs = NaN;
        end
        res.(scope).(met) = r;
    end
end

deltaNull = struct();
for s = 1:numel(scopeNames)
    for m = 1:numel(metrics)
        deltaNull.(scopeNames{s}).(metrics{m}) = nan(Nperm,1);
    end
end

max_shift = max(min_shift_frames+1,T-min_shift_frames);
for b = 1:Nperm
    shifts = randi([min_shift_frames,max_shift],N,1);
    FluoB = circshift_rows_independent(FluoV,shifts);
    PeaksB = circshift_rows_independent(double(PeaksV),shifts) > 0;
    EventB = widen_binary_events(PeaksB,width_frames);
    Mb = compute_similarity_matrices(FluoB,EventB);

    for s = 1:numel(scopeNames)
        scope = scopeNames{s};
        sameMask = poolMask.(scope) & sameLbl;
        diffMask = poolMask.(scope) & ~sameLbl;
        for m = 1:numel(metrics)
            met = metrics{m};
            if ~res.(scope).(met).enough, continue; end
            fm = isfinite(Mb.(met));
            vs = Mb.(met)(sameMask & fm);
            vd = Mb.(met)(diffMask & fm);
            if isempty(vs) || isempty(vd), continue; end
            deltaNull.(scope).(met)(b) = mean(vs,'omitnan')-mean(vd,'omitnan');
        end
    end
end

for s = 1:numel(scopeNames)
    scope = scopeNames{s};
    for m = 1:numel(metrics)
        met = metrics{m};
        dn = deltaNull.(scope).(met);
        ok = isfinite(dn);
        if res.(scope).(met).enough && any(ok)
            res.(scope).(met).p_value = (1+sum(dn(ok)>=res.(scope).(met).delta_obs))/(1+nnz(ok));
            res.(scope).(met).n_valid_perm = nnz(ok);
            if ismember(scope,keep_null_scopes)
                res.(scope).(met).delta_null = dn(ok);
            else
                res.(scope).(met).delta_null = [];
            end
        else
            res.(scope).(met).p_value = NaN;
            res.(scope).(met).n_valid_perm = 0;
            res.(scope).(met).delta_null = [];
        end
    end
end
end

function plot_session_shuffle_histograms(colorRes,temporalRes,metrics,scope_label, ...
    identifier,pathout,show_fig)
if show_fig
    vis = 'on';
else
    vis = 'off';
end
f = figure('Position',[50 50 900 1000],'Color','w','Visible',vis);
nMet = numel(metrics);
for m = 1:nMet
    met = metrics{m};

    subplot(nMet,2,(m-1)*2+1); hold on
    cr = colorRes.(met);
    if cr.enough && ~isempty(cr.delta_null)
        histogram(cr.delta_null,40,'Normalization','probability', ...
            'FaceColor',[0.55 0.70 0.85],'EdgeColor','none');
    end
    xline(0,'--k'); xline(cr.delta_obs,'-r','LineWidth',2);
    title(sprintf('%s | shuffle couleurs | p=%.4g (n=%d/%d)', ...
        upper(met),cr.p_value,cr.n_same,cr.n_diff),'FontSize',9);
    if m==nMet, xlabel('Delta (même - différente)'); end
    ylabel('Probabilité');
    set(gca,'Box','off','TickDir','out','FontSize',8);

    subplot(nMet,2,(m-1)*2+2); hold on
    tr = temporalRes.(met);
    if tr.enough && ~isempty(tr.delta_null)
        histogram(tr.delta_null,40,'Normalization','probability', ...
            'FaceColor',[0.85 0.70 0.55],'EdgeColor','none');
    end
    xline(0,'--k'); xline(tr.delta_obs,'-r','LineWidth',2);
    title(sprintf('%s | shuffle temporel | p=%.4g (n=%d/%d)', ...
        upper(met),tr.p_value,tr.n_same,tr.n_diff),'FontSize',9);
    if m==nMet, xlabel('Delta (même - différente)'); end
    set(gca,'Box','off','TickDir','out','FontSize',8);
end
sgtitle(sprintf('%s — portée : %s',identifier,scope_label),'Interpreter','none');
exportgraphics(f,pathout,'Resolution',200);
close(f);
end

function [Tobs,Tnull,pval] = hierarchical_color_shuffle(results_all,metric,scope, ...
    animal_ids,Nperm,min_same,min_diff)
nSessions = numel(results_all);
session_animal = strings(nSessions,1);
Sim = cell(nSessions,1);
pool = cell(nSessions,1);
labels = cell(nSessions,1);
finiteMask = cell(nSessions,1);
usable = false(nSessions,1);

for k = 1:nSessions
    session_animal(k) = string(results_all(k).summary.Animal_ID);
    out = results_all(k).out;
    if ~isfield(out,'M') || ~isfield(out,'pool') || ~isfield(out.pool,scope)
        continue
    end
    Sim{k} = out.M.(metric);
    pool{k} = out.pool.(scope);
    labels{k} = out.colors(:);
    finiteMask{k} = isfinite(Sim{k});
    usable(k) = true;
end

deltaObsSession = nan(nSessions,1);
for k = 1:nSessions
    if ~usable(k), continue; end
    sameLbl = labels{k}==labels{k}';
    sameMask = pool{k} & sameLbl & finiteMask{k};
    diffMask = pool{k} & ~sameLbl & finiteMask{k};
    if nnz(sameMask)<min_same || nnz(diffMask)<min_diff, continue; end
    deltaObsSession(k) = mean(Sim{k}(sameMask),'omitnan')-mean(Sim{k}(diffMask),'omitnan');
end

Tobs = animal_mean_of_session_delta(deltaObsSession,session_animal,animal_ids);

Tnull = nan(Nperm,1);
for b = 1:Nperm
    deltaB = nan(nSessions,1);
    for k = 1:nSessions
        if ~usable(k) || isnan(deltaObsSession(k)), continue; end
        N = numel(labels{k});
        perm = labels{k}(randperm(N));
        sameLblB = perm==perm';
        sameMaskB = pool{k} & sameLblB & finiteMask{k};
        diffMaskB = pool{k} & ~sameLblB & finiteMask{k};
        if nnz(sameMaskB)==0 || nnz(diffMaskB)==0, continue; end
        deltaB(k) = mean(Sim{k}(sameMaskB),'omitnan')-mean(Sim{k}(diffMaskB),'omitnan');
    end
    Tnull(b) = animal_mean_of_session_delta(deltaB,session_animal,animal_ids);
end

ok = isfinite(Tnull);
if isnan(Tobs) || ~any(ok)
    pval = NaN;
else
    pval = (1+sum(Tnull(ok)>=Tobs))/(1+nnz(ok));
end
end

function T = animal_mean_of_session_delta(deltaSession,session_animal,animal_ids)
animalDeltas = nan(numel(animal_ids),1);
for a = 1:numel(animal_ids)
    idx = session_animal==animal_ids(a);
    vals = deltaSession(idx);
    vals = vals(isfinite(vals));
    if isempty(vals), continue; end
    animalDeltas(a) = mean(vals);
end
T = mean(animalDeltas,'omitnan');
end
