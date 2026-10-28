clear; clc; close all

%% =========================================================
% CORRELATIONS / COACTIVATIONS PAR COULEUR
%
% Chargement : liste "matfile" contenue dans nocues_data_recording_all.mat.
% Pour chaque entrée, le script charge le fichier de données puis Fall.mat
% dans le dossier "path". Le masque iscell de Fall.mat est la référence.
%
% Analyse principale : phi sur événements binaires, sans filtre spatial.
% Analyses secondaires : mêmes calculs pour les paires <=100, <=200,
% et <=300 µm, si les nombres minimaux de paires sont disponibles.
% Robustesses : Spearman continu, Pearson continu et Jaccard binaire.
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
analysis.radii_um = [50 100 200 300];
analysis.min_same_pairs = 5;
analysis.min_diff_pairs = 5;
analysis.event_source = 'Raster';       % 'Raster' ou 'Race'
analysis.event_width_frames = 2;        % 1 = peak strict; 2 = tolérance d'une frame
analysis.min_events_per_cell = 1;
analysis.max_exact_animals = 20;
analysis.n_signflip_mc = 100000;
analysis.random_seed = 20260916;
analysis.stop_on_error = false;
rng(analysis.random_seed,'twister');

metrics = {'phi','spearman','pearson','jaccard'};
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

    % Ne ferme pas fig_log. Ferme seulement les figures MATLAB ordinaires.
    figs = findall(groot,'Type','figure');
    if ~isempty(figs), close(figs); end

    try
        historique = ecrire_log(txt_log,log_file,historique, ...
            sprintf('Traitement %d/%d : %s',file_num,numel(recording_files),filename));

        %% 4.1 Charge l'entrée data_recording
        [Session,recording_file] = load_recording(filename);
        identifier = string(get_required_field(Session,'identifier'));
        pathload = get_recording_path(Session,recording_file);

        if isfield(Session,'Animal_ID') && ~isempty(Session.Animal_ID) && string(Session.Animal_ID) ~= ""
            animal_id = string(Session.Animal_ID);
        else
            animal_id = identifier;
            historique = ecrire_log(txt_log,log_file,historique, ...
                sprintf('ATTENTION : Animal_ID absent pour %s; identifier utilisé.',identifier));
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
        % xpix/ypix sont préférés. +1 convertit les indices Python 0-based
        % de Suite2p vers les indices MATLAB 1-based.
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
            preprocessing_6(F,opts,speed,pathload,sampling_rate,namefull); %#ok<ASGLU>

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
        Fluo = double(Tr1b(:,WinRest));
        switch lower(analysis.event_source)
            case 'raster', EventRaw = Raster;
            case 'race', EventRaw = Race;
            otherwise, error('event_source doit être Raster ou Race.');
        end
        EventRaw = select_time_window(EventRaw,WinRest,size(Tr1b,2),'EventRaw');

        % Si Raster marque toute la durée d'un transient, chaque épisode
        % contigu est réduit à un seul peak au maximum local de Tr1b.
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

        %% 4.9 Similarités puis global / 100 / 200 / 300 µm
        M = compute_similarity_matrices(FluoV,EventV);
        out = struct();
        for m = 1:numel(metrics)
            met = metrics{m};
            out.(met) = summarize_same_diff_radii(M.(met),colorsV,DistMat, ...
                analysis.radii_um,analysis.min_same_pairs,analysis.min_diff_pairs);
        end
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

        %% 4.10 Résumé session
        row = struct('Animal_ID',animal_id,'identifier',identifier, ...
            'file_num',file_num,'n_cells',size(FluoV,1), ...
            'n_timepoints_rest',size(FluoV,2), ...
            'median_events_per_cell',median(eventsV,'omitnan'), ...
            'median_event_rate_hz',median(out.event_rate_hz,'omitnan'));

        for m = 1:numel(metrics)
            met = metrics{m};
            row.(['mean_same_' met '_global']) = out.(met).global.mean_same;
            row.(['mean_diff_' met '_global']) = out.(met).global.mean_diff;
            row.(['delta_' met '_global']) = out.(met).global.delta;
            row.(['n_same_' met '_global']) = out.(met).global.n_same;
            row.(['n_diff_' met '_global']) = out.(met).global.n_diff;

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
            end
        end

        summary_sessions = [summary_sessions;struct2table(row)]; %#ok<AGROW>
        kresult = numel(results_all)+1;
        results_all(kresult).identifier = identifier;
        results_all(kresult).filename = filename;
        results_all(kresult).out = out;
        results_all(kresult).summary = row;

        save(fullfile(PathSave,'results_all_backup.mat'), ...
            'results_all','summary_sessions','analysis','-v7.3');
        historique = ecrire_log(txt_log,log_file,historique, ...
            sprintf('OK %s terminé en %.2f min.',identifier,toc(temps_fichier)/60));

    catch ME
        historique = ecrire_log(txt_log,log_file,historique, ...
            sprintf('ERREUR fichier %d (%s) : %s',file_num,filename,ME.message));
        if analysis.stop_on_error, rethrow(ME); end
    end
end

if isempty(summary_sessions), error('Aucune session valide n''a été analysée.'); end

%% 5. AGRÉGATION NON PONDEREE PAR ANIMAL
animal_ids = unique(summary_sessions.Animal_ID,'stable');
summary_animals_true = table();
scopes = [{'global'},arrayfun(@(x)sprintf('R%d',x),analysis.radii_um,'UniformOutput',false)];

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
    summary_animals_true = [summary_animals_true;struct2table(row)]; %#ok<AGROW>
end

%% 6. TESTS AU NIVEAU ANIMAL
final_stats = table();
for m = 1:numel(metrics)
    met = metrics{m};
    for q = 1:numel(scopes)
        scope = scopes{q};
        D = summary_animals_true.(['delta_' met '_' scope]);
        st = signflip_test(D,'right',analysis.max_exact_animals,analysis.n_signflip_mc);
        try, p_wilcoxon = signrank(D,0,'tail','right');
        catch, p_wilcoxon = NaN;
        end
        newrow = table(string(met),string(scope),st.N,st.mean_delta, ...
            st.median_delta,st.p_value,p_wilcoxon,string(st.method), ...
            'VariableNames',{'Metric','Scope','N_animals','Mean_delta', ...
            'Median_delta','P_signflip','P_wilcoxon','Signflip_method'});
        final_stats = [final_stats;newrow]; %#ok<AGROW>
    end
end

final_stats.P_Holm_within_metric = nan(height(final_stats),1);
for m = 1:numel(metrics)
    idx = final_stats.Metric==string(metrics{m});
    final_stats.P_Holm_within_metric(idx) = holm_adjust(final_stats.P_signflip(idx));
end

disp(summary_animals_true)
disp(final_stats)

%% 7. FIGURES
for m = 1:numel(metrics)
    met = metrics{m};
    idx = final_stats.Metric==string(met) & final_stats.Scope=="global";
    plot_paired_animals(summary_animals_true.(['mean_same_' met '_global']), ...
        summary_animals_true.(['mean_diff_' met '_global']),metric_label(met), ...
        sprintf('%s global : p_{sign-flip}=%.4g, N=%d',upper(met), ...
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
save(fullfile(PathSave,'results_all_final_events_radii.mat'), ...
    'results_all','summary_sessions','summary_animals_true','final_stats','analysis','-v7.3');
writetable(summary_sessions,fullfile(PathSave,'summary_sessions_events_radii.csv'));
writetable(summary_animals_true,fullfile(PathSave,'summary_animals_events_radii.csv'));
writetable(final_stats,fullfile(PathSave,'final_stats_events_radii.csv'));
historique = ecrire_log(txt_log,log_file,historique, ...
    sprintf('=== FIN : durée totale %.2f min ===',toc(temps_global)/60));

%% =========================================================
% FONCTIONS LOCALES
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

function X = filter_suite2p_rows(X,mask,n_cells,name)
    % Accepte un vecteur ligne ou colonne pour colorcell et une matrice
    % [ROI x temps] pour F. Ne réapplique pas iscell si déjà filtré.
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
        T=mean(S.*D',2); st.method='exact'; mc=false;
    else
        S=2*(rand(nmc,numel(D))>=.5)-1;
        T=mean(S.*D',2); st.method='Monte-Carlo'; mc=true;
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
    f=figure('Position',[100 100 520 500],'Color','w'); hold on
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
end

function plot_radius_deltas(D,radii,ylab,tit,pathout)
    f=figure('Position',[100 100 620 500],'Color','w'); hold on
    for i=1:size(D,1), plot(radii,D(i,:),'-o','Color',[.75 .75 .75]); end
    n=sum(isfinite(D),1); mu=mean(D,1,'omitnan');
    sem=std(D,0,1,'omitnan')./sqrt(n);
    errorbar(radii,mu,sem,'k-o','LineWidth',2,'MarkerFaceColor','k','CapSize',8);
    yline(0,'k--'); xticks(radii); xlabel('Distance maximale (µm)');
    ylabel(['Delta ' ylab ' (même - différente)']); title(tit);
    set(gca,'Box','off','TickDir','out','FontName','Arial');
    exportgraphics(f,pathout,'Resolution',300);
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
