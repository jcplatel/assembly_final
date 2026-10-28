clear; clc; close all

%% ANALYSE COACTIVITE ET CORRELATIONS PAR COULEUR
% Analyse principale : phi sur événements binaires, sans filtre de distance.
% Analyses spatiales secondaires : mêmes calculs pour les paires <=100, <=200,
% et <=300 µm, uniquement si les nombres minimaux de paires sont atteints.
% Robustesses : Spearman continu, Pearson continu et Jaccard binaire.

%% 1. PARAMETRES
path_database = 'E:\Data\Aurelie\data\Global_Database.mat';
Pathsave = 'E:\Data\Aurelie\analysis\Final\correlation_events\';
if ~exist(Pathsave,'dir'), mkdir(Pathsave); end

analysis = struct();
analysis.pixel_size_um = 1.37;
analysis.radii_um = [100 200 300];
analysis.min_same_pairs = 5;
analysis.min_diff_pairs = 5;
analysis.event_source = 'Raster';       % 'Raster' ou 'Race'
analysis.event_width_frames = 2;        % fenêtre causale; 1 = pic strict
analysis.min_events_per_cell = 1;
analysis.max_exact_animals = 20;
analysis.n_signflip_mc = 100000;
analysis.random_seed = 20260916;
rng(analysis.random_seed,'twister');

metrics = {'phi','spearman','pearson','jaccard'};
sessions_a_traiter = recuperer_sessions_bdd(path_database,'nocues','Analysis_Correlation');
if isempty(sessions_a_traiter), error('Aucune session à traiter.'); end

results_all = struct([]);
summary_sessions = table();

%% 2. BOUCLE SUR LES SESSIONS
for file_num = 1:numel(sessions_a_traiter)
    t0 = tic;
    varname = sessions_a_traiter{file_num};
    fprintf('\nTraitement de %s (%d/%d)...\n',varname,file_num,numel(sessions_a_traiter));
    tmp = load(path_database,varname);
    Session = tmp.(varname);
    identifier = string(Session.identifier);

    if isfield(Session,'Animal_ID') && ~isempty(Session.Animal_ID) && string(Session.Animal_ID) ~= ""
        animal_id = string(Session.Animal_ID);
    else
        animal_id = identifier;
        warning('Animal_ID absent pour %s : identifier utilisé.',identifier);
    end

    F_all = Session.F;
    cell_mask = Session.iscell(:,1) > 0;
    F = F_all(cell_mask,:);
    speed = Session.speed;
    sampling_rate = Session.sampling_rate;

    if isfield(Session,'path') && ~isempty(Session.path)
        pathload = strcat(Session.path,'\');
    else
        pathload = strcat(Session.DossierSource,'\');
    end

    namefull = fullfile(Pathsave,char(identifier + "_" + datestr(now,'yy_mm_dd_HH_MM_SS')));
    if ~exist(namefull,'dir'), mkdir(namefull); end

    opts = struct('MinPeakDistancesce',5,'MinPeakDistance',3, ...
        'threshold_peak',3.09,'synchronous_frames',2, ...
        'sce_n_cells_threshold',10,'SG_window',7, ...
        'motion_correction',true,'colorsubstraction',false);

    [Tr1b,speedsm,Raster,SumAct,MAct,Race,RasterRace,WinRest,WinActive, ...
        TRace,Fzero,DFF0,th_detection,bad_frames,max_cells_allowed, ...
        n_transients_total,opts] = ...
        preprocessing_6(F,opts,speed,pathload,sampling_rate,namefull); %#ok<ASGLU>

    %% Couleurs
    if isfield(Session,'Colorcell_Version') && ~strcmp(Session.Colorcell_Version,'Non trouvé')
        colorfile = fullfile(pathload,Session.Colorcell_Version);
    else
        colorfile = fullfile(pathload,'colorcell_Maxwell.mat');
    end
    C = load(colorfile,'colorcell');
    colorcell = align_rows(C.colorcell,cell_mask,size(Tr1b,1),'colorcell');
    colorcell = double(colorcell(:));

    %% Coordonnées : inutiles pour le global, requises seulement pour les rayons
    coords = get_centroids(Session,cell_mask,size(Tr1b,1));
    if isempty(coords)
        warning(['Coordonnées absentes pour %s : analyse globale effectuée, ', ...
            'analyses 100/200/300 µm laissées à NaN.'],identifier);
    else
        coords = double(coords(:,1:2))*analysis.pixel_size_um;
    end

    %% Traces et événements pendant le repos
    Fluo = double(Tr1b(:,WinRest));
    switch lower(analysis.event_source)
        case 'raster', EventRaw = Raster;
        case 'race',   EventRaw = Race;
        otherwise, error('event_source doit être Raster ou Race.');
    end
    EventRaw = select_window(EventRaw,WinRest,size(Tr1b,2),'EventRaw');

    % Un seul événement au maximum de Tr1b pour chaque épisode détecté.
    EventPeaks = runs_to_peaks(EventRaw>0,Fluo);
    EventRaster = widen_events(EventPeaks,analysis.event_width_frames);
    n_events = sum(EventPeaks,2);

    %% Cellules valides
    valid = isfinite(colorcell) & colorcell>0 & ...
        n_events>=analysis.min_events_per_cell & std(EventRaster,0,2)>0;
    if nnz(valid)<3
        warning('%s ignorée : moins de 3 cellules valides.',identifier);
        continue
    end

    FluoV = Fluo(valid,:);
    EventV = EventRaster(valid,:);
    PeaksV = EventPeaks(valid,:);
    colorsV = colorcell(valid);
    eventsV = n_events(valid);

    if isempty(coords)
        DistMat = [];
    else
        coordsV = coords(valid,:);
        DistMat = squareform(pdist(coordsV));
    end

    %% Similarités puis résumés global / rayons
    M = similarity_matrices(FluoV,EventV);
    out = struct();
    for m = 1:numel(metrics)
        met = metrics{m};
        out.(met) = summarize_radii(M.(met),colorsV,DistMat, ...
            analysis.radii_um,analysis.min_same_pairs,analysis.min_diff_pairs);
    end
    out.event_peaks = PeaksV;
    out.event_raster = EventV;
    out.event_count = eventsV;
    out.event_rate_hz = eventsV/(size(FluoV,2)/sampling_rate);
    out.colors = colorsV;
    out.distance_matrix_um = DistMat;
    out.options = analysis;

    %% Ligne session
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
    results_all(end+1).identifier = identifier; %#ok<SAGROW>
    results_all(end).filename = varname;
    results_all(end).out = out;
    results_all(end).summary = row;
    save(fullfile(Pathsave,'results_all_backup.mat'), ...
        'results_all','summary_sessions','analysis','-v7.3');
    fprintf('=== %s terminé en %.2f min ===\n',identifier,toc(t0)/60);
end

if isempty(summary_sessions), error('Aucune session valide.'); end

%% 3. AGRÉGATION NON PONDEREE DES SESSIONS PAR ANIMAL
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

disp(summary_animals_true)

%% 4. TESTS AU NIVEAU ANIMAL
% Test directionnel préspécifié : même couleur > couleur différente.
% Global = analyse principale. R100/R200/R300 = secondaires et emboîtées.
final_stats = table();
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
        newrow = table(string(met),string(scope),st.N,st.mean_delta, ...
            st.median_delta,st.p_value,p_wilcoxon,string(st.method), ...
            'VariableNames',{'Metric','Scope','N_animals','Mean_delta', ...
            'Median_delta','P_signflip','P_wilcoxon','Signflip_method'});
        final_stats = [final_stats;newrow]; %#ok<AGROW>
    end
end

% Correction de Holm entre global/R100/R200/R300 pour chaque métrique.
% Si phi-global est préspécifié comme seul endpoint principal, garder sa p brute.
final_stats.P_Holm_within_metric = nan(height(final_stats),1);
for m = 1:numel(metrics)
    idx = final_stats.Metric==string(metrics{m});
    final_stats.P_Holm_within_metric(idx) = holm_adjust(final_stats.P_signflip(idx));
end
disp(final_stats)

%% 5. FIGURES GLOBALES
for m = 1:numel(metrics)
    met = metrics{m};
    idx = final_stats.Metric==string(met) & final_stats.Scope=="global";
    plot_paired(summary_animals_true.(['mean_same_' met '_global']), ...
        summary_animals_true.(['mean_diff_' met '_global']),metric_label(met), ...
        sprintf('%s global : p_{sign-flip}=%.4g, N=%d',upper(met), ...
        final_stats.P_signflip(idx),final_stats.N_animals(idx)), ...
        fullfile(Pathsave,sprintf('Figure_%s_global.png',met)));
end

%% 6. FIGURES DES RAYONS CUMULATIFS
for m = 1:numel(metrics)
    met = metrics{m};
    D = nan(height(summary_animals_true),numel(analysis.radii_um));
    for r = 1:numel(analysis.radii_um)
        D(:,r) = summary_animals_true.(['delta_' met sprintf('_R%d',analysis.radii_um(r))]);
    end
    plot_radii(D,analysis.radii_um,metric_label(met), ...
        sprintf('%s : paires dans le rayon indiqué',upper(met)), ...
        fullfile(Pathsave,sprintf('Figure_%s_radii.png',met)));
end

%% 7. SAUVEGARDE
save(fullfile(Pathsave,'results_all_final_events_radii.mat'), ...
    'results_all','summary_sessions','summary_animals_true','final_stats','analysis','-v7.3');
writetable(summary_sessions,fullfile(Pathsave,'summary_sessions_events_radii.csv'));
writetable(summary_animals_true,fullfile(Pathsave,'summary_animals_events_radii.csv'));
writetable(final_stats,fullfile(Pathsave,'final_stats_events_radii.csv'));
fprintf('\nTerminé. Primaire : phi global. Secondaires : <=100, <=200, <=300 µm.\n');

%% FONCTIONS LOCALES
function x = align_rows(x,mask,n,name)
    x = x(:);
    if numel(x)==numel(mask), x=x(mask);
    elseif numel(x)~=n, error('%s : taille incompatible.',name); end
end

function coords = get_centroids(Session,mask,n)
    names = {'gcamp_centroid','centroids','centroid','Coords','coords'};
    coords = [];
    for k=1:numel(names)
        if isfield(Session,names{k}) && ~isempty(Session.(names{k}))
            coords=Session.(names{k}); break
        end
    end
    if isempty(coords), return; end
    if size(coords,1)==numel(mask), coords=coords(mask,:);
    elseif size(coords,1)~=n, error('Taille des coordonnées incompatible.'); end
    if size(coords,2)<2, error('Deux coordonnées au minimum sont nécessaires.'); end
end

function X = select_window(X,win,n_total,name)
    if size(X,2)==n_total, X=X(:,win);
    elseif size(X,2)~=numel(win), error('%s : nombre de frames incompatible.',name); end
end

function Peaks = runs_to_peaks(B,Signal)
    B=logical(B); Peaks=false(size(B));
    if ~isequal(size(B),size(Signal)), error('Raster et Fluo doivent avoir la même taille.'); end
    for i=1:size(B,1)
        starts=find(B(i,:) & [true ~B(i,1:end-1)]);
        stops=find(B(i,:) & [~B(i,2:end) true]);
        for e=1:numel(starts)
            fr=starts(e):stops(e); [~,loc]=max(Signal(i,fr)); Peaks(i,fr(loc))=true;
        end
    end
end

function B = widen_events(Peaks,width)
    width=max(1,round(width));
    B=filter(ones(1,width),1,double(Peaks),[],2)>0;
end

function M = similarity_matrices(Fluo,Events)
    M.pearson=corr(Fluo','Type','Pearson','Rows','pairwise');
    M.spearman=corr(Fluo','Type','Spearman','Rows','pairwise');
    M.phi=corr(double(Events'),'Type','Pearson','Rows','pairwise');
    B=double(Events); inter=B*B'; act=sum(B,2); uni=act+act'-inter;
    M.jaccard=inter./uni; M.jaccard(uni==0)=NaN;
    M.jaccard(1:size(B,1)+1:end)=1;
end

function out = summarize_radii(Sim,labels,Dist,radii,min_same,min_diff)
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
    if isempty(Dist), return; end

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
    a=zeros(m,1); run=0;
    for i=1:m, run=max(run,(m-i+1)*ps(i)); a(i)=min(1,run); end
    tmp=nan(m,1); tmp(ord)=a; padj(ok)=tmp;
end

function txt = metric_label(m)
    switch lower(m)
        case 'phi', txt='Coefficient phi (événements binaires)';
        case 'spearman', txt='Corrélation de Spearman (traces continues)';
        case 'pearson', txt='Corrélation de Pearson (traces continues)';
        case 'jaccard', txt='Indice de Jaccard (événements binaires)';
    end
end

function plot_paired(a,b,ylab,tit,pathout)
    ok=isfinite(a)&isfinite(b); a=a(ok); b=b(ok);
    f=figure('Position',[100 100 520 500],'Color','w'); hold on
    for i=1:numel(a)
        plot([1 2],[a(i) b(i)],'-','Color',[.72 .72 .72]);
        plot(1,a(i),'o','MarkerFaceColor',[.15 .45 .75],'Color',[.15 .45 .75]);
        plot(2,b(i),'o','MarkerFaceColor',[.85 .35 .20],'Color',[.85 .35 .20]);
    end
    mu=[mean(a,'omitnan') mean(b,'omitnan')];
    se=[std(a,0,'omitnan') std(b,0,'omitnan')]/sqrt(numel(a));
    errorbar([1 2],mu,se,'k-','LineWidth',2,'CapSize',8);
    xlim([.5 2.5]); xticks([1 2]); xticklabels({'Même couleur','Différente'});
    ylabel(ylab); title(tit); set(gca,'Box','off','TickDir','out','FontName','Arial');
    exportgraphics(f,pathout,'Resolution',300);
end

function plot_radii(D,radii,ylab,tit,pathout)
    f=figure('Position',[100 100 620 500],'Color','w'); hold on
    for i=1:size(D,1), plot(radii,D(i,:),'-o','Color',[.75 .75 .75]); end
    n=sum(isfinite(D),1); mu=mean(D,1,'omitnan'); se=std(D,0,1,'omitnan')./sqrt(n);
    errorbar(radii,mu,se,'k-o','LineWidth',2,'MarkerFaceColor','k','CapSize',8);
    yline(0,'k--'); xticks(radii); xlabel('Distance maximale (µm)');
    ylabel(['Delta ' ylab ' (même - différente)']); title(tit);
    set(gca,'Box','off','TickDir','out','FontName','Arial');
    exportgraphics(f,pathout,'Resolution',300);
end
