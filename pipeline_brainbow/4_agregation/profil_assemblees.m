function lignes = profil_assemblees(fichier_resultats, sessions, config, dossier_agreg)
%PROFIL_ASSEMBLEES  Étape 4 (22 nov. 2026, demandé par l'utilisateur) : deux
%   analyses exploratoires, écrites dans <dossier_agreg>\profil_assemblees\ :
%
%   1. PROFIL DES ASSEMBLÉES BIAISÉES : les assemblées testables biaisées
%      (p_minp < config.clonal.seuilstat, test principal) ressemblent-elles
%      aux autres ? Caractéristiques aveugles à la couleur (taille,
%      compacité, étalement, activations, SCE de la session, fraction de
%      membres actifs, % du champ, cellules du champ, nombre de K et
%      co-appartenance du consensus, fréquence des membres au repos, part de
%      participations en gros événements, plan) + nombre et part de cellules
%      colorées (puissance du test). Wilcoxon rang-somme biaisées / autres,
%      seuil de Bonferroni indiqué. Exploratoire : peu d'assemblées.
%      -> profil_assemblees_biaisees.csv / _detail.csv, biased_assembly_profile.png
%
%   2. TAILLE DES ÉVÉNEMENTS (détection CASCADE en PA discrets seulement) :
%      événement = suite d'images consécutives avec >= 1 PA ; au repos ;
%      petit (type ripple) = <= 3 PA et <= 3 images (~300 ms), gros sinon.
%      Distribution pour tous les événements et pour les participations aux
%      SCE (événement qui contient la fenêtre du SCE), par animal.
%      PA comme utilisés à l'étape 1 (images mauvaises et gros événements
%      éventuellement retirés par opts.cascade_max_* exclus).
%      -> taille_evenements_par_session.csv, event_size.png
%
%   fichier_resultats : resultats_biais_clonal.csv de la variante principale
%   sessions          : sessions finies (dossier, dossier_source, dossier_suite2p, id_*)
%   lignes            : texte pour le résumé.

lignes = {};
dossier = fullfile(dossier_agreg, 'profil_assemblees');
if ~isfolder(dossier), mkdir(dossier); end
seuil = 0.05;
if isfield(config, 'clonal') && isfield(config.clonal, 'seuilstat'), seuil = config.clonal.seuilstat; end

T = readtable(fichier_resultats, 'TextType', 'string');
T = T(T.testable == 1, :);
for ch = {'animal', 'session', 'champ', 'condition'}
    T.(ch{1}) = string(T.(ch{1}));
end
cle_T = T.animal + "|" + T.session + "|" + T.champ + "|" + T.condition;
cle_S = string(sessions.id_animal) + "|" + string(sessions.id_session) + "|" + string(sessions.id_champ) + ...
    "|" + string(sessions.id_condition);

T.nK_presence = nan(height(T), 1);
T.consensus_moyen = nan(height(T), 1);
T.freq_membres = nan(height(T), 1);
T.freq_membres_rel = nan(height(T), 1);
T.part_gros = nan(height(T), 1);
E = table();                                 % taille des événements, par session
tailles_tous = []; durees_tous = []; tailles_sce = []; durees_sce = []; animal_tous = strings(0); animal_sce = strings(0);

for s = 1:height(sessions)
    ds = sessions.dossier{s};
    if ismember('dossier_source', sessions.Properties.VariableNames), ds = sessions.dossier_source{s}; end
    fpre = fullfile(ds, 'pretraitement.mat');
    if ~isfile(fpre), continue; end
    P = load(fpre, 'Race', 'TRace', 'Raster', 'WinRest', 'sampling_rate', 'opts', 'bad_frames');
    min_repos = numel(P.WinRest) / P.sampling_rate / 60;
    freq = sum(P.Raster > 0, 2) / max(eps, min_repos);
    med_act = median(freq(freq > 0));
    % événements CASCADE (si cette détection)
    ev = evenements_cascade(P, char(sessions.dossier_suite2p{s}));
    if ~isempty(ev)
        repos = false(1, size(P.Raster, 2)); repos(P.WinRest) = true;
        r = repos(ev.debut)';
        tailles_tous = [tailles_tous; ev.taille(r)]; durees_tous = [durees_tous; ev.duree(r)]; %#ok<AGROW>
        animal_tous = [animal_tous; repmat(string(sessions.id_animal{s}), sum(r), 1)]; %#ok<AGROW>
        [pm, ps] = find(P.Race > 0);
        e_p = arrayfun(@(k) ev_de(ev, pm(k), P.TRace(ps(k)), P.opts.synchronous_frames), 1:numel(pm));
        e_p = e_p(e_p > 0);
        tailles_sce = [tailles_sce; ev.taille(e_p)]; durees_sce = [durees_sce; ev.duree(e_p)]; %#ok<AGROW>
        animal_sce = [animal_sce; repmat(string(sessions.id_animal{s}), numel(e_p), 1)]; %#ok<AGROW>
        t_r = ev.taille(r);
        petit_t = t_r <= 3 & ev.duree(r) <= 3;
        petit_p = ev.taille(e_p) <= 3 & ev.duree(e_p) <= 3;
        E = [E; table(string(sessions.nom{s}), string(sessions.id_animal{s}), sum(r), mean(petit_t), ...
            sum(t_r(petit_t)) / max(1, sum(t_r)), numel(e_p), mean(petit_p), ...
            'VariableNames', {'session', 'animal', 'n_evenements_repos', 'part_petits', 'part_PA_petits', ...
            'n_participations_SCE', 'part_participations_petits'})]; %#ok<AGROW>
    end
    % assemblées testables de cette session
    lignes_T = find(cle_T == cle_S(s));
    if isempty(lignes_T), continue; end
    for i = lignes_T'
        K = T.K(i);
        fa = fichier_assemblees(sessions.dossier{s}, ds, K);
        if ~isfile(fa), continue; end
        noms_var = intersect({'assemblyortho', 'PCl_final', 'K_presence', 'consensus_moyen'}, who('-file', fa));
        A = load(fa, noms_var{:});
        c = T.assemblee(i);
        if c > numel(A.assemblyortho), continue; end
        membres = A.assemblyortho{c};
        if isfield(A, 'K_presence') && numel(A.K_presence) >= c, T.nK_presence(i) = numel(A.K_presence{c}); end
        if isfield(A, 'consensus_moyen') && numel(A.consensus_moyen) >= c, T.consensus_moyen(i) = A.consensus_moyen(c); end
        T.freq_membres(i) = median(freq(membres));
        T.freq_membres_rel(i) = T.freq_membres(i) / med_act;
        if ~isempty(ev) && isfield(A, 'PCl_final') && size(A.PCl_final, 1) >= c
            sce = find(A.PCl_final(c, :) > 0);
            [mm, ss] = find(P.Race(membres, sce) > 0);
            e = arrayfun(@(k) ev_de(ev, membres(mm(k)), P.TRace(sce(ss(k))), P.opts.synchronous_frames), 1:numel(mm));
            e = e(e > 0);
            if ~isempty(e), T.part_gros(i) = mean(ev.taille(e) > 3 | ev.duree(e) > 3); end
        end
    end
end

%% ---------- 1. Profil des assemblées biaisées ---------------------------
T.plan1 = double(T.champ == "plane1");
carac = {'n_total', 'Size (cells)'; 'n_colorees', 'Coloured cells (power)'; 'part_colorees', 'Share coloured'; ...
    'compacite_pct', 'Compactness (%)'; 'ratio_etalement', 'Spread / chance'; 'n_SCE_signif', 'Activations'; ...
    'n_SCE_total', 'Session SCEs'; 'frac_active_moy', 'Active members / activation'; 'pct_reseau', '% of field cells'; ...
    'n_cellules_champ', 'Field cells'; 'nK_presence', 'K where seen (consensus)'; ...
    'consensus_moyen', 'Mean co-membership (consensus)'; 'freq_membres', 'Member rate (ev/min rest)'; ...
    'freq_membres_rel', 'Member rate / session'; 'part_gros', 'Share from large events'; 'plan1', 'Plane 1'};
carac = carac(ismember(carac(:, 1), T.Properties.VariableNames), :);
b = T.p_minp < seuil;
P = table();
for k = 1:size(carac, 1)
    v = double(T.(carac{k, 1}));
    ok = ~isnan(v);
    p = NaN;
    if sum(b & ok) >= 2 && sum(~b & ok) >= 2, p = ranksum(v(b & ok), v(~b & ok)); end
    P = [P; table(string(carac{k, 1}), string(carac{k, 2}), median(v(b & ok)), median(v(~b & ok)), p, ...
        'VariableNames', {'variable', 'label', 'median_biased', 'median_other', 'p_ranksum'})]; %#ok<AGROW>
end
writetable(P, fullfile(dossier, 'profil_assemblees_biaisees.csv'));
garder = intersect([{'animal', 'session', 'champ', 'assemblee', 'K', 'p_minp', 'q'}, carac(:, 1)'], ...
    T.Properties.VariableNames, 'stable');
writetable(T(:, garder), fullfile(dossier, 'profil_assemblees_detail.csv'));
nk = height(P);
if any(isfinite(P.p_ranksum))
    [pmin, imin] = min(P.p_ranksum);
    lignes{end+1} = sprintf(['Profil des assemblées biaisées (exploratoire, %d biaisées / %d testables, %d caractéristiques ' ...
        'sans couleur) : la plus différente est « %s » (médiane %.3g contre %.3g, Wilcoxon p = %.3f ; seuil de ' ...
        'Bonferroni %.4f).'], sum(b), height(T), nk, P.label(imin), P.median_biased(imin), P.median_other(imin), ...
        pmin, 0.05 / nk);
else
    lignes{end+1} = sprintf('Profil des assemblées biaisées : trop peu d''assemblées (%d biaisées / %d testables).', ...
        sum(b), height(T));
end
fig = figure('Color', 'w', 'Position', [50 50 1600 900], 'Visible', 'off');
nc = 6; nr = ceil(nk / nc);
tl = tiledlayout(fig, nr, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:nk
    ax = nexttile(tl); hold(ax, 'on');
    v = double(T.(P.variable(k)));
    rng(k);
    scatter(ax, 1 + 0.12 * randn(sum(~b), 1), v(~b), 18, [0.6 0.6 0.6], 'filled');
    scatter(ax, 2 + 0.12 * randn(sum(b), 1), v(b), 24, [0.85 0.2 0.2], 'filled');
    plot(ax, [0.7 1.3], P.median_other(k) * [1 1], 'k-', 'LineWidth', 1.5);
    plot(ax, [1.7 2.3], P.median_biased(k) * [1 1], 'k-', 'LineWidth', 1.5);
    ax.XTick = [1 2]; ax.XTickLabel = {'other', 'biased'}; xlim(ax, [0.4 2.6]); box(ax, 'off');
    title(ax, sprintf('%s\np = %.3f', P.label(k), P.p_ranksum(k)), 'FontWeight', 'normal', 'FontSize', 9);
end
title(tl, sprintf(['Testable assemblies: %d biased (p\\_minp < %g) vs %d other. Colour-blind features ' ...
    '(+ coloured cells = power). Rank-sum p; Bonferroni threshold %.4f'], sum(b), seuil, sum(~b), 0.05 / nk));
exportgraphics(fig, fullfile(dossier, 'biased_assembly_profile.png'), 'Resolution', 120);
close(fig);

%% ---------- 2. Taille des événements ------------------------------------
if isempty(tailles_tous)
    lignes{end+1} = 'Taille des événements : seulement avec la détection CASCADE en PA discrets.';
    return
end
writetable(E, fullfile(dossier, 'taille_evenements_par_session.csv'));
petit_t = tailles_tous <= 3 & durees_tous <= 3;
petit_s = tailles_sce <= 3 & durees_sce <= 3;
lignes{end+1} = sprintf(['Taille des événements au repos (CASCADE, événement = images consécutives avec PA) : ' ...
    '%.1f %% sont petits (<= 3 PA et <= 300 ms, type ripple) et portent %.1f %% des PA ; %.1f %% des ' ...
    'participations aux SCE viennent d''un petit événement.'], 100 * mean(petit_t), ...
    100 * sum(tailles_tous(petit_t)) / sum(tailles_tous), 100 * mean(petit_s));
fig = figure('Color', 'w', 'Position', [80 80 1500 450], 'Visible', 'off');
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact');
bords = 0.5:1:15.5;
for ip = 1:2
    ax = nexttile(tl);
    if ip == 1, x = tailles_tous; t = 'All events at rest'; else, x = tailles_sce; t = 'Events of SCE participations'; end
    h = histcounts(min(x, 15), bords, 'Normalization', 'probability');
    br = bar(ax, 1:15, 100 * h, 1);
    br.FaceColor = 'flat'; br.CData = repmat([0.6 0.6 0.6], 15, 1); br.CData(1:3, :) = repmat([0.2 0.45 0.8], 3, 1);
    xlabel(ax, 'Spikes in the event (15 = 15 or more)'); ylabel(ax, '% of events');
    title(ax, sprintf('%s: %.0f %% with 1-3 spikes', t, 100 * mean(x <= 3)), 'FontWeight', 'normal');
    box(ax, 'off');
end
ax = nexttile(tl);
an = unique(animal_tous);
v = [arrayfun(@(a) 100 * mean(~(tailles_tous(animal_tous == a) <= 3 & durees_tous(animal_tous == a) <= 3)), an), ...
     arrayfun(@(a) 100 * mean(~(tailles_sce(animal_sce == a) <= 3 & durees_sce(animal_sce == a) <= 3)), an)];
bar(ax, categorical(an), v);
legend(ax, {'events at rest', 'SCE participations'}, 'Location', 'northoutside', 'Box', 'off');
ylabel(ax, '% large (> 3 spikes or > 300 ms)'); box(ax, 'off');
title(ax, 'Large events by animal', 'FontWeight', 'normal');
title(tl, 'CASCADE events at rest. Event = consecutive frames with spikes; small (ripple-like) = <= 3 spikes and <= 300 ms');
exportgraphics(fig, fullfile(dossier, 'event_size.png'), 'Resolution', 130);
close(fig);
end


function ev = evenements_cascade(P, dossier_s2p)
% événements (images consécutives avec PA) d'une session détectée avec
% CASCADE en PA discrets ; [] sinon. PA comme à l'étape 1 (images mauvaises
% mises à 0 ; gros événements retirés si opts.cascade_max_* existe).
ev = [];
o = P.opts;
if ~isfield(o, 'detection') || ~strcmpi(o.detection, 'cascade_discrets'), return; end
modele = 'GC8_EXC_10Hz_smoothing150ms_high_noise';
if isfield(o, 'cascade_modele'), modele = o.cascade_modele; end
hz = [];
if isfield(o, 'cascade_reechantillonnage_Hz'), hz = o.cascade_reechantillonnage_Hz; end
suffixe = '';
if ~isempty(hz), suffixe = sprintf('_%gHz', hz); end
fc = fullfile(dossier_s2p, ['cascade_discrets_' modele suffixe '.mat']);
if ~isfile(fc), return; end
C = load(fc, 'n_pa');
n = double(C.n_pa);
[NC, Nz] = size(P.Raster);
if ~isempty(hz)
    t0 = (0:Nz - 1) / P.sampling_rate;
    t_r = 0:1 / hz:t0(end);
    io = min(Nz, round(t_r * P.sampling_rate) + 1);
    n = full(n * sparse(1:numel(io), io, 1, numel(io), Nz));
end
if ~isequal(size(n), [NC Nz]), return; end
if isfield(P, 'bad_frames') && numel(P.bad_frames) == Nz, n(:, logical(P.bad_frames)) = 0; end
ev = decouper(n);
max_pa = Inf; max_im = Inf;
if isfield(o, 'cascade_max_pa_evenement') && ~isempty(o.cascade_max_pa_evenement), max_pa = o.cascade_max_pa_evenement; end
if isfield(o, 'cascade_max_duree_evenement_s') && ~isempty(o.cascade_max_duree_evenement_s)
    max_im = round(o.cascade_max_duree_evenement_s * P.sampling_rate);
end
if isfinite(max_pa) || isfinite(max_im)
    gros = find(ev.taille > max_pa | ev.duree > max_im);
    for e = gros(:)', n(ev.ci(e), ev.debut(e):ev.fin(e)) = 0; end
    ev = decouper(n);
end
end

function ev = decouper(n)
[NC, Nz] = size(n);
d = diff([false(NC, 1), n > 0, false(NC, 1)], 1, 2);
[debut, ci] = find(d' == 1);          % ordre : cellule puis temps (débuts et fins)
[fin, ~] = find(d' == -1);
fin = fin - 1;
cs = cumsum([zeros(NC, 1), n], 2);
ev.ci = ci; ev.debut = debut; ev.fin = fin;
ev.taille = cs(sub2ind(size(cs), ci, fin + 1)) - cs(sub2ind(size(cs), ci, debut));
ev.duree = fin - debut + 1;
ev.id = zeros(NC, Nz);
for e = 1:numel(ci), ev.id(ci(e), debut(e):fin(e)) = e; end
end

function e = ev_de(ev, cellule, t_sce, sync)
% événement de la cellule qui contient la fenêtre du SCE (0 si aucun)
f = t_sce + (0:sync);
f = f(f >= 1 & f <= size(ev.id, 2));
x = ev.id(cellule, f); x = x(x > 0);
if isempty(x), e = 0; else, e = x(1); end
end
