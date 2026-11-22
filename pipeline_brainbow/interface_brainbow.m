function interface_brainbow()
%INTERFACE_BRAINBOW  Graphical interface of the brainbow pipeline.
%
%   interface_brainbow
%
%   Left: settings for the run, buttons, session table (state of each
%   step), log. Right, for the selected session:
%     - speed, with the time of each SCE (triangle in the colour of the
%       active assembly, grey otherwise)
%     - activity raster, full width: transient peaks over the whole session
%       (running included), cells sorted by assembly; coloured: events of an
%       assembly's cells during the SCEs where that assembly is active;
%       zoom and pan linked with the speed
%     - below: the SCE raster sorted by assembly, the GCaMP image with the
%       cell masks coloured by assembly, and the brainbow colours of the
%       cells (or a colour image you choose) with the assembly outlines
%     - at the bottom: session panel (displayed K, assembly, summary, score of each K)
%   During an analysis, the display follows the session that just finished.
%
%   Folders: "Step 1 folder" holds the preprocessing and the assemblies of
%   each K (one subfolder per session) and the main analysis (steps 2-4).
%   "Analysis" chooses where steps 2-4 are written and read: main, or
%   <step 1 folder>\analyses\<name> ("New..."), on the same step 1.
%   "Analyses table..." (22 nov.): every analysis of ...\nocues and ...\cues
%   with a code (A001...), its settings and results (registre_analyses);
%   choose one there. The analysis shown is restored at the next opening.
%
%   The menu next to "Overall setting" chooses the condition shown: no cues
%   (sans_indice) or cues (avec_indice); the step 1 folder follows
%   (...\nocues\<folder> or ...\cues\<folder>).
%
%   Colours are shown only once step 3 is done for the session: the choice
%   of the settings and of K must stay blind to the colours.
%   Interface settings apply to the run; they do not modify
%   config_brainbow.m (a run with other settings redoes step 1).
%   Default values come from config_brainbow.m, read each time the
%   interface opens.

dossier_pipeline = fileparts(mfilename('fullpath'));
addpath(genpath(dossier_pipeline));

config0 = config_brainbow();
% 19 nov. (demandé par l'utilisateur) : condition affichée (sessions de la
% table), au choix dans l'interface ; config_brainbow.m n'est pas modifié
condition_courante = 'sans_indice';
if numel(config0.conditions) == 1, condition_courante = config0.conditions{1}; end
dossier_resultats = dossier_condition();   % step 1 folder (one subfolder per session); main analysis there too
nom_analyse = '';                       % steps 2-4: '' = main analysis, otherwise <step 1 folder>\analyses\<name>
noms_analyses = {};                     % analyses found in <step 1 folder>\analyses
minimum_SCE_calcule = 5;      % absolute minimum of cells in an SCE when the threshold is computed
                              % (used only if config_brainbow uses 'fixe')
taille_point_gris    = 1.5;   % activity raster: size of the grey points
taille_point_couleur = 2.5;   % activity raster: size of the coloured points
sessions = table();
selection = [];
en_cours = false;
arret_demande = false;
cache_geo = containers.Map('KeyType', 'char', 'ValueType', 'any');
images_couleur = containers.Map('KeyType', 'char', 'ValueType', 'any');
cache_act = containers.Map('KeyType', 'char', 'ValueType', 'any');
vue_coul = struct('nA', 0, 'cc', []);   % 3rd image: one assembly in brainbow colours (8 nov.)
idx_coul = 1;                           % assembly shown in the 3rd image
cle_vue_coul = '';                      % session | K of that view
cache_p = containers.Map('KeyType', 'char', 'ValueType', 'any');   % p and p_minp computed for a session | K

%% ======================= WINDOW =======================================
% Window size: the screen, at most 1860 x 1000
ecran = get(groot, 'ScreenSize');
taille_fenetre = [min(1860, ecran(3) - 60), min(1000, ecran(4) - 110)];
fig = uifigure('Name', 'Brainbow: assemblies and clonal analysis', 'Position', [30 50 taille_fenetre]);
g = uigridlayout(fig, [1 2]);
g.ColumnWidth = {440, '1x'};

%% ---------- Left column ----------
gauche = uigridlayout(g, [5 1]);
gauche.RowHeight = {351, 196, 62, '1x', 70};    % settings, actions, progress, sessions, log
gauche.Padding = [0 0 0 0];

% Settings
p_reg = uipanel(gauche, 'Title', 'Settings for this run');
gr = uigridlayout(p_reg, [13 4]);
gr.RowHeight = repmat({21}, 1, 13);
gr.ColumnWidth = {'1x', 95, 90, 60};
gr.RowSpacing = 3;
placer(uilabel(gr, 'Text', 'Overall setting'), 1, 1);
dd_reg = placer(uidropdown(gr, 'Items', {'validated', 'old + rank'}, ...
    'ItemsData', {'valide', 'ancien_rang'}, 'Value', config0.reglage, ...
    'ValueChangedFcn', @(~, ~) changer_reglage(), ...
    'Tooltip', sprintf(['validated: curveball, rank, weighted RACEORTHO, global FDR (no false ' ...
    'assembly on shuffled 444119)\nold + rank: old null models + rank (4 assemblies on ' ...
    '444119 at K = 9, but 1.3 false assemblies per shuffle)\nResults in separate folders'])), 1, 2);
dd_cond = placer(uidropdown(gr, 'Items', {'no cues', 'cues'}, ...
    'ItemsData', {'sans_indice', 'avec_indice'}, 'Value', condition_courante, ...
    'ValueChangedFcn', @(src, ~) changer_condition(src.Value), ...
    'Tooltip', sprintf(['Sessions shown and analysed: no cues (sans_indice) or cues (avec_indice).\n' ...
    'Switches the step 1 folder to the default folder of that condition\n' ...
    '(...\\nocues\\%s or ...\\cues\\%s). config_brainbow.m is not changed.'], ...
    nom_dossier(config0.PathSave), nom_dossier(config0.PathSave))), 1, [3 4]);
placer(uilabel(gr, 'Text', 'SCE threshold'), 2, 1);
regle_initiale = 'percentile';
if isfield(config0.pre, 'sce_regle') && ismember(config0.pre.sce_regle, {'percentile', 'fixe', 'fdr'})
    regle_initiale = config0.pre.sce_regle;
end
% Value of each rule, kept when switching from one to another:
%   percentile: percentile of the shuffles; fixe: number of cells; fdr: share of chance SCEs
seuils_sce = struct('percentile', 99, 'fixe', 15, 'fdr', 0.05);
if isfield(config0.pre, 'sce_percentile'), seuils_sce.percentile = config0.pre.sce_percentile; end
if isfield(config0.pre, 'sce_fdr'), seuils_sce.fdr = config0.pre.sce_fdr; end
if strcmp(regle_initiale, 'fixe'), seuils_sce.fixe = config0.pre.sce_n_cells_threshold; end
regle_affichee = regle_initiale;
dd_sce = placer(uidropdown(gr, 'Items', {'percentile', 'fixed', 'FDR'}, ...
    'ItemsData', {'percentile', 'fixe', 'fdr'}, 'Value', regle_initiale, ...
    'ValueChangedFcn', @(~, ~) changer_regle_sce(), ...
    'Tooltip', sprintf(['percentile: threshold = percentile of the shuffles (jitter +/- 10 frames)\n' ...
    'fixed: threshold = given number of cells\n' ...
    'FDR: smallest threshold where chance SCEs stay below this share'])), 2, 2);
lbl_sce = placer(uilabel(gr, 'Text', '', 'HorizontalAlignment', 'right'), 2, 3);
ef_sce = placer(uieditfield(gr, 'numeric', 'Value', 1), 2, 4);
placer(uilabel(gr, 'Text', 'Cluster threshold'), 3, 1);
seuil_cl0 = 'rang';
if isfield(config0.asm, 'seuil_clusters'), seuil_cl0 = config0.asm.seuil_clusters; end
dd_cl = placer(uidropdown(gr, 'Items', {'rank', 'max', 'all'}, 'ItemsData', {'rang', 'max', 'tous'}, ...
    'Value', seuil_cl0, ...
    'Tooltip', sprintf(['Significance threshold of the SCE clusters (step 1)\n' ...
    'rank: i-th cluster against the i-th surrogate cluster\n' ...
    'max: all against the best surrogate cluster (strict)\nall: old rule (all clusters pooled)'])), 3, 2);
pct0 = 95;
if isfield(config0.asm, 'percentile_seuil'), pct0 = config0.asm.percentile_seuil; end
placer(uilabel(gr, 'Text', 'percentile:', 'HorizontalAlignment', 'right'), 3, 3);
ef_pct = placer(uieditfield(gr, 'numeric', 'Value', pct0, 'Limits', [50 99.9], ...
    'Tooltip', ['Percentile of the surrogates used as cluster threshold (95 recommended; ' ...
    '99 needs at least 500 surrogates)']), 3, 4);
placer(uilabel(gr, 'Text', 'K computed, from'), 4, 1);
ef_kmin = placer(uieditfield(gr, 'numeric', 'Value', min(config0.liste_K), 'Limits', [2 60], ...
    'RoundFractionalValues', 'on'), 4, 2);
placer(uilabel(gr, 'Text', 'to', 'HorizontalAlignment', 'right'), 4, 3);
ef_kmax = placer(uieditfield(gr, 'numeric', 'Value', max(config0.liste_K), 'Limits', [2 60], ...
    'RoundFractionalValues', 'on'), 4, 4);
regleK0 = 'exces';
if isfield(config0.choixK, 'regle') && ismember(config0.choixK.regle, {'exces', 'consensus', 'score'})
    regleK0 = config0.choixK.regle;
end
placer(uilabel(gr, 'Text', 'Step 2: K choice', 'FontWeight', 'bold'), 5, 1);
dd_regleK = placer(uidropdown(gr, 'Items', {'best K', 'consensus of the K (no single K)', 'old score'}, ...
    'ItemsData', {'exces', 'consensus', 'score'}, 'Value', regleK0, ...
    'ValueChangedFcn', @(~, ~) changer_regleK(), ...
    'Tooltip', sprintf(['What step 2 does, with the K above (from "K computed, from" to "to"):\n' ...
    'best K: keeps the K whose significant clusters exceed chance the most\n' ...
    'consensus: assemblies whose cells stay together across all these K (no single K)\n' ...
    'old score: weighted score (silhouette, recall, specificity, purity)\n' ...
    'Choosing an analysis restores the choice it was made with.'])), 5, [2 4]);
rep_surr0 = config0.asm.kmean_iter;
if isfield(config0.asm, 'kmean_iter_surrogate'), rep_surr0 = config0.asm.kmean_iter_surrogate; end
placer(uilabel(gr, 'Text', 'k-means replicates'), 6, 1);
ef_rep = placer(uieditfield(gr, 'numeric', 'Value', config0.asm.kmean_iter, 'Limits', [1 2000], ...
    'RoundFractionalValues', 'on', 'Tooltip', 'Replicates of the k-means on the data'), 6, 2);
placer(uilabel(gr, 'Text', 'surrogates:', 'HorizontalAlignment', 'right'), 6, 3);
ef_nsur = placer(uieditfield(gr, 'numeric', 'Value', config0.asm.kmeans_surrogate, 'Limits', [1 5000], ...
    'RoundFractionalValues', 'on', 'Tooltip', 'Number of surrogates per K (at least 100 recommended)'), 6, 4);
placer(uilabel(gr, 'Text', 'Surrogate replicates'), 7, 1);
ef_rep_sur = placer(uieditfield(gr, 'numeric', 'Value', rep_surr0, 'Limits', [1 2000], ...
    'RoundFractionalValues', 'on', 'Tooltip', 'Replicates of the k-means on each surrogate'), 7, 2);
placer(uilabel(gr, 'Text', 'Cell FDR:', 'HorizontalAlignment', 'right'), 7, 3);
ef_fdr = placer(uieditfield(gr, 'numeric', 'Value', config0.asm.FDR_Q, 'Limits', [0.001 0.5], ...
    'Tooltip', ['RACEORTHO FDR (step 1): which cells belong to an assembly ' ...
    '(and which SCEs are activations, unless config.asm.FDR_Q_activations exists)']), 7, 4);
cb_e1 = placer(uicheckbox(gr, 'Text', 'Step 1: assemblies', 'Value', config0.faire_etape1), 8, [1 2]);
cb_e2 = placer(uicheckbox(gr, 'Text', 'Step 2: K choice', 'Value', config0.faire_etape2), 8, [3 4]);
cb_e3 = placer(uicheckbox(gr, 'Text', 'Step 3: clonal', 'Value', config0.faire_etape3), 9, [1 2]);
cb_e3b = placer(uicheckbox(gr, 'Text', 'Step 3b: pairs', 'Value', valeur_3b(config0), ...
    'Tooltip', sprintf(['Same-colour vs different-colour neighbours: synchrony at rest (STTC) by distance,\n' ...
    'without the assemblies. Per session: <session>\\paires ; all sessions: paires_agregation'])), 9, [3 4]);
cb_e4 = placer(uicheckbox(gr, 'Text', 'Step 4: aggregation', 'Value', config0.faire_etape4), 10, [1 2]);
cb_e5 = placer(uicheckbox(gr, 'Text', 'Step 5: sequences', 'Value', valeur_5(config0), ...
    'Tooltip', sprintf(['Cells activated in sequence during running (method of sequence_2026_01_27),\n' ...
    'order learnt on half of the runs and checked on the other half;\n' ...
    'then distance in the sequence, same colour vs different colour.\n' ...
    'Per session: <session>\\sequences ; all sessions: sequences_agregation'])), 10, [3 4]);
cb_fig = placer(uicheckbox(gr, 'Text', 'PNG figures', 'Value', config0.faire_figures), 11, 1);
cb_forcer = placer(uicheckbox(gr, 'Text', 'Recompute all (ignore cache)', ...
    'Value', config0.forcer_pretraitement || config0.forcer_assemblees), 11, [2 4]);
placer(uilabel(gr, 'Text', 'Step 1 folder', 'Tooltip', ...
    sprintf(['Preprocessing and assemblies of each K (one subfolder per session).\n' ...
    'The main analysis (steps 2-4) is in the same folder.'])), 12, 1);
ef_dossier = placer(uieditfield(gr, 'text', 'Value', dossier_resultats, 'Tooltip', dossier_resultats, ...
    'ValueChangedFcn', @(src, ~) changer_dossier(src.Value)), 12, [2 3]);
g_dos = placer(uigridlayout(gr, [1 2], 'Padding', 0, 'ColumnSpacing', 3), 12, 4);
placer(uibutton(g_dos, 'Text', '...', 'Tooltip', 'Choose the step 1 folder (one subfolder per session)', ...
    'ButtonPushedFcn', @(~, ~) choisir_dossier()), 1, 1);
placer(uibutton(g_dos, 'Text', char(8634), 'Tooltip', 'Back to the default folder of the overall setting', ...
    'ButtonPushedFcn', @(~, ~) dossier_defaut()), 1, 2);
placer(uilabel(gr, 'Text', 'Analysis', 'Tooltip', ...
    sprintf(['Where steps 2, 3, 3b and 4 are written and read.\n' ...
    'main: in the step 1 folder.\nOther: <step 1 folder>\\analyses\\<name>, on the same step 1\n' ...
    '(e.g. another K rule, without redoing step 1 or losing the main analysis).'])), 13, 1);
dd_ana = placer(uidropdown(gr, 'Items', {'main (in the step 1 folder)'}, 'ItemsData', 0, ...
    'ValueChangedFcn', @(src, ~) changer_analyse(src.Value)), 13, [2 3]);
placer(uibutton(gr, 'Text', 'New...', 'Tooltip', ...
    'New analysis on the same step 1: steps 2-4 in <step 1 folder>\analyses\<name>', ...
    'ButtonPushedFcn', @(~, ~) nouvelle_analyse()), 13, 4);

% Buttons
texte_arret = 'Stop after this session';
p_btn = uipanel(gauche, 'Title', 'Actions');
gb = uigridlayout(p_btn, [6 2]);
gb.RowHeight = repmat({24}, 1, 6);
gb.RowSpacing = 4;
b_tout   = uibutton(gb, 'Text', 'Analyse all sessions', 'ButtonPushedFcn', @(~, ~) lancer('tout'));
b_sel    = uibutton(gb, 'Text', 'Analyse selection', 'ButtonPushedFcn', @(~, ~) lancer('selection'));
b_ajout  = uibutton(gb, 'Text', 'Add sessions...', 'ButtonPushedFcn', @(~, ~) ajouter());
b_excel  = uibutton(gb, 'Text', 'Open table (Excel)', 'ButtonPushedFcn', @(~, ~) ouvrir_table());
b_doss   = uibutton(gb, 'Text', 'Open session folder', 'ButtonPushedFcn', @(~, ~) ouvrir_dossier());
b_rafr   = uibutton(gb, 'Text', 'Refresh', 'ButtonPushedFcn', @(~, ~) rafraichir_table());
b_diagk  = uibutton(gb, 'Text', 'k-means diagnostic', 'ButtonPushedFcn', @(~, ~) diagnostic_k(), ...
    'Tooltip', 'k-means replicates and online phase, on the selected session and the displayed K');
b_diagk.Layout.Row = 4; b_diagk.Layout.Column = 1;
b_arret  = uibutton(gb, 'Text', texte_arret, 'Enable', 'off', ...
    'FontColor', [0.7 0 0], 'ButtonPushedFcn', @(~, ~) demander_arret(), ...
    'Tooltip', 'The analysis stops after the current session');
b_arret.Layout.Row = 4; b_arret.Layout.Column = 2;
b_resume = uibutton(gb, 'Text', 'Results summary', 'ButtonPushedFcn', @(~, ~) voir_resume(), ...
    'Tooltip', 'Summary sentences of the last step 4 of this analysis (tests and numbers)');
b_resume.Layout.Row = 5; b_resume.Layout.Column = 1;
b_figs = uibutton(gb, 'Text', 'Results figures', 'ButtonPushedFcn', @(~, ~) voir_figures(), ...
    'Tooltip', 'Figures of the last steps 2, 3b and 4 of this analysis (images; Redraw = MATLAB figures)');
b_figs.Layout.Row = 5; b_figs.Layout.Column = 2;
b_param = uibutton(gb, 'Text', 'Analysis parameters', 'ButtonPushedFcn', @(~, ~) voir_parametres(), ...
    'Tooltip', 'All the parameters of this analysis, steps 1 to 4, read in the result files');
b_param.Layout.Row = 6; b_param.Layout.Column = 1;
% 22 nov. (demandé par l'utilisateur) : tableau de toutes les analyses, avec
% un code, leurs réglages et leurs résultats ; on y choisit l'analyse
b_registre = uibutton(gb, 'Text', 'Analyses table...', 'ButtonPushedFcn', @(~, ~) ouvrir_registre(), ...
    'Tooltip', sprintf(['All the analyses (no cues and cues) with a code, their date, settings,\n' ...
    'steps done, sessions, SCEs, assemblies and clonal bias. Choose one to display it.']));
b_registre.Layout.Row = 6; b_registre.Layout.Column = 2;
boutons = [b_tout, b_sel, b_ajout, b_excel, b_doss, b_rafr, b_diagk, b_resume, b_figs, b_param, b_registre];

% Progress
p_prog = uipanel(gauche, 'Title', 'Progress');
gp = uigridlayout(p_prog, [2 1]);
gp.RowHeight = {16, 18};
gp.RowSpacing = 2;
gp.Padding = [6 2 6 2];
lbl_etat = uilabel(gp, 'Text', 'No analysis running', 'FontSize', 11);
lbl_barre = uilabel(gp, 'Text', '', 'FontName', 'Consolas', 'FontSize', 12);
t_debut = [];

% Session table
tab = uitable(gauche, 'ColumnName', {'Session', 'Step 1', 'K', 'Assemblies', 'Clonal'}, ...
    'ColumnWidth', {200, 55, 45, 75, 50}, 'RowName', {}, ...
    'CellSelectionCallback', @selection_table);
try, tab.SelectionType = 'row'; catch, end

% Log
txt_log = uitextarea(gauche, 'Editable', 'off', 'FontName', 'Consolas', 'FontSize', 11);

%% ---------- Right part ----------
% Three rows: speed + activity raster (full width) | SCE raster, GCaMP
% image, one assembly in brainbow colours (with < > buttons) | session panel
droite = uigridlayout(g, [4 1]);
droite.RowHeight = {90, '1x', '1.35x', 195};   % speed, activity, 3 figures, session panel
droite.RowSpacing = 4;
droite.Padding = [0 0 0 0];
ax_vit = uiaxes(droite);
ax_vit.Layout.Row = 1; ax_vit.Layout.Column = 1;
ax_act = uiaxes(droite);
ax_act.Layout.Row = 2; ax_act.Layout.Column = 1;
zone_bas = uigridlayout(droite, [1 4]);            % SCE raster, field colours, one assembly, pies
zone_bas.Layout.Row = 3; zone_bas.Layout.Column = 1;
zone_bas.ColumnWidth = {'0.7x', '1x', '1x', 130};   % raster narrower (9 nov.) to make room for the pies
zone_bas.ColumnSpacing = 6;
zone_bas.Padding = [0 0 0 0];
ax_raster = uiaxes(zone_bas);
ax_raster.Layout.Row = 1; ax_raster.Layout.Column = 1;
ax_gcamp = uiaxes(zone_bas);
ax_gcamp.Layout.Row = 1; ax_gcamp.Layout.Column = 2;
g_coul = uigridlayout(zone_bas, [1 2], 'Padding', 0, 'ColumnSpacing', 2);   % image | assembly buttons on the side
g_coul.Layout.Row = 1; g_coul.Layout.Column = 3;
g_coul.ColumnWidth = {'1x', 34};
ax_coul = uiaxes(g_coul);
ax_coul.Layout.Row = 1; ax_coul.Layout.Column = 1;
g_nav = placer(uigridlayout(g_coul, [5 1], 'Padding', 0, 'RowSpacing', 3, ...
    'RowHeight', {'1x', 30, 22, 30, '1x'}), 1, 2);
placer(uibutton(g_nav, 'Text', char(9650), 'Tooltip', 'Previous assembly', ...
    'ButtonPushedFcn', @(~, ~) changer_asm_coul(-1)), 2, 1);
lbl_asm_coul = placer(uilabel(g_nav, 'Text', '', 'HorizontalAlignment', 'center', 'FontSize', 10), 3, 1);
placer(uibutton(g_nav, 'Text', char(9660), 'Tooltip', 'Next assembly', ...
    'ButtonPushedFcn', @(~, ~) changer_asm_coul(+1)), 4, 1);
% Pies (9 nov., request of the user): colours of the coloured cells of the
% field (middle image) and of the assembly of the 3rd image
g_pie = uigridlayout(zone_bas, [2 1], 'Padding', 0, 'RowSpacing', 4);
g_pie.Layout.Row = 1; g_pie.Layout.Column = 4;
ax_pie_champ = uiaxes(g_pie);
ax_pie_champ.Layout.Row = 1; ax_pie_champ.Layout.Column = 1;
ax_pie_asm = uiaxes(g_pie);
ax_pie_asm.Layout.Row = 2; ax_pie_asm.Layout.Column = 1;
for ax_p = [ax_pie_champ, ax_pie_asm]
    try, ax_p.Toolbar.Visible = 'off'; catch, end
    try, disableDefaultInteractivity(ax_p); catch, end
end
p_info = uipanel(droite, 'Title', 'Session');
p_info.Layout.Row = 4; p_info.Layout.Column = 1;
% Session panel: controls | results text | K choice
gi = uigridlayout(p_info, [6 4]);
gi.RowHeight = {22, 22, 22, 22, 22, 24};
gi.ColumnWidth = {70, 150, '1x', 220};   % controls | results text | K choice
gi.RowSpacing = 4;
gi.ColumnSpacing = 8;
placer(uilabel(gi, 'Text', 'Displayed K'), 1, 1);
g_K = placer(uigridlayout(gi, [1 3], 'Padding', 0, 'ColumnSpacing', 3, ...
    'ColumnWidth', {24, '1x', 24}), 1, 2);
placer(uibutton(g_K, 'Text', char(8722), 'FontWeight', 'bold', 'Tooltip', 'Previous K', ...
    'ButtonPushedFcn', @(~, ~) changer_K(-1)), 1, 1);
dd_K = placer(uidropdown(g_K, 'Items', {'-'}, 'ItemsData', NaN, 'ValueChangedFcn', @(~, ~) afficher()), 1, 2);
placer(uibutton(g_K, 'Text', '+', 'FontWeight', 'bold', 'Tooltip', 'Next K', ...
    'ButtonPushedFcn', @(~, ~) changer_K(+1)), 1, 3);
placer(uilabel(gi, 'Text', 'Assembly'), 2, 1);
dd_asm = placer(uidropdown(gi, 'Items', {'all'}, 'ItemsData', 0, 'ValueChangedFcn', @(~, ~) afficher(), ...
    'Tooltip', 'Show a single assembly in all views'), 2, 2);
cb_coul = placer(uicheckbox(gi, 'Text', 'Show brainbow colours', 'Enable', 'off', ...
    'Value', true, 'Tooltip', 'Available once step 3 (clonal) is done for the session', ...
    'ValueChangedFcn', @(~, ~) afficher()), 3, [1 2]);
cb_cont = placer(uicheckbox(gi, 'Text', 'Assembly outlines (brainbow)', 'Value', true, ...
    'ValueChangedFcn', @(~, ~) afficher()), 4, [1 2]);
cb_sce = placer(uicheckbox(gi, 'Text', 'Other SCEs (white lines)', 'Value', false, ...
    'Tooltip', 'Activity raster: a white line at each SCE where no assembly is significantly active', ...
    'ValueChangedFcn', @(~, ~) afficher()), 5, [1 2]);
b_img = placer(uibutton(gi, 'Text', 'Background...', 'Enable', 'off', ...
    'Tooltip', 'Colour background image for the assembly views', ...
    'ButtonPushedFcn', @(~, ~) choisir_image()), 6, 1);
placer(uibutton(gi, 'Text', 'Composite viewer...', ...
    'Tooltip', 'Brainbow composite image of the field, channel by channel, with the GCaMP masks and their numbers', ...
    'ButtonPushedFcn', @(~, ~) ouvrir_composite()), 6, 2);
txt_info = placer(uitextarea(gi, 'Editable', 'off', 'FontSize', 11), [1 6], 3);
ax_K = uiaxes(gi);
ax_K.Layout.Row = [1 6]; ax_K.Layout.Column = 4;

afficher_regle_sce();
vider_axes();
remplir_analyses('');
restaurer_choix();          % 22 nov. : dernière analyse choisie (session précédente)
verifier_fichiers();
journal('Interface ready. Select a session in the table to display it.');
rafraichir_table();


%% ======================= RUNS =========================================
    function lancer(mode)
        if en_cours, return; end
        cfg = config_interface();
        switch mode
            case 'tout'
                cfg.seulement = {};
                texte = 'all sessions';
            case 'selection'
                if isempty(selection) || isempty(sessions)
                    uialert(fig, 'Select one or more sessions in the table.', 'No session');
                    return
                end
                cfg.seulement = sessions.nom(selection)';
                texte = strjoin(cfg.seulement, ', ');
        end
        % Step 3b needs the new main_brainbow.m (an older copy ignores it silently)
        f_main = which('main_brainbow');
        if cfg.faire_etape3b
            texte_main = '';
            try, texte_main = fileread(f_main); catch, end
            if ~contains(texte_main, 'faire_etape3b') || isempty(which('etape3b_paires'))
                uialert(fig, sprintf(['Step 3b is ticked, but the main_brainbow.m used by MATLAB does ' ...
                    'not run it:\n%s\n\nReplace it with the new main_brainbow.m (and check that the ' ...
                    'folder 3b_paires is in the pipeline).'], f_main), 'Step 3b');
                return
            end
        end
        cfg.rappel_session = @rappel;
        cfg.rappel_progression = @progression;
        cfg.rappel_resume = @afficher_resume;     % summary sentences at the end of step 4
        debut_travail(sprintf('Analysis: %s', texte));
        if cfg.faire_etape5 && (isempty(which('etape5_sequences')) || ~contains(fileread(f_main), 'faire_etape5'))
            uialert(fig, sprintf(['Step 5 is ticked, but the main_brainbow.m used by MATLAB does not run it:\n' ...
                '%s\n\nReplace it with the new main_brainbow.m (and check that the folder 5_sequences is ' ...
                'in the pipeline).'], f_main), 'Step 5');
            return
        end
        etapes = {'1', '2', '3', '3b', '4', '5'};
        cochees = [cfg.faire_etape1, cfg.faire_etape2, cfg.faire_etape3, cfg.faire_etape3b, cfg.faire_etape4, ...
            cfg.faire_etape5];
        journal(sprintf('Steps: %s | step 1 folder: %s | analysis: %s | pipeline: %s', ...
            strjoin(etapes(cochees), ', '), cfg.PathEtape1, texte_analyse(), fileparts(f_main)));
        if cfg.faire_etape2
            journal(sprintf('Step 2: %s, K %d to %d', ...
                dd_regleK.Items{strcmp(dd_regleK.ItemsData, dd_regleK.Value)}, ef_kmin.Value, ef_kmax.Value));
        end
        if cfg.faire_etape1 && ~isempty(nom_analyse)
            journal(['Step 1 writes in the step 1 folder: it is shared by the main analysis ' ...
                'and all the analyses of this folder.']);
        end
        try
            main_brainbow(cfg);
            progression(1, 'Analysis finished');
            journal('Analysis finished.');
        catch err
            if strcmp(err.identifier, 'brainbow:arret')
                journal('Analysis stopped at your request.');
            else
                journal(['ERROR: ' err.message]);
            end
        end
        fin_travail();
    end

    function rappel(nom, etape, msg)
        % called by main_brainbow after each session (steps 1 and 3)
        journal(msg);
        try
            rafraichir_table();
            s = find(strcmp(sessions.nom, nom), 1);
            if ~isempty(s) && etape == 1
                selection = s;
                remplir_K(s);
                afficher();
            end
        catch err
            journal(['Display failed: ' err.message]);
        end
        drawnow;
        if arret_demande
            error('brainbow:arret', 'Stop requested');
        end
    end

    function afficher_resume(lignes, fichier)
        % summary sentences of step 4 in a separate window (text can be copied)
        journal(sprintf('Results summary: %s', fichier));
        f2 = uifigure('Name', ['Results summary: ' fichier], 'Position', [120 120 1050 560]);
        gg = uigridlayout(f2, [1 1]);
        uitextarea(gg, 'Value', lignes(:), 'Editable', 'off', 'FontSize', 13);
    end

    function voir_resume()
        % last summary of step 4 in the current analysis
        % + the summaries of step 3b (pairs) and step 5 (sequences) when there are some
        d = fullfile(dossier_analyse(), 'agregation');
        L = dir(fullfile(d, 'resume_resultats*.txt'));
        da = dossier_analyse();
        autres = {fullfile(da, 'paires_agregation', 'resume_paires.txt'), ...
            fullfile(da, 'soeurs_agregation', 'resume_soeurs.txt'), ...
            fullfile(da, 'soeurs_agregation', 'resume_coeurs.txt'), ...
            fullfile(da, 'hierarchique_agregation', 'resume_hierarchique.txt'), ...
            fullfile(da, 'clonal_tous_K', 'resume_clonal_tous_K.txt'), ...
            fullfile(da, 'paires_assemblees', 'resume_paires_assemblees.txt'), ...
            fullfile(da, 'diagnostic_SCE', 'resume_diagnostic_SCE.txt'), ...
            fullfile(da, 'sequences_agregation', 'resume_sequences.txt')};
        existe = cellfun(@isfile, autres);
        if isempty(L) && ~any(existe)
            uialert(fig, sprintf('No summary yet in\n%s\n(run step 4, 3b, 3d or 5).', d), 'Results summary');
            return
        end
        lignes = {};
        f = '';
        if any(existe), f = autres{find(existe, 1)}; end
        if ~isempty(L)
            [~, i] = max([L.datenum]);
            f = fullfile(d, L(i).name);
            lignes = lire_lignes(f);
        end
        for fx = autres(existe)
            if isfile(fx{1})
                if ~isempty(lignes), lignes = [lignes, {''}]; end
                lignes = [lignes, lire_lignes(fx{1})];
            end
        end
        afficher_resume(lignes, f);
    end

    function lignes = lire_lignes(f)
        lignes = {};
        fid = fopen(f, 'r', 'n', 'UTF-8');
        if fid < 0
            lignes = {['Cannot read ' f]};
            return
        end
        while true
            t = fgetl(fid);
            if ~ischar(t), break; end
            lignes{end+1} = t; %#ok<AGROW>
        end
        fclose(fid);
    end

    function voir_parametres()
        % all the parameters actually used (step 1 to step 4), read in the files
        if isempty(sessions)
            uialert(fig, 'No session in the table.', 'Analysis parameters');
            return
        end
        journal('Reading the parameters of this analysis (a few seconds)...');
        drawnow;
        f = fullfile(dossier_analyse(), 'parametres_analyse.txt');
        try
            lignes = parametres_analyse(sessions, config_courante(), f);
        catch err
            journal(['Parameters: ' err.message]);
            return
        end
        f2 = uifigure('Name', ['Analysis parameters: ' texte_analyse()], 'Position', [100 80 1150 820]);
        gg = uigridlayout(f2, [1 1]);
        uitextarea(gg, 'Value', lignes(:), 'Editable', 'off', 'FontName', 'Consolas', 'FontSize', 11);
        journal(sprintf('Analysis parameters: %s (lines starting with ≠ differ between sessions)', f));
    end

    function voir_figures()
        % images of the last steps 2, 3b and 4 of the current analysis, in a window
        da = dossier_analyse();
        fichiers = {}; noms = {}; dossiers_v = {};
        % step 4: one folder per variant, main test first
        L = dir(fullfile(da, 'agregation'));
        L = L([L.isdir] & ~ismember({L.name}, {'.', '..', 'controles_champ'}));
        cfg = config_courante();
        principal = [cfg.clonal.variantes{1, 1} '_' cfg.clonal.variantes{1, 2}];
        [~, ordre] = sort(~startsWith({L.name}, principal));
        L = L(ordre);
        for i = 1:numel(L)
            dv = fullfile(L(i).folder, L(i).name);
            P = dir(fullfile(dv, '*.png'));
            for j = 1:numel(P)
                fichiers{end+1} = fullfile(dv, P(j).name); %#ok<AGROW>
                noms{end+1} = sprintf('%s | %s', nom_court(L(i).name), P(j).name); %#ok<AGROW>
                dossiers_v{end+1} = dv; %#ok<AGROW>
            end
        end
        % step 2 and step 3b
        autres = {fullfile(da, 'agregation', 'controles_champ', '*.png'), 'field controls'; ...
            fullfile(da, '*.png'), 'K choice'; fullfile(da, 'paires_agregation', '*.png'), 'pairs'; ...
            fullfile(da, 'soeurs_agregation', '*.png'), 'sisters / clonal cores'; ...
            fullfile(da, 'hierarchique_agregation', '*.png'), 'hierarchical'; ...
            fullfile(da, 'clonal_tous_K', '*.png'), 'clonal bias, all K'; ...
            fullfile(da, 'paires_assemblees', '*.png'), 'same-colour pairs in assemblies'; ...
            fullfile(da, 'diagnostic_SCE', '*.png'), 'SCE diagnostic'; ...
            fullfile(da, 'sequences_agregation', '*.png'), 'sequences'};
        for i = 1:size(autres, 1)
            P = dir(autres{i, 1});
            for j = 1:numel(P)
                fichiers{end+1} = fullfile(P(j).folder, P(j).name); %#ok<AGROW>
                noms{end+1} = sprintf('%s | %s', autres{i, 2}, P(j).name); %#ok<AGROW>
                dossiers_v{end+1} = ''; %#ok<AGROW>
            end
        end
        if isempty(fichiers)
            uialert(fig, sprintf('No figure yet in\n%s\n(run step 4).', da), 'Results figures');
            return
        end
        f2 = uifigure('Name', ['Results figures: ' texte_analyse()], 'Position', [70 50 1550 920]);
        gg = uigridlayout(f2, [2 2]);
        gg.ColumnWidth = {340, '1x'};
        gg.RowHeight = {'1x', 28};
        lb = uilistbox(gg, 'Items', noms, 'ItemsData', 1:numel(noms), 'Value', 1);
        lb.Layout.Row = 1; lb.Layout.Column = 1;
        im = uiimage(gg, 'ImageSource', fichiers{1}, 'ScaleMethod', 'fit');
        im.Layout.Row = [1 2]; im.Layout.Column = 2;
        lb.ValueChangedFcn = @(src, ~) set(im, 'ImageSource', fichiers{src.Value});
        gbt = uigridlayout(gg, [1 3], 'Padding', 0, 'ColumnSpacing', 4);
        gbt.Layout.Row = 2; gbt.Layout.Column = 1;
        uibutton(gbt, 'Text', 'Open', 'Tooltip', 'Open this image in the system viewer', ...
            'ButtonPushedFcn', @(~, ~) ouvrir(fichiers{lb.Value}));
        uibutton(gbt, 'Text', 'Folder', 'Tooltip', 'Open the folder of this image', ...
            'ButtonPushedFcn', @(~, ~) ouvrir(fileparts(fichiers{lb.Value})));
        uibutton(gbt, 'Text', 'Redraw', 'Tooltip', ['Redraw the step 4 figures of this variant as ' ...
            'MATLAB figures (zoom, data tips), without recomputing'], ...
            'ButtonPushedFcn', @(~, ~) redessiner(dossiers_v{lb.Value}));
        journal(sprintf('Results figures: %d image(s) in %s', numel(fichiers), da));
    end

    function redessiner(dv)
        if isempty(dv)
            journal('Redraw: only for the step 4 figures (variant folders).');
            return
        end
        try
            redessiner_figures_etape4(dv, config_courante());
            journal(['Figures redrawn: ' dv]);
        catch err
            journal(['Redraw failed: ' err.message]);
        end
    end

    function demander_arret()
        arret_demande = true;
        b_arret.Text = 'Stop requested...';
        journal('Stop requested: the analysis will stop after the current session.');
    end

    function diagnostic_k()
        s = une_session();
        if isempty(s) || en_cours, return; end
        d = dossier_source(s);          % step 1 results (preprocessing)
        K = dd_K.Value;
        if ~isfile(fullfile(d, 'pretraitement.mat')) || isnan(K) || K == 0
            uialert(fig, ['Needs the preprocessing of the session and a K chosen in "Displayed K" ' ...
                '(step 1 done).'], 'k-means diagnostic');
            return
        end
        debut_travail(sprintf('k-means diagnostic: %s, K = %d (a few minutes)', sessions.nom{s}, K));
        try
            T = diagnostic_kmeans(d, K);
            f2 = uifigure('Name', sprintf('k-means replicates: %s, K = %d', sessions.nom{s}, K), ...
                'Position', [100 100 1250 260]);
            gg = uigridlayout(f2, [1 1]);
            uitable(gg, 'Data', T, 'RowName', {});
            journal('k-means diagnostic finished (estimated time per session in the command window).');
        catch err
            journal(['ERROR: ' err.message]);
        end
        fin_travail();
    end

    function ajouter()
        if en_cours, return; end
        [f, p] = uigetfile({'*.mat', 'data_recording files (*.mat)'}, ...
            'Sessions to add', 'MultiSelect', 'on');
        figure(fig);
        if isequal(f, 0), return; end
        if ischar(f), f = {f}; end
        try
            ajouter_sessions(config_courante(), fullfile(p, f));
            journal(sprintf('%d file(s) proposed; check the new rows in the Excel table.', numel(f)));
        catch err
            journal(['ERROR: ' err.message]);
        end
        rafraichir_table();
    end

    function changer_reglage()
        % overall setting: every field takes the values of
        % config_brainbow(setting), and the table shows its results
        if en_cours
            dd_reg.Value = config0.reglage;
            return
        end
        config0 = config_brainbow(dd_reg.Value);
        dossier_resultats = dossier_condition();   % default folder of this setting (and condition)
        ef_dossier.Value = dossier_resultats;
        ef_dossier.Tooltip = dossier_resultats;
        remplir_analyses('');                      % main analysis
        remplir_reglages();
        rafraichir_table();
        if ~isempty(selection), remplir_K(selection(end)); afficher(); end
        journal(sprintf('Overall setting: %s (step 1 folder %s)', dd_reg.Value, dossier_resultats));
    end

    function cfg = config_courante()
        % config_brainbow of the overall setting, with the folders chosen in
        % the interface (read again: takes your edits of the file into account)
        cfg = config_brainbow(dd_reg.Value);
        cfg.conditions = {condition_courante};     % condition chosen in the interface (19 nov.)
        cfg.PathSave = dossier_resultats;
        cfg.analyse = nom_analyse;
        % 22 nov. : la version des assemblées suit l'analyse choisie (sinon,
        % avec version_assemblees = 'v2' dans config_brainbow, une analyse
        % v1 était cherchée dans <nom>_v2 ou analyses\v2 : « K choice not
        % done », couleurs cachées). Nom en v2 (v2, v2_..., ..._v2) : v2 ;
        % autre nom : v1 ; principale : v1 si elle a déjà ses étapes 2-4,
        % sinon la version de config_brainbow (lancement : comme avant).
        if ~isempty(nom_analyse)
            if est_v2(fullfile(dossier_resultats, 'analyses', nom_analyse)), cfg.version_assemblees = 'v2';
            else, cfg.version_assemblees = 'v1'; end
        elseif isfile(fullfile(dossier_resultats, 'best_k_par_session.mat'))
            cfg.version_assemblees = 'v1';
        end
        if isfield(cfg, 'PathSource'), cfg = rmfield(cfg, 'PathSource'); end
        if isfield(cfg, 'PathEtape1'), cfg = rmfield(cfg, 'PathEtape1'); end
        cfg = dossiers_analyse(cfg);    % PathEtape1 = step 1, PathSave = steps 2-4
    end

    function ds = dossier_source(s)
        % folder of the session where the step 1 results are read
        if ismember('dossier_source', sessions.Properties.VariableNames)
            ds = sessions.dossier_source{s};
        else
            ds = sessions.dossier{s};
        end
    end

    function d = dossier_analyse()
        % folder of steps 2-4 for all sessions
        if isempty(nom_analyse)
            d = dossier_resultats;
        else
            d = fullfile(dossier_resultats, 'analyses', nom_analyse);
        end
    end

    function t = texte_analyse()
        if isempty(nom_analyse)
            t = 'main';
        else
            t = sprintf('%s (%s)', nom_analyse, dossier_analyse());
        end
    end

    function remplir_analyses(nom_voulu)
        % "Analysis" menu: main + the subfolders of <step 1 folder>\analyses;
        % selects nom_voulu if it exists in this folder, otherwise main
        L = dir(fullfile(dossier_resultats, 'analyses'));
        L = L([L.isdir]);
        noms = {L.name};
        noms_analyses = sort(noms(~ismember(noms, {'.', '..'})));
        i = find(strcmp(noms_analyses, nom_voulu), 1);
        if isempty(i), i = 0; end
        % 22 nov. : code de l'analyse (tableau des analyses) devant son nom
        textes = [{'main (in the step 1 folder)'}, noms_analyses];
        chemins = [{dossier_resultats}, cellfun(@(n) fullfile(dossier_resultats, 'analyses', n), ...
            noms_analyses, 'UniformOutput', false)];
        C = table();
        try, C = registre_analyses(racine_analyses(), 'codes'); catch, end
        for k = 1:numel(textes)
            if isempty(C), break; end
            j = find(strcmpi(C.chemin, chemins{k}), 1);
            if ~isempty(j), textes{k} = [C.code{j} '  ' textes{k}]; end
        end
        set(dd_ana, 'Items', textes, 'ItemsData', 0:numel(noms_analyses));
        dd_ana.Value = i;
        if i == 0, nom_analyse = ''; else, nom_analyse = noms_analyses{i}; end
        dd_ana.Tooltip = dossier_analyse();
        mettre_titre();
    end

    function c = code_analyse(chemin)
        % code A001... du tableau des analyses ('' si pas encore inventoriée)
        c = '';
        try
            C = registre_analyses(racine_analyses(), 'codes');
            i = find(strcmpi(C.chemin, chemin), 1);
            if ~isempty(i), c = C.code{i}; end
        catch
        end
    end

    function r = racine_analyses()
        % dossier qui contient nocues\ et cues\ (celui de config_brainbow)
        r = fileparts(fileparts(config0.PathSave));
    end

    function mettre_titre()
        c = code_analyse(dossier_analyse());
        t = 'Brainbow: assemblies and clonal analysis';
        if ~isempty(c), t = sprintf('%s  [%s]', t, c); end
        fig.Name = t;
    end

    function memoriser_choix()
        % 22 nov. : l'analyse choisie est gardée pour la prochaine ouverture
        try
            setpref('brainbow_interface', 'derniere_analyse', struct('condition', condition_courante, ...
                'dossier', dossier_resultats, 'analyse', nom_analyse));
        catch
        end
    end

    function restaurer_choix()
        % au démarrage : dernière analyse choisie, si son dossier existe encore
        if ~ispref('brainbow_interface', 'derniere_analyse'), return; end
        try
            ch = getpref('brainbow_interface', 'derniere_analyse');
            if ~isfolder(ch.dossier), return; end
            aller_a(ch.condition, ch.dossier, ch.analyse);
        catch err
            journal(['Last analysis not restored: ' err.message]);
        end
    end

    function aller_a(condition, dossier, nom)
        % affiche l'analyse <dossier>\(analyses\<nom>) de cette condition
        if ismember(condition, dd_cond.ItemsData)
            condition_courante = condition;
            dd_cond.Value = condition;
        end
        dossier_resultats = dossier;
        ef_dossier.Value = dossier;
        ef_dossier.Tooltip = dossier;
        remplir_analyses(nom);
        selection = [];
        rafraichir_table();
        journal(sprintf('Analysis %s: step 1 %s | %s', code_analyse(dossier_analyse()), dossier, texte_analyse()));
        reglage_de_l_analyse();
        memoriser_choix();
    end

    function ouvrir_registre()
        % 22 nov. (demandé par l'utilisateur) : tableau de toutes les analyses
        ecran2 = get(groot, 'ScreenSize');
        f2 = uifigure('Name', 'Analyses', 'Position', [20 60 min(1860, ecran2(3) - 40) min(820, ecran2(4) - 120)]);
        g2 = uigridlayout(f2, [4 1]);
        g2.RowHeight = {34, '1x', 60, 28};
        lbl2 = uilabel(g2, 'Text', 'Reading the analyses...', 'WordWrap', 'on');
        t2 = uitable(g2, 'RowName', {}, 'ColumnSortable', true);
        try, t2.SelectionType = 'row'; catch, end
        txt2 = uitextarea(g2, 'Editable', 'off', 'FontSize', 11, ...
            'Value', {'Select a row: its full settings and folders are shown here.'});
        gb2 = uigridlayout(g2, [1 5], 'Padding', 0, 'ColumnSpacing', 6);
        gb2.ColumnWidth = {130, 130, 130, '1x', 110};
        uibutton(gb2, 'Text', 'Choose', 'FontWeight', 'bold', 'ButtonPushedFcn', @(~, ~) choisir());
        uibutton(gb2, 'Text', 'Open folder', 'ButtonPushedFcn', @(~, ~) ouvrir());
        uibutton(gb2, 'Text', 'Rescan all', 'Tooltip', 'Read every analysis again (slower; otherwise only changed ones)', ...
            'ButtonPushedFcn', @(~, ~) remplir('rescanner'));
        uilabel(gb2, 'Text', '');
        uibutton(gb2, 'Text', 'Close', 'ButtonPushedFcn', @(~, ~) delete(f2));
        T2 = table();
        ligne_choisie = [];
        t2.CellSelectionCallback = @(~, e) noter(e);
        try, t2.DoubleClickedFcn = @(~, ~) choisir(); catch, end
        drawnow;
        remplir('');

        function remplir(mode)
            lbl2.Text = 'Reading the analyses...';
            drawnow;
            try
                T2 = registre_analyses(racine_analyses(), mode);
            catch err
                lbl2.Text = ['Error: ' err.message];
                return
            end
            vues = {'Code', 'Condition', 'Date', 'Step1_folder', 'Analysis', 'Version', ...
                'Steps', 'Sessions', 'Last_run', 'SCE', 'Assemblies', 'Testable', 'Pct_p05', ...
                'Pct_pminp05', 'N_q05', 'N_q10', 'Settings'};
            t2.Data = table2cell(T2(:, vues));     % cellules : ColumnFormat (22.2 et non 22.2000)
            t2.ColumnName = {'Code', 'Condition', 'Date', 'Step 1 folder', 'Analysis', 'Version', ...
                'Steps (sessions)', 'Sessions ok', 'Last run', 'SCEs', 'Assemblies', ...
                'Testable', '% p<0.05', '% p_minp<0.05', 'q<0.05', 'q<0.10', 'Settings'};
            t2.ColumnWidth = {45, 60, 105, 160, 200, 45, 110, 70, 150, 55, 75, 65, 70, 100, 58, 58, 'auto'};
            try
                t2.ColumnFormat = [repmat({'char'}, 1, 9), repmat({'shortG'}, 1, 7), {'char'}];
            catch
            end
            % analyse affichée en surbrillance
            i = find(strcmpi(T2.Analysis_path, dossier_analyse()), 1);
            try
                removeStyle(t2);
                if ~isempty(i), addStyle(t2, uistyle('BackgroundColor', [1 0.95 0.7]), 'row', i); end
            catch
            end
            lbl2.Text = sprintf(['%d analyses (newest code first; click a column title to sort). ' ...
                'Steps: number of sessions done at each step. Bias: field / global test, testable ' ...
                'assemblies; p = hypergeometric p of the most enriched colour (uncorrected); q = BH ' ...
                'FDR of step 4. Highlighted: the analysis shown. Double-click or Choose to display one.'], height(T2));
            ligne_choisie = [];
        end

        function noter(e)
            ligne_choisie = [];
            if isempty(e.Indices), return; end
            r = e.Indices(end, 1);
            % ligne affichée -> ligne du tableau (tri possible)
            try
                if ~isempty(t2.DisplayData), code = t2.DisplayData{r, 1}; else, code = t2.Data{r, 1}; end
                ligne_choisie = find(strcmp(T2.Code, code), 1);
            catch
                ligne_choisie = r;
            end
            if ~isempty(ligne_choisie) && ligne_choisie <= height(T2)
                L = T2(ligne_choisie, :);
                txt2.Value = {sprintf('%s  %s | %s | %s', L.Code{1}, L.Condition{1}, L.Step1_folder{1}, L.Analysis{1}); ...
                    ['Settings: ' L.Settings{1}]; ['Steps 2-4: ' L.Analysis_path{1}]};
            end
        end

        function choisir()
            if isempty(ligne_choisie) || ligne_choisie > height(T2)
                uialert(f2, 'Select an analysis in the table first.', 'Analyses');
                return
            end
            if en_cours
                uialert(f2, 'An analysis is running: wait until it ends.', 'Analyses');
                return
            end
            L = T2(ligne_choisie, :);
            cond = 'sans_indice';
            if strcmp(L.Condition{1}, 'cues'), cond = 'avec_indice'; end
            nom = L.Analysis{1};
            if strcmp(nom, 'main'), nom = ''; end
            aller_a(cond, L.Step1_path{1}, nom);
            figure(fig);
            delete(f2);
        end

        function ouvrir()
            if isempty(ligne_choisie), return; end
            winopen(T2.Analysis_path{ligne_choisie});
        end
    end

    function changer_analyse(v)
        if en_cours
            remplir_analyses(nom_analyse);      % no change during a run
            return
        end
        if v == 0, nom_analyse = ''; else, nom_analyse = noms_analyses{v}; end
        dd_ana.Tooltip = dossier_analyse();
        rafraichir_table();
        if ~isempty(selection), remplir_K(selection(end)); afficher(); end
        journal(sprintf('Analysis (steps 2-4): %s', texte_analyse()));
        reglage_de_l_analyse();
        mettre_titre();
        memoriser_choix();
    end

    function reglage_de_l_analyse()
        % the K choice and K range this analysis was made with (best_k_par_session.mat
        % of its last step 2) are shown in the settings
        f = fullfile(dossier_analyse(), 'best_k_par_session.mat');
        if ~isfile(f), return; end
        try
            B = load(f, 'regle_choix_K');
            r = B.regle_choix_K;
            regle = 'score';
            if isfield(r, 'regle'), regle = r.regle; end
            if ismember(regle, {'exces', 'consensus', 'score'}), dd_regleK.Value = regle; end
            if isfield(r, 'liste_K') && ~isempty(r.liste_K)
                kmin = min(r.liste_K);
                if isfield(r, 'k_start'), kmin = max(kmin, r.k_start); end
                ef_kmin.Value = kmin;
                ef_kmax.Value = max(r.liste_K);
            end
            journal(sprintf('This analysis was made with: %s, K %d to %d', ...
                dd_regleK.Items{strcmp(dd_regleK.ItemsData, dd_regleK.Value)}, ef_kmin.Value, ef_kmax.Value));
        catch
        end
    end

    function changer_regleK()
        % consensus: proposes a separate analysis, so that the current results are kept
        if en_cours, return; end
        if strcmp(dd_regleK.Value, 'consensus') && isempty(nom_analyse)
            choix = uiconfirm(fig, sprintf(['Step 2 will build the consensus of K %d to %d.\n\n' ...
                'Put its results (steps 2-4) in a separate analysis named "consensus"? ' ...
                'The current results (main) are then kept.'], ef_kmin.Value, ef_kmax.Value), ...
                'Consensus', 'Options', {'Yes, analysis "consensus"', 'No, replace main'}, ...
                'DefaultOption', 1, 'CancelOption', 1);
            if startsWith(choix, 'Yes')
                d = fullfile(dossier_resultats, 'analyses', 'consensus');
                if ~isfolder(d), mkdir(d); end
                remplir_analyses('consensus');
                dd_ana.Tooltip = dossier_analyse();
                rafraichir_table();
                if ~isempty(selection), remplir_K(selection(end)); afficher(); end
                journal(sprintf('Analysis (steps 2-4): %s', texte_analyse()));
            end
        end
        journal(sprintf('Step 2 K choice: %s, K %d to %d. Tick step 2 (and 3, 4), untick step 1, then run.', ...
            dd_regleK.Items{strcmp(dd_regleK.ItemsData, dd_regleK.Value)}, ef_kmin.Value, ef_kmax.Value));
    end

    function nouvelle_analyse()
        if en_cours, return; end
        r = inputdlg({sprintf(['Name of the new analysis.\nSteps 2-4 will be written in\n' ...
            '%s\\<name>,\nusing the step 1 of %s.'], fullfile(dossier_resultats, 'analyses'), ...
            dossier_resultats)}, 'New analysis', [1 70], {['analyse_' datestr(now, 'yyyy_mm_dd')]});
        figure(fig);
        if isempty(r), return; end
        nom = strtrim(r{1});
        if isempty(nom) || ~isempty(regexp(nom, '[\\/:*?"<>|]', 'once')) || any(strcmp(nom, {'.', '..'}))
            uialert(fig, 'Invalid name (no \ / : * ? " < > |).', 'New analysis');
            return
        end
        % 22 nov. : la version vient du nom ; avec version_assemblees = 'v2'
        % (config_brainbow), '_v2' est ajouté comme le fait dossiers_analyse
        if strcmpi(char(config0.version_assemblees), 'v2') && ...
                ~est_v2(fullfile(dossier_resultats, 'analyses', nom))
            nom = [nom '_v2'];
        end
        d = fullfile(dossier_resultats, 'analyses', nom);
        if ~isfolder(d), mkdir(d); end
        remplir_analyses(nom);
        changer_analyse(dd_ana.Value);
        journal('Tick steps 2, 3 (3b, 4) and run: step 1 is read in the step 1 folder.');
    end

    function verifier_fichiers()
        % the functions used must be those of this pipeline folder: an old
        % copy elsewhere on the path (or a file missing here) breaks the
        % folders silently
        utiles = {'main_brainbow', 'dossiers_analyse', 'lire_table_sessions', ...
            'analyse_assembly_results6', 'etape3_clonal', 'etape3b_paires', 'agregation_stats_clonales3', ...
            'etape5_sequences'};
        mauvais = {};
        for i = 1:numel(utiles)
            w = which(utiles{i});
            if isempty(w) || ~startsWith(lower(w), lower(dossier_pipeline))
                if isempty(w), w = 'missing'; end
                mauvais{end+1} = sprintf('%s: %s', utiles{i}, w); %#ok<AGROW>
            end
        end
        if ~isempty(mauvais)
            msg = sprintf(['These files do not come from the pipeline folder %s:\n%s\n\n' ...
                'Replace the whole folder with the content of the zip.'], dossier_pipeline, ...
                strjoin(mauvais, '\n'));
            journal(['WARNING: ' strrep(msg, sprintf('\n'), ' ')]);
            uialert(fig, msg, 'Old files', 'Icon', 'warning');
        end
    end

    function choisir_dossier()
        if en_cours, return; end
        depart = dossier_resultats;
        if ~isfolder(depart), depart = fileparts(depart); end
        if ~isfolder(depart), depart = pwd; end
        d = uigetdir(depart, 'Step 1 folder (one subfolder per session)');
        figure(fig);
        if isequal(d, 0), return; end
        changer_dossier(d);
    end

    function dossier_defaut()
        % default folder of the overall setting (config0 read at click time)
        changer_dossier(dossier_condition());
    end

    function d = dossier_condition()
        % default step 1 folder of the overall setting for the chosen
        % condition: config PathSave for no cues; for cues, the folder of the
        % same name in ...\cues (next to ...\nocues)
        d = config0.PathSave;
        if strcmp(condition_courante, 'avec_indice')
            racine_nocues = fileparts(d);
            d = fullfile(fileparts(racine_nocues), 'cues', nom_dossier(d));
        end
    end

    function n = nom_dossier(d)
        [~, n, e] = fileparts(regexprep(char(d), '[\\/]+$', ''));
        n = [n e];
    end

    function changer_condition(v)
        % 19 nov. : no cues / cues ; table, folder and analyses follow
        if en_cours
            dd_cond.Value = condition_courante;
            return
        end
        condition_courante = v;
        selection = [];
        journal(sprintf('Condition: %s', v));
        changer_dossier(dossier_condition());
    end

    function changer_dossier(d)
        % new step 1 folder: the table and the display show its results;
        % runs write there. The session table (Excel) does not change.
        if en_cours
            ef_dossier.Value = dossier_resultats;
            return
        end
        d = strtrim(char(d));
        if isempty(d), d = dossier_condition(); end
        if ~isfolder(d)
            choix = uiconfirm(fig, sprintf(['This folder does not exist:\n%s\nCreate it? ' ...
                '(step 1 must then be run there)'], d), ...
                'Step 1 folder', 'Options', {'Create', 'Cancel'}, 'DefaultOption', 1, 'CancelOption', 2);
            if ~strcmp(choix, 'Create')
                ef_dossier.Value = dossier_resultats;
                return
            end
            mkdir(d);
        end
        dossier_resultats = d;
        ef_dossier.Value = d;
        ef_dossier.Tooltip = d;
        remplir_analyses(nom_analyse);      % same analysis if it exists here, otherwise main
        rafraichir_table();
        if ~isempty(selection), remplir_K(selection(end)); afficher(); end
        journal(sprintf('Step 1 folder: %s | analysis: %s', d, texte_analyse()));
        memoriser_choix();
        sans_etape1 = isempty(dir(fullfile(d, '*', 'pretraitement.mat')));
        if sans_etape1 && ~isempty(dir(fullfile(d, '*', 'choix_K.mat')))
            % steps 2-3 without step 1: a separate analysis folder (old way)
            msg = sprintf(['This folder has K choices but no step 1: it is an analysis folder.\n\n' ...
                'Move it into <step 1 folder>\\analyses\\ (the step 1 it was computed on), ' ...
                'then choose that step 1 folder here and the analysis in the "Analysis" menu.']);
            journal(['WARNING: ' strrep(msg, sprintf('\n'), ' ')]);
            uialert(fig, msg, 'Analysis folder', 'Icon', 'warning');
        elseif sans_etape1
            journal('No step 1 result in this folder yet.');
        end
    end

    function remplir_reglages()
        c = config0;
        regle = 'percentile';
        if isfield(c.pre, 'sce_regle') && ismember(c.pre.sce_regle, {'percentile', 'fixe', 'fdr'})
            regle = c.pre.sce_regle;
        end
        seuils_sce = struct('percentile', 99, 'fixe', 15, 'fdr', 0.05);
        if isfield(c.pre, 'sce_percentile'), seuils_sce.percentile = c.pre.sce_percentile; end
        if isfield(c.pre, 'sce_fdr'), seuils_sce.fdr = c.pre.sce_fdr; end
        if strcmp(regle, 'fixe'), seuils_sce.fixe = c.pre.sce_n_cells_threshold; end
        regle_affichee = regle;
        dd_sce.Value = regle;
        afficher_regle_sce();
        if isfield(c.asm, 'seuil_clusters'), dd_cl.Value = c.asm.seuil_clusters; end
        if isfield(c.asm, 'percentile_seuil'), ef_pct.Value = c.asm.percentile_seuil; else, ef_pct.Value = 95; end
        ef_kmin.Value = min(c.liste_K);
        ef_kmax.Value = max(c.liste_K);
        if isfield(c.choixK, 'regle') && ismember(c.choixK.regle, {'exces', 'consensus', 'score'})
            dd_regleK.Value = c.choixK.regle;
        end
        ef_rep.Value  = c.asm.kmean_iter;
        if isfield(c.asm, 'kmean_iter_surrogate')
            ef_rep_sur.Value = c.asm.kmean_iter_surrogate;
        else
            ef_rep_sur.Value = c.asm.kmean_iter;
        end
        ef_nsur.Value = c.asm.kmeans_surrogate;
        ef_fdr.Value  = c.asm.FDR_Q;
        cb_e1.Value = c.faire_etape1;
        cb_e2.Value = c.faire_etape2;
        cb_e3.Value = c.faire_etape3;
        cb_e3b.Value = valeur_3b(c);
        cb_e4.Value = c.faire_etape4;
        cb_e5.Value = valeur_5(c);
        cb_forcer.Value = c.forcer_pretraitement || c.forcer_assemblees;
        cb_fig.Value = c.faire_figures;
    end

    function changer_regle_sce()
        % keeps the value of the rule left, shows the value of the new one
        seuils_sce.(regle_affichee) = ef_sce.Value;
        regle_affichee = dd_sce.Value;
        afficher_regle_sce();
    end

    function afficher_regle_sce()
        switch regle_affichee
            case 'percentile'
                texte = 'percentile:'; limites = [50 99.9]; arrondi = 'off';
                aide = 'Percentile of the shuffles (99 recommended; 95 = lower threshold)';
            case 'fixe'
                texte = 'cells:'; limites = [1 500]; arrondi = 'on';
                aide = 'Minimum number of co-active cells for an SCE';
            case 'fdr'
                texte = 'SCE FDR:'; limites = [0.001 0.5]; arrondi = 'off';
                aide = sprintf(['Maximum share of SCEs expected by chance (0.05 = 5 %%).\n' ...
                    'Chance: shifts of +/- sce_decalage_fdr_s (config_brainbow) during rest']);
        end
        lbl_sce.Text = texte;
        v = min(max(seuils_sce.(regle_affichee), limites(1)), limites(2));
        set(ef_sce, 'Limits', [-Inf Inf], 'RoundFractionalValues', 'off');
        set(ef_sce, 'Value', v);
        set(ef_sce, 'Limits', limites, 'RoundFractionalValues', arrondi, 'Tooltip', aide);
    end

    function cfg = config_interface()
        cfg = config_courante();
        seuils_sce.(regle_affichee) = ef_sce.Value;
        cfg.pre.sce_regle      = dd_sce.Value;
        cfg.pre.sce_percentile = seuils_sce.percentile;
        cfg.pre.sce_fdr        = seuils_sce.fdr;
        if strcmp(dd_sce.Value, 'fixe')
            cfg.pre.sce_n_cells_threshold = seuils_sce.fixe;
        elseif ~isfield(cfg.pre, 'sce_minimum') || isempty(cfg.pre.sce_minimum)
            % règle calculée : minimum séparé du seuil fixe (sce_minimum, 8 nov.)
            cfg.pre.sce_minimum = minimum_SCE_calcule;
        end
        cfg.asm.seuil_clusters = dd_cl.Value;
        cfg.asm.FDR_Q = ef_fdr.Value;
        kmin = min(ef_kmin.Value, ef_kmax.Value);
        kmax = max(ef_kmin.Value, ef_kmax.Value);
        cfg.liste_K = kmin:kmax;
        cfg.asm.kmean_iter           = ef_rep.Value;
        cfg.asm.kmeans_surrogate     = ef_nsur.Value;
        cfg.asm.kmean_iter_surrogate = ef_rep_sur.Value;
        cfg.asm.percentile_seuil     = ef_pct.Value;
        cfg.choixK.k_start = kmin;
        cfg.choixK.regle = dd_regleK.Value;
        cfg.faire_etape1 = cb_e1.Value;
        cfg.faire_etape2 = cb_e2.Value;
        cfg.faire_etape3 = cb_e3.Value;
        cfg.faire_etape3b = cb_e3b.Value;
        cfg.faire_etape4 = cb_e4.Value;
        cfg.faire_etape5 = cb_e5.Value;
        cfg.forcer_pretraitement = cb_forcer.Value;
        cfg.forcer_assemblees    = cb_forcer.Value;
        cfg.faire_figures        = cb_fig.Value;
    end

    function progression(fraction, message)
        % Progress bar: called by main_brainbow (session, K, step)
        fraction = min(max(fraction, 0), 1);
        n = 34;
        plein = round(fraction * n);
        barre = [repmat(char(9608), 1, plein), repmat(char(9617), 1, n - plein)];
        reste = '';
        if fraction > 0.02 && fraction < 1 && ~isempty(t_debut)
            ecoule = toc(t_debut);
            reste = ['  ~' duree_texte(ecoule * (1 - fraction) / fraction) ' left'];
        end
        lbl_barre.Text = sprintf('%s %3.0f %%%s', barre, 100 * fraction, reste);
        lbl_etat.Text = message;
        drawnow limitrate;
    end

    function t = duree_texte(secondes)
        if secondes < 90
            t = sprintf('%.0f s', secondes);
        elseif secondes < 3600
            t = sprintf('%.0f min', secondes / 60);
        else
            t = sprintf('%d h %02.0f', floor(secondes / 3600), mod(secondes, 3600) / 60);
        end
    end

    function debut_travail(texte)
        en_cours = true;
        arret_demande = false;
        set(boutons, 'Enable', 'off');
        b_arret.Enable = 'on';
        b_arret.Text = texte_arret;
        t_debut = tic;
        progression(0, texte);
        journal(texte);
        drawnow;
    end

    function fin_travail()
        en_cours = false;
        set(boutons, 'Enable', 'on');
        b_arret.Enable = 'off';
        b_arret.Text = texte_arret;
        rafraichir_table();
        if ~isempty(selection), afficher(); end
    end


%% ======================= SESSION TABLE ================================
    function rafraichir_table()
        try
            sessions = lire_table_sessions(config_courante());
        catch err
            sessions = table();
            tab.Data = cell(0, 5);
            journal(['Session table: ' err.message]);
            return
        end
        n = height(sessions);
        D = cell(n, 5);
        for s = 1:n
            D(s, :) = etat_session(s);
        end
        tab.Data = D;
        selection = selection(selection <= n);
    end

    function ligne = etat_session(s)
        d = sessions.dossier{s};        % results of steps 2-4
        ds = dossier_source(s);          % step 1 results
        e1 = '-'; K = ''; nasm = ''; cl = '';
        if isfile(fullfile(ds, 'pretraitement.mat'))
            e1 = sprintf('%d K', numel(K_disponibles(ds)));
        end
        f = fullfile(d, 'choix_K.mat');
        if isfile(f)
            C = load(f, 'best_K', 'statut_choix');
            if strcmp(C.statut_choix, 'ok')
                K = texte_K(C.best_K);
                fa = fichier_asm(s, C.best_K);
                if isfile(fa)
                    A = load(fa, 'assemblyortho');
                    nasm = sprintf('%d', numel(A.assemblyortho));
                end
            else
                K = 'excluded';
            end
        end
        if isfile(fullfile(d, 'clonal', 'stats_clonales.mat')), cl = 'done'; end
        ligne = {sessions.nom{s}, e1, K, nasm, cl};
    end

    function selection_table(~, evt)
        if isempty(evt.Indices), return; end
        selection = unique(evt.Indices(:, 1))';
        s = selection(end);
        try
            remplir_K(s);
        catch err
            journal(sprintf('K list of %s: %s', sessions.nom{s}, err.message));
        end
        afficher();
    end

    function s = une_session()
        s = [];
        if isempty(selection) || isempty(sessions)
            uialert(fig, 'Select a session in the table.', 'No session');
            return
        end
        s = selection(end);
    end

    function f = fichier_asm(s, K)
        % assemblies of the session for K (0 = consensus of the K, in the analysis folder)
        f = fichier_assemblees(sessions.dossier{s}, dossier_source(s), K);   % v1 ou v2 selon l'analyse (11 nov.)
    end

    function Ks = K_disponibles(d)
        L = dir(fullfile(d, 'k*'));
        L = L([L.isdir]);
        Ks = str2double(regexprep({L.name}, '^k', ''));
        Ks = sort(Ks(~isnan(Ks)));
        Ks = Ks(arrayfun(@(k) isfile(fullfile(d, sprintf('k%d', k), 'assemblees.mat')), Ks));
    end

    function remplir_K(s)
        d = sessions.dossier{s};
        Ks = K_disponibles(dossier_source(s));
        best = NaN;
        f = fullfile(d, 'choix_K.mat');
        if isfile(f)
            C = load(f, 'best_K', 'statut_choix');
            if strcmp(C.statut_choix, 'ok'), best = C.best_K; end
        end
        if isfile(fullfile(d, 'consensus', 'assemblees.mat'))
            Ks = [0, Ks];                   % 0 = consensus of the K (step 2)
        end
        if isempty(Ks)
            set(dd_K, 'Items', {'-'}, 'ItemsData', NaN);
            return
        end
        noms = arrayfun(@(k) texte_K(k, true), Ks, 'UniformOutput', false);
        noms(Ks == best) = strcat(noms(Ks == best), ' (chosen)');
        set(dd_K, 'Items', noms, 'ItemsData', Ks);
        if any(Ks == best), dd_K.Value = best; else, dd_K.Value = Ks(1); end
    end

    function changer_K(pas)
        % - and + buttons: previous or next K among those computed
        Ks = dd_K.ItemsData;
        if isempty(Ks) || all(isnan(Ks)), return; end
        i = find(Ks == dd_K.Value, 1);
        if isempty(i), return; end
        i = i + pas;
        if i < 1 || i > numel(Ks), return; end
        dd_K.Value = Ks(i);
        afficher();
    end


%% ======================= DISPLAY ======================================
    function afficher()
        % Display of the selected session. An error is written in the panel
        % and in the log (with the file and line) instead of stopping silently.
        try
            afficher_session();
        catch err
            if isempty(selection) || isempty(sessions), return; end
            s = selection(end);
            ou = '';
            if ~isempty(err.stack), ou = sprintf(' (%s, line %d)', err.stack(1).name, err.stack(1).line); end
            journal(sprintf('Display of %s failed: %s%s', sessions.nom{s}, err.message, ou));
            txt_info.Value = {sessions.nom{s}; ['Display failed: ' err.message ou]; ''; ...
                ['Step 1 folder: ' dossier_source(s)]; ['Steps 2-4: ' sessions.dossier{s}]};
        end
    end

    function afficher_session()
        % Unlink the two images BEFORE clearing them: otherwise the message
        % written in one of them (limits 0-1) imposes its limits on the other,
        % and the GCaMP image shrinks to a black corner.
        try
            linkaxes([ax_gcamp, ax_coul], 'off');
        catch
        end
        vider_axes();
        if isempty(selection) || isempty(sessions), return; end
        s = selection(end);
        d = sessions.dossier{s};        % results of steps 2-4 (K choice, clonal)
        ds = dossier_source(s);          % step 1 results (preprocessing, assemblies of each K)
        info = {sessions.nom{s}};
        dossiers_info = {''; ['Step 1: ' ds]};
        if ~strcmp(fileparts(d), fileparts(ds))
            dossiers_info{end+1} = ['Steps 2-4 (analysis ' nom_analyse '): ' d];
        end
        fpre = fullfile(ds, 'pretraitement.mat');
        if ~isfile(fpre)
            txt_info.Value = [info; {'No step 1 result for this session in the step 1 folder.'}; dossiers_info];
            message_axes(ax_act, {'No step 1 result for this session'});
            return
        end
        P = load(fpre, 'Race', 'opts');
        Race = double(P.Race);
        [NCell, NRace] = size(Race);
        seuil = NaN;
        if isfield(P.opts, 'sce_n_cells_threshold_calcule'), seuil = P.opts.sce_n_cells_threshold_calcule; end
        info{end+1} = sprintf('%d cells, %d SCEs (threshold %g cells)', NCell, NRace, seuil);

        % Assemblies of the displayed K
        K = dd_K.Value;
        A = struct('assemblyortho', {{}});
        if ~isnan(K)
            fa = fichier_asm(s, K);
            noms_var = intersect({'assemblyortho', 'PCl_final', 'x1', 'RList', 'NClOK', 'sCl', ...
                'threshold', 'IDX2', 'cluster_origine', 'K_presence', 'K_utilises'}, who('-file', fa));
            A = load(fa, noms_var{:});
        end
        if ~isfield(A, 'assemblyortho'), A.assemblyortho = {}; end
        nA = numel(A.assemblyortho);
        cols = palette(nA);
        % assembly(ies) shown: all, or a single one ("Assembly" menu)
        ancienne = dd_asm.Value;
        set(dd_asm, 'Items', [{'all'}, arrayfun(@(c) sprintf('%d', c), 1:nA, 'UniformOutput', false)], ...
            'ItemsData', 0:nA);
        if ancienne <= nA, dd_asm.Value = ancienne; else, dd_asm.Value = 0; end
        if dd_asm.Value == 0, quelles = 1:nA; else, quelles = dd_asm.Value; end
        n_sce_asm = NaN;
        if isfield(A, 'PCl_final') && size(A.PCl_final, 2) == NRace
            n_sce_asm = sum(any(A.PCl_final > 0, 1));
        end
        info{end+1} = sprintf('%s: %d assembl%s, %d of %d SCEs with an active assembly', ...
            texte_K(K, true), nA, pluriel(nA), n_sce_asm, NRace);
        if K == 0 && isfield(A, 'K_presence') && numel(A.K_presence) == nA
            info{end+1} = sprintf('Consensus of K %s; each assembly seen at: %s', ...
                mat2str(A.K_utilises), strjoin(cellfun(@(k) mat2str(k), A.K_presence, ...
                'UniformOutput', false), ' | '));
        end
        if ~isfile(fullfile(d, 'choix_K.mat'))
            info{end+1} = 'K choice: not done in this analysis yet (step 2).';
        end
        % Clonal results (only once step 3 is done): p_minp and biased colour
        % on each assembly's line
        fcl = fullfile(d, 'clonal', 'stats_clonales.mat');
        clonal_fait = isfile(fcl);
        suffixes = repmat({''}, 1, nA);
        lignes_clonal = {};
        if clonal_fait
            try
                [suffixes, lignes_clonal] = resume_clonal(fcl, d, K, nA);
            catch err
                lignes_clonal = {['Clonal: ' err.message]};
            end
        end
        % other K (8 nov.): p and p_minp from main_clonal_tous_K, or computed
        % here for this session and K (only once step 3 is done: colours)
        if nA > 0 && all(cellfun(@isempty, suffixes))
            try
                st = stats_assemblees(s, d, K, nA, A, clonal_fait);
                if ~isempty(st)
                    for c = 1:nA
                        if st.testable(c)
                            etoile = '';
                            if st.p_minp(c) < 0.05, etoile = ' *'; end
                            suffixes{c} = sprintf(' | p = %s, p_minp = %s (%s)%s', texte_p_court(st.p_initial(c)), ...
                                texte_p_court(st.p_minp(c)), st.couleur{c}, etoile);
                        else
                            suffixes{c} = ' | not testable';
                        end
                    end
                    lignes_clonal = {'', sprintf(['p: hypergeometric p of the most enriched colour (uncorrected); ' ...
                        'p_minp: corrected for the 6 colours (%s); * p_minp < 0.05, before FDR. Exploratory at K ' ...
                        'other than the one chosen at step 2.'], st.source)};
                end
            catch err
                lignes_clonal = {['Clonal (all K): ' err.message]};
            end
        end
        for c = 1:nA
            n_act = NaN;
            if isfield(A, 'PCl_final') && size(A.PCl_final, 1) >= c
                n_act = sum(A.PCl_final(c, :) > 0);
            end
            info{end+1} = sprintf('  %d: %d cells, %d activations%s', c, numel(A.assemblyortho{c}), ...
                n_act, suffixes{c}); %#ok<AGROW>
        end

        dessiner_raster(Race, A, cols, K, quelles);
        try
            dessiner_activite(ds, A, cols, NCell, quelles, char(sessions.dossier_suite2p{s}));
        catch err
            message_axes(ax_act, {'Activity raster unavailable:', err.message});
        end
        dessiner_courbe_K(d);

        % Geometry (GCaMP image and masks: no brainbow colour)
        G = [];
        try
            G = geometrie(sessions.dossier_suite2p{s});
        catch err
            info{end+1} = ['GCaMP image: ' err.message];
        end
        if ~isempty(G) && G.n ~= NCell
            info{end+1} = sprintf('WARNING: %d cell ROIs in Fall.mat, %d in the data.', G.n, NCell);
            G = [];
        end
        % brainbow colours: only once step 3 is done, and if the box is ticked
        colorcell = [];
        err_coul = '';
        if clonal_fait && cb_coul.Value && ~isempty(G)
            try
                colorcell = charger_couleurs(avec_separateur(sessions.dossier_suite2p{s}), ...
                    texte_cellule(sessions.fichier_colorcell, s));
                colorcell = couleurs_sans_doublons(colorcell, struct('dossier', d, 'dossier_source', ds));
                colorcell = colorcell(:);
            catch err
                err_coul = err.message;
            end
        end

        % Middle image (8 nov., request of the user): every cell of the field
        % filled with its brainbow colour. While colours are hidden (before
        % step 3, or box unticked): the assembly cells, as before.
        if ~isempty(G)
            if ~isempty(colorcell) && numel(colorcell) == G.n
                dessiner_champ_couleur(G, colorcell, char(sessions.dossier_suite2p{s}));
            else
                dessiner_gcamp(G, A, cols, quelles);
            end
        end

        % 3rd image (8 nov., request of the user): ONE assembly at a time on
        % the background image, each cell filled with its brainbow colour,
        % with its p (uncorrected), p_minp and the colour in excess; < > to
        % change assembly. Colours only once step 3 is done (shown straight
        % away, unless the box was unticked).
        couleurs_dessinees = false;
        cb_coul.Enable = onoff(clonal_fait);
        b_img.Enable = onoff(clonal_fait);
        vue_coul = struct('G', G, 'A', A, 'nA', nA, 'cc', [], 'stats', [], ...
            'cle_img', char(sessions.dossier_suite2p{s}), 'K', K);
        cle = sprintf('%s|%g', sessions.nom{s}, K);
        if ~strcmp(cle, cle_vue_coul), idx_coul = 1; cle_vue_coul = cle; end
        if dd_asm.Value > 0, idx_coul = dd_asm.Value; end
        lbl_asm_coul.Text = '';
        if ~clonal_fait
            message_axes(ax_coul, {'Colours hidden', 'until step 3 (clonal)'});
        elseif ~cb_coul.Value
            message_axes(ax_coul, {'Tick "Show brainbow colours"'});
        elseif isempty(G)
            message_axes(ax_coul, {'GCaMP image / masks unavailable'});
        elseif nA == 0
            message_axes(ax_coul, {'No assembly at this K'});
        elseif ~isempty(err_coul)
            message_axes(ax_coul, {'Colours unavailable:', err_coul});
        else
            try
                vue_coul.cc = colorcell;
                vue_coul.stats = stats_assemblees(s, d, K, nA, A, clonal_fait);
                couleurs_dessinees = dessiner_assemblee_couleur();
            catch err
                message_axes(ax_coul, {'Colours unavailable:', err.message});
            end
        end
        if ~isempty(lignes_clonal)
            info = [info, lignes_clonal];
        end
        txt_info.Value = [info(:); dossiers_info];
        if ~isempty(G) && couleurs_dessinees      % same zoom, only if both are images
            try
                linkaxes([ax_gcamp, ax_coul], 'xy');
            catch
            end
        end
    end

    function t = pluriel(n)
        if n == 1, t = 'y'; else, t = 'ies'; end
    end

    function dessiner_raster(Race, A, cols, K, quelles)
        [NCell, NRace] = size(Race);
        ordre_c = 1:NCell; ordre_s = 1:NRace;
        if isfield(A, 'x1') && numel(A.x1) == NCell, ordre_c = A.x1(:)'; end
        if isfield(A, 'RList') && numel(A.RList) == NRace, ordre_s = A.RList(:)'; end
        actif = Race > 0;
        img = repmat(0.55 * double(actif), 1, 1, 3);        % grey on black
        for c = quelles
            if ~isfield(A, 'PCl_final') || size(A.PCl_final, 1) < c || size(A.PCl_final, 2) ~= NRace
                continue
            end
            m = false(NCell, NRace);
            m(A.assemblyortho{c}, A.PCl_final(c, :) > 0) = true;
            m = m & actif;
            for k = 1:3
                plan = img(:, :, k);
                plan(m) = cols(c, k);
                img(:, :, k) = plan;
            end
        end
        image(ax_raster, img(ordre_c, ordre_s, :));
        ax_raster.XLim = [0.5 NRace + 0.5];
        ax_raster.YLim = [0.5 NCell + 0.5];
        ax_raster.YDir = 'reverse';
        ax_raster.Color = 'k';
        xlabel(ax_raster, 'SCEs (sorted)');
        ylabel(ax_raster, 'Cells (sorted)');
        title(ax_raster, sprintf('SCEs by assembly (%s)', texte_K(K, true)), 'FontWeight', 'normal');
    end

    function dessiner_activite(d, A, cols, NCell, quelles, dossier_s2p)
        % Activity raster (cells x time, whole session) sorted like the SCE
        % raster (assemblies at the bottom), with the speed above.
        act = activite_session(d, dossier_s2p);
        if act.NCell ~= NCell
            message_axes(ax_act, {'Different number of cells between files'});
            return
        end
        ordre_c = 1:NCell;
        if isfield(A, 'x1') && numel(A.x1) == NCell, ordre_c = A.x1(:)'; end
        rang = zeros(NCell, 1);
        rang(ordre_c) = 1:NCell;
        t = act.ti / act.fs;
        y = rang(act.ci);
        % Coloured: event of an assembly cell DURING an SCE where this
        % assembly is significantly active (SCE window: TRace to
        % TRace + synchronous_frames, as for Race). The rest in grey.
        couleur_ev = zeros(numel(act.ti), 1);
        if isfield(A, 'PCl_final') && size(A.PCl_final, 2) == numel(act.TRace)
            for c = quelles(quelles <= size(A.PCl_final, 1))
                debut = act.TRace(A.PCl_final(c, :) > 0);
                fenetre = false(1, act.Nz);
                for k = 0:act.sync
                    idx = debut + k;
                    fenetre(idx(idx >= 1 & idx <= act.Nz)) = true;
                end
                membre = false(NCell, 1);
                cel = A.assemblyortho{c};
                membre(cel(cel >= 1 & cel <= NCell)) = true;
                dans = membre(act.ci) & fenetre(act.ti)';
                couleur_ev(dans & couleur_ev == 0) = c;
            end
        end
        cla(ax_act, 'reset');
        hold(ax_act, 'on');
        % Optional: a white line (full height, behind the points) at each SCE
        % where no assembly is significantly active
        if cb_sce.Value && ~isempty(act.TRace)
            autres = true(1, numel(act.TRace));
            if isfield(A, 'PCl_final') && size(A.PCl_final, 2) == numel(act.TRace) && ~isempty(A.PCl_final)
                autres = ~any(A.PCl_final > 0, 1);
            end
            ts = (act.TRace(autres) + act.sync / 2) / act.fs;      % middle of the SCE window
            X = [ts; ts; nan(1, numel(ts))];
            Y = repmat([0.5; NCell + 0.5; NaN], 1, numel(ts));
            plot(ax_act, X(:), Y(:), '-', 'Color', [1 1 1 0.4], 'LineWidth', 0.5);
        end
        m = couleur_ev == 0;
        plot(ax_act, t(m), y(m), '.', 'Color', [0.55 0.55 0.55], 'MarkerSize', taille_point_gris);
        for c = quelles
            m = couleur_ev == c;
            if any(m)
                plot(ax_act, t(m), y(m), '.', 'Color', cols(c, :), 'MarkerSize', taille_point_couleur);
            end
        end
        hold(ax_act, 'off');
        ax_act.Color = 'k';
        ax_act.YDir = 'reverse';
        ax_act.YLim = [0.5 NCell + 0.5];
        ax_act.XLim = [0 act.Nz / act.fs];
        ylabel(ax_act, 'Cells (sorted)');
        xlabel(ax_act, 'Time (s)');
        title(ax_act, sprintf('Activity: %s (coloured: during the SCEs of their assembly)', act.source), ...
            'FontWeight', 'normal');

        % Speed, and SCE times (triangle in the colour of the active assembly)
        ax_vit.Visible = 'on';
        if isempty(act.vitesse) || all(isnan(act.vitesse))
            message_axes(ax_vit, {'No speed'});
        else
            tv = (1:numel(act.vitesse)) / act.fs;
            plot(ax_vit, tv, act.vitesse, '-', 'Color', [0.25 0.25 0.25]);
            vmax = max(act.vitesse(:), [], 'omitnan');
            vmin = min(act.vitesse(:), [], 'omitnan');
            haut = vmax + 0.15 * max(vmax - vmin, eps);
            hold(ax_vit, 'on');
            if ~isempty(act.TRace)
                plot(ax_vit, act.TRace / act.fs, repmat(haut, 1, numel(act.TRace)), 'v', ...
                    'Color', [0.75 0.75 0.75], 'MarkerSize', 3);
            end
            if isfield(A, 'PCl_final') && size(A.PCl_final, 2) == numel(act.TRace)
                for c = quelles(quelles <= size(A.PCl_final, 1))
                    j = find(A.PCl_final(c, :) > 0);
                    plot(ax_vit, act.TRace(j) / act.fs, repmat(haut, 1, numel(j)), 'v', ...
                        'Color', cols(c, :), 'MarkerFaceColor', cols(c, :), 'MarkerSize', 4);
                end
            end
            hold(ax_vit, 'off');
            ax_vit.YLim = [vmin, haut + 0.1 * max(vmax - vmin, eps)];
            ax_vit.XLim = ax_act.XLim;
            ax_vit.XTickLabel = {};
            ylabel(ax_vit, 'Speed');
        end
        try
            linkaxes([ax_vit, ax_act], 'x');
        catch
        end
    end

    function act = activite_session(d, dossier_s2p)
        % Peaks of all transients (running included), with each cell's
        % threshold used by the detection; kept in memory.
        % 20 nov. (request of the user): CASCADE discrete spikes have no
        % threshold (th_detection = NaN), so findpeaks gave an empty raster.
        % For opts.detection = 'cascade_discrets': frames with >= 1 spike, read
        % from the spike cache of the suite2p folder (whole session, running
        % included); without the cache: Raster of step 1 (rest only).
        if nargin < 2, dossier_s2p = ''; end
        if isKey(cache_act, d)
            act = cache_act(d);
            return
        end
        f = fullfile(d, 'pretraitement.mat');
        noms_var = intersect({'Tr1b', 'th_detection', 'opts', 'speedsm', 'speed', 'sampling_rate', ...
            'TRace', 'bad_frames'}, who('-file', f));
        P = load(f, noms_var{:});
        [NC, Nz] = size(P.Tr1b);
        mpd = 2;
        if isfield(P, 'opts') && isfield(P.opts, 'MinPeakDistance'), mpd = P.opts.MinPeakDistance; end
        sync = 2;
        if isfield(P, 'opts') && isfield(P.opts, 'synchronous_frames'), sync = P.opts.synchronous_frames; end
        cascade = isfield(P, 'opts') && isfield(P.opts, 'detection') && ...
            strcmpi(P.opts.detection, 'cascade_discrets');
        source = 'transient peaks';
        if cascade
            modele = 'GC8_EXC_10Hz_smoothing150ms_high_noise';
            if isfield(P.opts, 'cascade_modele'), modele = P.opts.cascade_modele; end
            n_pa = [];
            fc = fullfile(dossier_s2p, ['cascade_discrets_' modele '.mat']);
            if ~isempty(dossier_s2p) && isfile(fc)
                C = load(fc, 'n_pa');
                if isequal(size(C.n_pa), [NC, Nz]), n_pa = C.n_pa; end
            end
            if ~isempty(n_pa)
                [ci, ti] = find(n_pa > 0);
                source = 'CASCADE spikes (frames with >= 1 spike)';
            else
                R = load(f, 'Raster');
                [ci, ti] = find(R.Raster > 0);
                source = 'CASCADE spikes, rest only (spike cache not found)';
            end
            ci = ci(:);
            ti = ti(:);
        else
            ci = cell(NC, 1);
            ti = cell(NC, 1);
            for i = 1:NC
                if ~(P.th_detection(i) > 0), continue; end
                [~, locs] = findpeaks(double(P.Tr1b(i, :)), 'MinPeakProminence', P.th_detection(i), ...
                    'MinPeakDistance', mpd);
                ti{i} = locs(:);
                ci{i} = repmat(i, numel(locs), 1);
            end
            ti = vertcat(ti{:});
            ci = vertcat(ci{:});
        end
        if isfield(P, 'bad_frames') && numel(P.bad_frames) == Nz && ~isempty(ti)
            ok = ~logical(P.bad_frames(ti));
            ti = ti(ok);
            ci = ci(ok);
        end
        vitesse = [];
        if isfield(P, 'speedsm'), vitesse = double(P.speedsm(:)');
        elseif isfield(P, 'speed'), vitesse = double(P.speed(:)'); end
        TRace = [];
        if isfield(P, 'TRace'), TRace = double(P.TRace(:)'); end
        act = struct('ci', ci, 'ti', ti, 'fs', P.sampling_rate, 'Nz', Nz, 'NCell', NC);
        act.vitesse = vitesse;
        act.TRace = TRace;
        act.sync = sync;
        act.source = source;
        cache_act(d) = act;
    end

    function dessiner_gcamp(G, A, cols, quelles)
        % Cells of the assemblies shown: filled with the colour of their
        % assembly. The others: thin grey outline.
        n = G.n;
        couleur = zeros(n, 3);
        alpha = zeros(n, 1);
        membre = false(n, 1);
        for c = quelles
            cel = cellules_valides(A.assemblyortho{c}, n);
            couleur(cel, :) = repmat(cols(c, :), numel(cel), 1);
            alpha(cel) = 0.8;
            membre(cel) = true;
        end
        image(ax_gcamp, superposer(G, couleur, alpha, 1));
        hold(ax_gcamp, 'on');
        contours(ax_gcamp, G, find(~membre), [0.45 0.45 0.45], 0.5);
        hold(ax_gcamp, 'off');
        image_axes(ax_gcamp, G);
        title(ax_gcamp, 'GCaMP: assembly cells', 'FontWeight', 'normal');
    end

    function dessiner_champ_couleur(G, cc, cle_img)
        % Middle image (8 nov., request of the user): every cell of the
        % field filled with its brainbow colour (grey: no colour, including
        % duplicates and ROIs removed for their size), nothing outlined, on
        % the same background as the 3rd image
        n = G.n;
        cc(~(cc >= 1 & cc <= 8)) = 8;
        couleur = rgb_brainbow();
        couleur = couleur(cc, :);
        alpha = 0.85 * ones(n, 1);
        if isKey(images_couleur, cle_img)
            image(ax_gcamp, remplir(images_couleur(cle_img), G, couleur, alpha));
        else
            image(ax_gcamp, superposer(G, couleur, alpha, 0.6));
        end
        image_axes(ax_gcamp, G);
        ax_gcamp.Box = 'on';
        title(ax_gcamp, sprintf('Field: brainbow colours (%d of %d cells coloured)', sum(cc <= 6), n), ...
            'FontWeight', 'normal', 'FontSize', 9);
        camembert_couleurs(ax_pie_champ, accumarray(cc(cc <= 6), 1, [6 1]), 'Field');
    end

    function camembert_couleurs(ax, nb, titre)
        % Pie of the coloured cells (R G B RG RB GB), drawn with patches;
        % percentage written in the slices of at least 10 %
        cla(ax, 'reset');
        cacher_axes(ax);
        hold(ax, 'on');
        nb = nb(:)';
        tot = sum(nb);
        rgb = rgb_brainbow();
        noms = {'R', 'G', 'B', 'RG', 'RB', 'GB'};
        if tot > 0
            a0 = pi / 2;                                  % from the top, clockwise
            for k = 1:6
                if nb(k) == 0, continue; end
                a1 = a0 - 2 * pi * nb(k) / tot;
                t = linspace(a0, a1, max(3, ceil(60 * nb(k) / tot)));
                patch(ax, [0, cos(t), 0], [0, sin(t), 0], rgb(k, :), 'EdgeColor', 'w', 'LineWidth', 0.5);
                if nb(k) / tot >= 0.1
                    am = (a0 + a1) / 2;
                    text(ax, 0.62 * cos(am), 0.62 * sin(am), sprintf('%s\n%d%%', noms{k}, round(100 * nb(k) / tot)), ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 8);
                end
                a0 = a1;
            end
        end
        hold(ax, 'off');
        ax.XLim = [-1.05 1.05]; ax.YLim = [-1.05 1.05];
        ax.DataAspectRatio = [1 1 1];
        title(ax, sprintf('%s: %d col.', titre, tot), 'FontWeight', 'normal', 'FontSize', 9);
    end

    function cacher_axes(ax)
        % no axes, ticks or background (the title stays visible)
        ax.XColor = 'none'; ax.YColor = 'none'; ax.Color = 'none';
        ax.XTick = []; ax.YTick = [];
    end

    function rgb = rgb_brainbow()
        % display colours of R G B RG RB GB, white (7), none (8)
        rgb = [1 0.15 0.15; 0.15 1 0.15; 0.25 0.5 1; 1 1 0.1; 1 0.2 1; 0.1 1 1; 1 1 1; 0.6 0.6 0.6];
    end

    function changer_asm_coul(pas)
        % < > under the 3rd image: previous / next assembly
        if ~isfield(vue_coul, 'nA') || vue_coul.nA == 0 || isempty(vue_coul.cc), return; end
        idx_coul = mod(idx_coul - 1 + pas, vue_coul.nA) + 1;
        try
            dessiner_assemblee_couleur();
        catch err
            message_axes(ax_coul, {'Colours unavailable:', err.message});
        end
    end

    function ok = dessiner_assemblee_couleur()
        % One assembly on the background image (GCaMP, or the colour image
        % chosen with "Colour background image..."): each member filled
        % with its brainbow colour, nothing outlined (8 nov., request of the
        % user). If p_minp < 0.05: frame of the image in the colour in
        % excess. Title on one line: p of the most enriched colour
        % (uncorrected), p_minp (corrected for the 6 colours).
        ok = false;
        G = vue_coul.G; A = vue_coul.A; n = G.n;
        cc = vue_coul.cc;
        if numel(cc) ~= n
            message_axes(ax_coul, {sprintf('colorcell: %d values for %d cells', numel(cc), n)});
            return
        end
        c = min(max(idx_coul, 1), vue_coul.nA);
        idx_coul = c;
        rgb_bb = rgb_brainbow();
        noms_bb = {'R', 'G', 'B', 'RG', 'RB', 'GB'};
        cc(~(cc >= 1 & cc <= 8)) = 8;
        mem = cellules_valides(A.assemblyortho{c}, n);
        st = vue_coul.stats;
        exces = '';
        if ~isempty(st) && numel(st.p_minp) >= c, exces = st.couleur{c}; end
        k_exc = find(strcmp(noms_bb, exces));
        couleur = zeros(n, 3);
        alpha = zeros(n, 1);
        couleur(mem, :) = rgb_bb(cc(mem), :);
        alpha(mem) = 0.85;
        if isKey(images_couleur, vue_coul.cle_img)
            image(ax_coul, remplir(images_couleur(vue_coul.cle_img), G, couleur, alpha));
        else
            image(ax_coul, superposer(G, couleur, alpha, 0.6));
        end
        image_axes(ax_coul, G);
        % frame in the colour in excess when p_minp < 0.05 (before FDR)
        signif = ~isempty(st) && numel(st.p_minp) >= c && st.testable(c) && st.p_minp(c) < 0.05 && ~isempty(k_exc);
        ax_coul.Box = 'on';
        if signif
            ax_coul.XColor = min(rgb_bb(k_exc, :), 0.9); ax_coul.YColor = ax_coul.XColor;
            ax_coul.LineWidth = 4;
        else
            ax_coul.XColor = [0.15 0.15 0.15]; ax_coul.YColor = [0.15 0.15 0.15];
            ax_coul.LineWidth = 0.5;
        end
        n_col = sum(cc(mem) <= 6);
        n_exc = 0;
        if ~isempty(k_exc), n_exc = sum(cc(mem) == k_exc); end
        if isempty(st) || numel(st.p_minp) < c
            fin = 'p not computed';
        elseif ~st.testable(c)
            fin = sprintf('%s %d/%d, not testable', exces, n_exc, n_col);
        else
            fin = sprintf('%s %d/%d, p = %s, p_{minp} = %s', exces, n_exc, n_col, ...
                texte_p_court(st.p_initial(c)), texte_p_court(st.p_minp(c)));
        end
        title(ax_coul, sprintf('%d cells (%d col.) | %s', numel(mem), n_col, fin), 'Interpreter', 'tex', ...
            'FontWeight', 'normal', 'FontSize', 9);
        lbl_asm_coul.Text = sprintf('%d/%d', c, vue_coul.nA);
        ccm = cc(mem);
        camembert_couleurs(ax_pie_asm, accumarray(ccm(ccm <= 6), 1, [6 1]), sprintf('Assembly %d', c));
        ok = true;
    end

    function t = texte_p_court(p)
        if isnan(p), t = 'NA'; elseif p < 0.001, t = '< 0.001'; else, t = sprintf('%.3f', p); end
    end

    function st = stats_assemblees(s, d, K, nA, A, calculer)
        % p (uncorrected), p_minp, colour in excess and testability of each
        % assembly at K: step 3 if it used this K, otherwise main_clonal_tous_K,
        % otherwise computed here (same test as step 3, main variant) if
        % calculer (step 3 done: colours allowed); [] otherwise
        st = [];
        if nA == 0, return; end
        cle_p = sprintf('%s|%s|%g', nom_analyse, sessions.nom{s}, K);
        if isKey(cache_p, cle_p)
            st = cache_p(cle_p);
            return
        end
        fcl = fullfile(d, 'clonal', 'stats_clonales.mat');
        if isfile(fcl)
            noms_var = intersect({'res', 'id_K'}, who('-file', fcl));
            C = load(fcl, noms_var{:});
            if isfield(C, 'id_K') && isequal(C.id_K, K) && isfield(C, 'res') && istable(C.res) && height(C.res) >= nA
                st = stats_depuis(C.res);
                st.source = 'step 3';
                cache_p(cle_p) = st;
                return
            end
        end
        if isempty(nom_analyse)
            base = dossier_resultats;
        else
            base = fullfile(dossier_resultats, 'analyses', nom_analyse);
        end
        ft = fullfile(base, 'clonal_tous_K', 'clonal_tous_K.mat');
        if isfile(ft)
            S = load(ft, 'T');
            R = S.T(strcmp(S.T.nom_session, sessions.nom{s}) & S.T.K == K, :);
            if height(R) >= nA
                st = stats_depuis(sortrows(R, 'assemblee'));
                st.source = 'main_clonal_tous_K';
                cache_p(cle_p) = st;
                return
            end
        end
        if ~calculer, return; end
        st = calculer_p(s, d, K, A);
        if ~isempty(st), cache_p(cle_p) = st; end
    end

    function st = calculer_p(s, d, K, A)
        % same test as step 3 (main variant of config.clonal.variantes), for
        % the session and K displayed; the slow descriptors (compactness) are
        % skipped. A few seconds.
        st = [];
        journal(sprintf('%s, %s: computing p and p_minp of each assembly...', sessions.nom{s}, texte_K(K, true)));
        drawnow;
        cfg = config_courante();
        ds = dossier_source(s);
        pre = load(fullfile(ds, 'pretraitement.mat'), 'Race');
        Race = double(pre.Race);
        path = avec_separateur(sessions.dossier_suite2p{s});
        [outline_gcampx, ~, gcamp_centroid] = charger_geometrie(path);
        [colorcell, ~, source_couleur] = charger_couleurs(path, texte_cellule(sessions.fichier_colorcell, s));
        colorcell = couleurs_sans_doublons(colorcell, struct('dossier', d, 'dossier_source', ds));
        verifier_alignement(size(Race, 1), colorcell, outline_gcampx, sessions.nom{s});
        um_px = sessions.um_par_px(s);
        if isnan(um_px), um_px = cfg.um_par_px_defaut; end
        taille_bloc_px = cfg.clonal.taille_bloc_px;
        if ~isnan(cfg.clonal.taille_bloc_um), taille_bloc_px = cfg.clonal.taille_bloc_um / um_px; end
        marge_um = 10;
        if isfield(cfg.clonal, 'marge_territoire_um'), marge_um = cfg.clonal.marge_territoire_um; end
        D = struct('marge_territoire_px', marge_um / um_px, 'colorcell', colorcell, ...
            'assembly', {A.assemblyortho}, 'gcamp_centroid', gcamp_centroid, 'Race', Race, ...
            'PCl_final', [], 'score_sequence', nan(numel(A.assemblyortho), 1), ...
            'id_animal', sessions.id_animal{s}, 'id_champ', '', 'id_session', sessions.nom{s}, ...
            'id_condition', '', 'id_K', K, 'source_couleur', source_couleur, ...
            'taille_bloc_px', taille_bloc_px, 'activite', sum(Race > 0, 2), 'p_detection', NaN);
        if isfield(A, 'PCl_final'), D.PCl_final = A.PCl_final; end
        opts = cfg.clonal;
        opts.descripteurs = false;
        variante = cfg.clonal.variantes(1, :);
        res = stats_clonales_champ3(D, opts, variante{1}, variante{2}, sessions.graine(s) + 10);
        st = stats_depuis(res);
        st.source = sprintf('computed here, %s_%s', variante{1}, variante{2});
        journal(sprintf('%s, %s: p and p_minp computed.', sessions.nom{s}, texte_K(K, true)));
    end

    function st = stats_depuis(R)
        st = struct();
        if ismember('p_initial', R.Properties.VariableNames)
            st.p_initial = R.p_initial;
        elseif ismember('p_brut', R.Properties.VariableNames)
            st.p_initial = min(R.p_brut, [], 2);
        else
            st.p_initial = nan(height(R), 1);
        end
        st.p_minp = R.p_minp;
        st.couleur = cellstr(string(R.couleur_dominante));
        st.testable = true(height(R), 1);
        if ismember('testable', R.Properties.VariableNames), st.testable = logical(R.testable); end
    end

    function cel = cellules_valides(cel, n)
        cel = cel(:);
        cel = cel(cel >= 1 & cel <= n);
    end

    function dessiner_courbe_K(d)
        f = fullfile(d, 'choix_K.mat');
        if ~isfile(f), return; end
        noms_var = intersect({'best_K', 'tableau_K', 'statut_choix', 'regle_choix_K'}, who('-file', f));
        C = load(f, noms_var{:});
        T = C.tableau_K;
        if ~istable(T) || isempty(T) || ~ismember('K', T.Properties.VariableNames)
            return
        end
        % plot the criterion that chose K: excess silhouette (rule 'exces') or the old score
        colonne = 'Score'; nom_y = 'Score';
        if isfield(C, 'regle_choix_K') && (~isfield(C.regle_choix_K, 'regle') || ...
                ~strcmp(C.regle_choix_K.regle, 'score')) && ismember('Exces', T.Properties.VariableNames)
            colonne = 'Exces'; nom_y = 'Excess silhouette';
        end
        if ~ismember(colonne, T.Properties.VariableNames), return; end
        plot(ax_K, T.K, T.(colonne), '-o', 'Color', [0.2 0.4 0.8], 'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerSize', 4);
        if strcmp(C.statut_choix, 'ok') && C.best_K > 0
            xline(ax_K, C.best_K, '--r');
        end
        xlabel(ax_K, 'K'); ylabel(ax_K, nom_y);
        title(ax_K, 'K choice', 'FontWeight', 'normal');
        ax_K.Visible = 'on';
    end

    function [suffixes, lignes] = resume_clonal(fcl, d, K, nA)
        % For each assembly: coloured cells, p_minp and most enriched colour
        % (main test), then the p_minp of the other tests; and the field
        % test. Only for the K chosen at step 2 (the clonal step uses it).
        suffixes = repmat({''}, 1, nA);
        lignes = {};
        C = load(fcl);
        if ~isfield(C, 'res') || ~istable(C.res) || isempty(C.res)
            if isfield(C, 'statut'), lignes = {['Clonal: ' char(string(C.statut))]}; end
            return
        end
        best = NaN;
        if isfile(fullfile(d, 'choix_K.mat'))
            B = load(fullfile(d, 'choix_K.mat'), 'best_K');
            best = B.best_K;
        end
        if ~isequal(K, best)
            lignes = {sprintf('Clonal: computed for %s only', texte_K(best, true))};
            return
        end
        liste = {C.res};
        noms = {'main'};
        if isfield(C, 'res_variantes') && isfield(C, 'noms_variantes') && ~isempty(C.res_variantes)
            liste = C.res_variantes;
            noms = C.noms_variantes;
        end
        R = liste{1};
        testable = true(height(R), 1);
        if ismember('testable', R.Properties.VariableNames), testable = logical(R.testable); end
        for i = 1:min(nA, height(R))
            n_col = NaN;
            if ismember('n_colorees', R.Properties.VariableNames), n_col = R.n_colorees(i); end
            if ~testable(i)
                txt = sprintf(' | %d coloured: not testable', n_col);
            else
                etoile = '';
                if R.p_minp(i) < 0.05, etoile = ' *'; end
                if ismember('p_initial', R.Properties.VariableNames)
                    p_init = R.p_initial(i);
                else
                    p_init = min(R.p_brut(i, :), [], 2);     % anciens résultats de l'étape 3
                end
                txt = sprintf(' | %d coloured: p = %s, p_minp = %.4f (%s)%s', n_col, texte_p_court(p_init), ...
                    R.p_minp(i), char(string(R.couleur_dominante(i))), etoile);
            end
            autres = {};
            for v = 2:numel(liste)
                Rv = liste{v};
                if height(Rv) < i, continue; end
                if ismember('testable', Rv.Properties.VariableNames) && ~Rv.testable(i), continue; end
                autres{end+1} = sprintf('%s %.4f', nom_court(noms{v}), Rv.p_minp(i)); %#ok<AGROW>
            end
            if ~isempty(autres)
                txt = [txt '   [' strjoin(autres, ', ') ']']; %#ok<AGROW>
            end
            suffixes{i} = txt;
        end
        lignes{end+1} = '';
        lignes{end+1} = sprintf(['p: hypergeometric p of the most enriched colour (uncorrected); p_minp: main ' ...
            'test (%s), corrected for the 6 colours; * p_minp < 0.05. [ ]: other tests (blocks = spatial ' ...
            'control). Per session, before the FDR of step 4.'], nom_court(noms{1}));
    end

    function t = nom_court(nom)
        switch nom
            case 'champ_globale',       t = 'field, global';
            case 'champ_blocs',         t = 'blocks';
            case 'champ_activite',      t = 'activity';
            case 'champ_territoire',    t = 'territory';
            case 'assemblees_blocs',    t = 'blocks (assemblies)';
            case 'assemblees_activite', t = 'activity';
            case 'assemblees_globale',  t = 'main';
            otherwise,                  t = nom;
        end
    end

    function ouvrir_composite()
        % separate window (10 nov., request of the user): visionneuse_composite
        if isempty(selection) || isempty(sessions), return; end
        s = selection(end);
        try
            visionneuse_composite(sessions.dossier_suite2p{s}, sessions.nom{s}, ...
                texte_cellule(sessions.fichier_colorcell, s));
        catch err
            uialert(fig, err.message, 'Composite viewer');
        end
    end

    function choisir_image()
        if isempty(selection) || isempty(sessions), return; end
        s = selection(end);
        [f, p] = uigetfile({'*.tif;*.tiff;*.png;*.jpg', 'Images'}, 'Colour image of this field', ...
            avec_separateur(sessions.dossier_suite2p{s}));
        figure(fig);
        if isequal(f, 0), return; end
        I = imread(fullfile(p, f));
        if size(I, 3) == 1, I = repmat(I, 1, 1, 3); end
        I = double(I(:, :, 1:3));
        for k = 1:3          % contrast per channel
            plan = I(:, :, k);
            lo = prctile(plan(:), 1); hi = prctile(plan(:), 99.7);
            I(:, :, k) = min(max((plan - lo) / max(hi - lo, eps), 0), 1);
        end
        images_couleur(char(sessions.dossier_suite2p{s})) = I;
        cb_coul.Value = true;
        afficher();
    end


%% ======================= TOOLS ========================================
    function G = geometrie(dossier_s2p)
        cle = char(dossier_s2p);
        if isKey(cache_geo, cle)
            G = cache_geo(cle);
            return
        end
        f = fullfile(cle, 'Fall.mat');
        S2P = load(f, 'stat', 'iscell');
        roi = find(S2P.iscell(:, 1) == 1);
        img = [];
        if isfile(fullfile(cle, 'calcium.tif'))
            img = imread(fullfile(cle, 'calcium.tif'));
        else
            try
                O = load(f, 'ops');
                img = O.ops.meanImg;
            catch
            end
        end
        n = numel(roi);
        xs = cell(1, n); ys = cell(1, n); cx = cell(1, n); cy = cell(1, n);
        for m = 1:n
            if iscell(S2P.stat), st = S2P.stat{roi(m)}; else, st = S2P.stat(roi(m)); end
            xs{m} = double(st.xpix(:)) + 1;  ys{m} = double(st.ypix(:)) + 1;
            cx{m} = double(st.xext(:)) + 1;  cy{m} = double(st.yext(:)) + 1;
        end
        if isempty(img)
            img = zeros(max(cellfun(@max, ys)), max(cellfun(@max, xs)));
        end
        img = double(img(:, :, 1));
        [Ly, Lx] = size(img);
        L = zeros(Ly, Lx, 'uint16');
        for m = 1:n
            ok = xs{m} >= 1 & xs{m} <= Lx & ys{m} >= 1 & ys{m} <= Ly;
            L(sub2ind([Ly Lx], ys{m}(ok), xs{m}(ok))) = m;
        end
        centre = zeros(n, 2);
        for m = 1:n
            centre(m, :) = [mean(xs{m}), mean(ys{m})];
        end
        G = struct('img', img, 'L', L, 'cx', {cx}, 'cy', {cy}, 'n', n, 'centre', centre);
        cache_geo(cle) = G;
    end

    function rgb = remplir(I, G, couleur, alpha)
        % colour background image (0-1) with the cells filled (same as superposer)
        rgb = double(I);
        if max(rgb(:)) > 1, rgb = rgb / max(rgb(:)); end
        if size(rgb, 1) ~= size(G.L, 1) || size(rgb, 2) ~= size(G.L, 2)    % same size as the masks
            iy = round(linspace(1, size(rgb, 1), size(G.L, 1)));
            ix = round(linspace(1, size(rgb, 2), size(G.L, 2)));
            rgb = rgb(iy, ix, :);
        end
        dedans = G.L > 0;
        idx = double(G.L(dedans));
        a = alpha(idx);
        for k = 1:3
            plan = rgb(:, :, k);
            plan(dedans) = plan(dedans) .* (1 - a) + couleur(idx, k) .* a;
            rgb(:, :, k) = plan;
        end
    end

    function rgb = superposer(G, couleur, alpha, luminosite)
        b = G.img;
        lo = prctile(b(:), 1); hi = prctile(b(:), 99.7);
        b = luminosite * min(max((b - lo) / max(hi - lo, eps), 0), 1);
        rgb = repmat(b, 1, 1, 3);
        dedans = G.L > 0;
        idx = double(G.L(dedans));
        a = alpha(idx);
        for k = 1:3
            plan = rgb(:, :, k);
            plan(dedans) = plan(dedans) .* (1 - a) + couleur(idx, k) .* a;
            rgb(:, :, k) = plan;
        end
    end

    function contours(ax, G, quelles, couleur, epaisseur)
        if isempty(quelles), return; end
        X = cellfun(@(x) [x; x(1); NaN], G.cx(quelles), 'UniformOutput', false);
        Y = cellfun(@(y) [y; y(1); NaN], G.cy(quelles), 'UniformOutput', false);
        plot(ax, vertcat(X{:}), vertcat(Y{:}), '-', 'Color', couleur, 'LineWidth', epaisseur);
    end

    function image_axes(ax, G)
        ax.XLim = [0.5 size(G.L, 2) + 0.5];
        ax.YLim = [0.5 size(G.L, 1) + 0.5];
        ax.YDir = 'reverse';
        ax.DataAspectRatio = [1 1 1];
        ax.XTick = []; ax.YTick = [];
    end

    function vider_axes()
        for ax = [ax_vit, ax_act, ax_raster, ax_gcamp, ax_coul, ax_K, ax_pie_champ, ax_pie_asm]
            cla(ax, 'reset');
        end
        cacher_axes(ax_pie_champ);
        cacher_axes(ax_pie_asm);
        message_axes(ax_act, {'Select a session'});
        ax_vit.Visible = 'off';
        ax_K.Visible = 'off';
        txt_info.Value = {''};
    end

    function message_axes(ax, lignes)
        cla(ax);
        text(ax, 0.5, 0.5, lignes, 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
            'Color', [0.5 0.5 0.5], 'FontSize', 12);
        ax.XTick = []; ax.YTick = [];
        ax.XLim = [0 1]; ax.YLim = [0 1];
        ax.Box = 'on';
    end

    function cols = palette(n)
        % Assembly colours (utilitaires/couleurs_assemblees.m): bright hues
        % between the brainbow colours, the same as in the PNG figures
        cols = couleurs_assemblees(n, 'noir');
    end

    function t = texte_cellule(colonne, s)
        t = '';
        if iscell(colonne) && ischar(colonne{s}), t = colonne{s}; end
    end

    function journal(msg)
        txt_log.Value = [{sprintf('%s  %s', datestr(now, 'HH:MM:SS'), msg)}; txt_log.Value(:)];
        drawnow limitrate;
    end

    function ouvrir_table()
        cfg = config_courante();
        ouvrir(cfg.fichier_sessions);
    end

    function ouvrir_dossier()
        s = une_session();
        if isempty(s), return; end
        d = sessions.dossier{s};                 % steps 2-4 of the analysis
        if ~isfolder(d), d = dossier_source(s); end   % otherwise step 1
        if ~isfolder(d)
            uialert(fig, 'This folder does not exist yet (session not analysed yet).', 'Folder');
            return
        end
        ouvrir(d);
    end

    function ouvrir(chemin)
        if ~(isfile(chemin) || isfolder(chemin))
            uialert(fig, ['Not found: ' chemin], 'Open');
            return
        end
        if ispc
            winopen(chemin);
        elseif ismac
            system(['open "' chemin '"']);
        else
            system(['xdg-open "' chemin '"']);
        end
    end
end


function t = texte_K(K, long)
% 'consensus' for K = 0 (step 2, rule 'consensus'), otherwise the K
if nargin < 2, long = false; end
if K == 0
    t = 'consensus';
elseif long
    t = sprintf('K = %g', K);
else
    t = sprintf('%g', K);
end
end


function v = valeur_3b(c)
% étape 3b demandée dans la config (false si le champ n'existe pas)
v = isfield(c, 'faire_etape3b') && c.faire_etape3b;
end


function v = valeur_5(c)
% étape 5 demandée dans la config (false si le champ n'existe pas)
v = isfield(c, 'faire_etape5') && c.faire_etape5;
end


function t = onoff(vrai)
if vrai, t = 'on'; else, t = 'off'; end
end


function h = placer(h, ligne, colonne)
% places a component in a grid
h.Layout.Row = ligne;
h.Layout.Column = colonne;
end
