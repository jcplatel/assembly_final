function fig = interface_couleurs(fichier_sessions)
%% =========================================================================
%  INTERFACE_COULEURS : attribution et vérification des couleurs Brainbow
%  (color_maxwellV2, 17 nov., demandée par l'utilisateur)
%  -------------------------------------------------------------------------
%  À gauche : toutes les sessions de la table des sessions (sessions_brainbow,
%  seulement les lignes à inclure = 1 depuis le 17 nov., toutes conditions), avec les fichiers
%  de couleur présents dans leur dossier suite2p.
%  Au centre : l'image Brainbow avec les masques des cellules (Fall.mat,
%  iscell = 1, même ordre que colorcell et que le pipeline), colorés selon un
%  fichier de couleur au choix ; comparaison avec un 2e fichier (accord,
%  confusions, cellules en désaccord entourées). Souris sur l'image : molette
%  = zoom autour du pointeur, clic droit maintenu = déplacer, double-clic =
%  image entière ; le zoom est gardé quand l'image est redessinée. Clic sur une cellule : ses
%  couleurs dans chaque fichier, ses valeurs R G B, et boutons pour changer
%  sa couleur dans la copie de travail.
%  À droite : les étapes du traitement de color_maxwellV2, une par bouton,
%  dans l'ordre de color_predictor_model (fonctions de l'utilisateur, non
%  modifiées) ; puis enregistrement dans colorcell_MaxwellAuto.mat.
%  17 nov. : méthode MaxwellAuto (comme colorcell_maxwellauto) : les étapes 4
%  et 6 proposent au départ max(R,V,B) < 0,20 puis < 0,15 (modifiable à la
%  main) ; l'étape 7 utilise l'échelle 'ancien', pas de classe blanche, et
%  met en noir les couleurs proches d'une frontière (marge < champ « margin »).
%
%  Lancer :  interface_couleurs                    (table de config_brainbow si
%                                                   le pipeline est dans le chemin,
%                                                   sinon elle est demandée)
%            interface_couleurs('E:\...\sessions_brainbow.xlsx')
%
%  Le dossier color_maxwellV2 doit être en tête du chemin MATLAB (addpath
%  ..., '-begin'), avec les autres fonctions de l'utilisateur
%  (brainbow_image_read, interactive_threshold_selection, triangle_maxwell,
%  interactive_color_correction...).
%
%  Fichiers écrits (dossier suite2p de la session) :
%    colorcell_MaxwellAuto.mat (variable colorcell ; l'ancien, s'il existe,
%                               est d'abord renommé colorcell_MaxwellAuto_<date>.mat)
%    couleurs_details_Auto.mat (tout ce qui a servi : valeurs, proportions,
%                               seuils, fuite, recalage, corrections à la main)
%    id_low_intensity.mat      (comme color_predictor_model, aux étapes 4 et 6)
%  Aucun autre fichier de couleur n'est modifié.
%  =========================================================================

if nargin < 1, fichier_sessions = ''; end

NOMS = {'R', 'G', 'B', 'RG', 'RB', 'GB', 'white', 'none'};
RGB_C = [1 0.15 0.15; 0.15 1 0.15; 0.3 0.55 1; 1 1 0.1; 1 0.2 1; 0.1 1 1; 1 1 1; 0.55 0.55 0.55];
TRAVAIL = 'working copy (not saved)';
IMG_CORR = 'corrected files (aligned_*_corrected.tif)';
IMG_AUTO = 'corrected as MaxwellAuto (GCaMP leak and background removed)';

% ---------------- état ----------------
sessions = table();
S = struct();               % session courante (géométrie, fichiers, traitement)
selection = 0;              % cellule sélectionnée
vue = struct('cle', '', 'xlim', [], 'ylim', []);   % zoom gardé entre deux dessins
glisse = [];                % déplacement au clic droit : point de départ
historique = {};            % pour « Undo »
nouvelle_session = false;   % vrai pendant l'ouverture d'une session : choisir le fichier par défaut

% ---------------- fenêtre ----------------
fig = uifigure('Name', 'Brainbow colours (color_maxwellV2)', 'Position', [40 40 1560 940]);
G = uigridlayout(fig, [1 3]);
G.ColumnWidth = {390, '1x', 300};      % 17 nov. : sessions 390 (300 avant), droite 300 (340 avant)

% --- gauche : sessions
GL = uigridlayout(G, [5 2]);
GL.RowHeight = {24, 24, 24, '1x', 60};
GL.ColumnWidth = {'1x', 80};
uilabel(GL, 'Text', 'Sessions', 'FontWeight', 'bold');
b_table = uibutton(GL, 'Text', 'Table...', 'ButtonPushedFcn', @(~, ~) choisir_table());
lbl_table = uilabel(GL, 'Text', '', 'FontSize', 9);
lbl_table.Layout.Row = 2; lbl_table.Layout.Column = [1 2];
ef_filtre = uieditfield(GL, 'text', 'Placeholder', 'filter (animal, date, plane...)', ...
    'ValueChangedFcn', @(~, ~) remplir_liste());
ef_filtre.Layout.Row = 3; ef_filtre.Layout.Column = [1 2];
lb_sessions = uilistbox(GL, 'Items', {}, 'FontName', 'Consolas', 'FontSize', 10, ...
    'ValueChangedFcn', @(~, ~) ouvrir_session());
lb_sessions.Layout.Row = 4; lb_sessions.Layout.Column = [1 2];
lbl_leg = uilabel(GL, 'Text', ['Before |: file read by the pipeline ((t) = forced in the table). ' ...
    'After |: files present. GT, new = Maxwellnew, old = Maxwell, V2 = MaxwellV2, Auto = MaxwellAuto.'], ...
    'WordWrap', 'on', 'FontSize', 9);
lbl_leg.Layout.Row = 5; lbl_leg.Layout.Column = [1 2];

% --- centre : image
GC = uigridlayout(G, [3 1]);
GC.RowHeight = {56, '1x', 150};
GO = uigridlayout(GC, [2 10]);
GO.ColumnWidth = {90, 170, 95, 170, 70, 70, 70, 70, 60, '1x'};
GO.RowHeight = {22, 22};
GO.Padding = [0 0 0 0];
uilabel(GO, 'Text', 'Colours of:');
dd_affiche = uidropdown(GO, 'Items', {'(outlines only)'}, 'ValueChangedFcn', @(~, ~) changer_affichage());
uilabel(GO, 'Text', 'Compare with:');
dd_compare = uidropdown(GO, 'Items', {'(none)'}, 'ValueChangedFcn', @(~, ~) changer_affichage());
cb_R = uicheckbox(GO, 'Text', 'Red', 'Value', true, 'ValueChangedFcn', @(~, ~) dessiner());
cb_V = uicheckbox(GO, 'Text', 'Green', 'Value', true, 'ValueChangedFcn', @(~, ~) dessiner());
cb_B = uicheckbox(GO, 'Text', 'Blue', 'Value', true, 'ValueChangedFcn', @(~, ~) dessiner());
cb_Ca = uicheckbox(GO, 'Text', 'GCaMP', 'Value', false, 'ValueChangedFcn', @(~, ~) dessiner());
uilabel(GO, 'Text', 'Gain');
sp_gain = uispinner(GO, 'Value', 1, 'Limits', [0.2 8], 'Step', 0.2, 'ValueChangedFcn', @(~, ~) dessiner());
uilabel(GO, 'Text', 'Image:');
dd_image = uidropdown(GO, 'Items', {'composite file'}, 'ValueChangedFcn', @(~, ~) dessiner());
uilabel(GO, 'Text', 'Show:');
% 17 nov. : une seule couleur à la fois (comme l'interface de color_predictor_model)
dd_masques = uidropdown(GO, 'Items', {'all cells', 'disagreements only', 'no masks', 'only R', 'only G', ...
    'only B', 'only RG (yellow)', 'only RB (magenta)', 'only GB (cyan)', 'only white', 'only none'}, ...
    'ValueChangedFcn', @(~, ~) dessiner());
cb_num = uicheckbox(GO, 'Text', 'Numbers', 'Value', false, 'ValueChangedFcn', @(~, ~) dessiner());
cb_num.Layout.Row = 2; cb_num.Layout.Column = 5;
b_img = uibutton(GO, 'Text', 'Image file...', 'ButtonPushedFcn', @(~, ~) choisir_image());
b_img.Layout.Row = 2; b_img.Layout.Column = [6 7];
b_fig = uibutton(GO, 'Text', 'Triangle','ButtonPushedFcn', @(~, ~) montrer_figures(), ...
    'Tooltip', 'Maxwell triangle and colour proportions vs reference, for the file in "Colours of"');
b_fig.Layout.Row = 2; b_fig.Layout.Column = [8 10];
% 17 nov. : image à gauche, petit triangle de Maxwell à droite (cellule cliquée mise en évidence)
GI = uigridlayout(GC, [1 2]);
GI.ColumnWidth = {'1x', 230}; GI.Padding = [0 0 0 0]; GI.ColumnSpacing = 4;
ax = uiaxes(GI);
GTri = uigridlayout(GI, [2 1]);
GTri.RowHeight = {215, '1x'}; GTri.Padding = [0 0 0 0];
axT = uiaxes(GTri);
axT.XTick = []; axT.YTick = []; axT.XColor = 'none'; axT.YColor = 'none';
disableDefaultInteractivity(axT);
lbl_tri = uilabel(GTri, 'Text', '', 'WordWrap', 'on', 'FontSize', 9, 'VerticalAlignment', 'top');
tri = struct('cle', '', 'D', [], 'source', '');     % coordonnées du triangle gardées en mémoire
ax.XTick = []; ax.YTick = [];
disableDefaultInteractivity(ax);
fig.WindowScrollWheelFcn = @(~, e) molette(e);
fig.WindowButtonMotionFcn = @(~, ~) rien();   % pour que ax.CurrentPoint suive la souris (molette)
ta_comp = uitextarea(GC, 'Editable', 'off', 'FontName', 'Consolas', 'FontSize', 10);

% --- droite : étapes et cellule
GR = uigridlayout(G, [24 4]);
GR.RowHeight = [{22}, repmat({24}, 1, 11), {22, 22, 24, 24, 24, 24, 22, 22, 24, 24, '1x', 22}];
GR.ColumnWidth = {28, '1x', 60, 60};
t1 = uilabel(GR, 'Text', 'Processing (color_maxwellV2)', 'FontWeight', 'bold');
t1.Layout.Column = [1 4];
noms_etapes = {'1. Read and normalise images', '2. Align (multimodal)', '3. Extract ROI values', ...
    '4. Low-intensity threshold', '5. Bleed-through correction', '6. Threshold after correction', ...
    '7. Classify (Maxwell)', '8. Maxwell triangle', '9. Manual correction (your tool)', ...
    '10. Save colorcell_MaxwellAuto'};
lbl_ok = gobjects(1, 10); b_etape = gobjects(1, 10);
for k = 1:10
    lbl_ok(k) = uilabel(GR, 'Text', '', 'FontColor', [0 0.55 0], 'FontWeight', 'bold');
    lbl_ok(k).Layout.Row = 1 + k; lbl_ok(k).Layout.Column = 1;
    b_etape(k) = uibutton(GR, 'Text', noms_etapes{k}, 'HorizontalAlignment', 'left', ...
        'ButtonPushedFcn', @(~, ~) lancer_etape(k));
    b_etape(k).Layout.Row = 1 + k; b_etape(k).Layout.Column = [2 4];
end
b_123 = uibutton(GR, 'Text', 'Run 1 to 3', 'ButtonPushedFcn', @(~, ~) lancer_1_3());
b_123.Layout.Row = 12; b_123.Layout.Column = [2 3];
l_ep = uilabel(GR, 'Text', 'margin (step 7)');
l_ep.Layout.Row = 13; l_ep.Layout.Column = [1 2];
ef_ep = uieditfield(GR, 'numeric', 'Value', 0, 'Limits', [0 0.5], 'Visible', 'off', ...
    'Tooltip', 'ep: width of the uncertain band (not used since 17 Nov.)');
ef_ep.Layout.Row = 13; ef_ep.Layout.Column = 3;
% marge (17 nov., méthode MaxwellAuto) : noire si la couleur est à moins de cette distance d'une frontière
ef_marge = uieditfield(GR, 'numeric', 'Value', 0.05, 'Limits', [0 0.5], ...
    'Tooltip', 'margin: cells closer than this to a border of the triangle become black (0 = keep all)');
ef_marge.Layout.Row = 13; ef_marge.Layout.Column = 4;
l_seuil = uilabel(GR, 'Text', 'steps 4, 6 start: max(R,G,B) < 0.20 / 0.15', 'FontSize', 9, ...
    'Tooltip', 'cells proposed as black at the start of steps 4 and 6; the slider starts at 0 and can be moved');
l_seuil.Layout.Row = 14; l_seuil.Layout.Column = [1 2];
ef_seuil = uieditfield(GR, 'numeric', 'Value', 0.25, 'Limits', [0 1], ...
    'Tooltip', 'slider threshold used only if no cell is proposed at the start');
ef_seuil.Layout.Row = 14; ef_seuil.Layout.Column = 3;
b_edit = uibutton(GR, 'Text', 'Edit a copy of the displayed file', 'ButtonPushedFcn', @(~, ~) copie_travail());
b_edit.Layout.Row = 15; b_edit.Layout.Column = [1 4];

t2 = uilabel(GR, 'Text', 'Selected cell (click on the image)', 'FontWeight', 'bold');
t2.Layout.Row = 16; t2.Layout.Column = [1 4];
lbl_cellule = uilabel(GR, 'Text', 'none', 'WordWrap', 'on', 'FontSize', 10, 'VerticalAlignment', 'top');
lbl_cellule.Layout.Row = [17 18]; lbl_cellule.Layout.Column = [1 4];
l_set = uilabel(GR, 'Text', 'Set colour in the working copy:');
l_set.Layout.Row = 19; l_set.Layout.Column = [1 4];
GB = uigridlayout(GR, [2 4]);
GB.Layout.Row = [20 21]; GB.Layout.Column = [1 4];
GB.Padding = [0 0 0 0]; GB.RowSpacing = 2; GB.ColumnSpacing = 2;
for c = 1:8
    uibutton(GB, 'Text', NOMS{c}, 'BackgroundColor', RGB_C(c, :), ...
        'FontColor', ifelse(c == 3 || c == 8, [1 1 1], [0 0 0]), 'ButtonPushedFcn', @(~, ~) mettre_couleur(c));
end
b_undo = uibutton(GR, 'Text', 'Undo', 'ButtonPushedFcn', @(~, ~) annuler());
b_undo.Layout.Row = 22; b_undo.Layout.Column = [3 4];
lbl_etat = uilabel(GR, 'Text', '', 'WordWrap', 'on', 'FontSize', 9, 'VerticalAlignment', 'top');
lbl_etat.Layout.Row = [23 24]; lbl_etat.Layout.Column = [1 4];

% ---------------- démarrage ----------------
etat_boutons();
if isempty(fichier_sessions)
    try
        cfg = config_brainbow();
        fichier_sessions = cfg.fichier_sessions;
    catch
    end
end
if isempty(fichier_sessions) || ~isfile(fichier_sessions)
    choisir_table();
else
    charger_table(fichier_sessions);
end


%% ======================= sessions =======================
    function choisir_table()
        [f, p] = uigetfile({'*.xlsx;*.xls;*.csv', 'Session table'}, 'Session table (sessions_brainbow)');
        figure(fig);
        if isequal(f, 0), return; end
        charger_table(fullfile(p, f));
    end

    function charger_table(f)
        try
            sessions = lire_sessions(f);
        catch err
            uialert(fig, err.message, 'Session table');
            return
        end
        lbl_table.Text = f;
        remplir_liste();
    end

    function remplir_liste()
        if isempty(sessions), lb_sessions.Items = {}; return; end
        filtre = strtrim(ef_filtre.Value);
        items = {}; donnees = [];
        for s = 1:height(sessions)
            if ~isempty(filtre) && isempty(strfind(lower(sessions.nom{s}), lower(filtre))), continue; end %#ok<STREMP>
            items{end+1} = sprintf('%-23s%s', regexprep(sessions.nom{s}, '_sans_indice$', ''), resume_fichiers(s)); %#ok<AGROW>
            donnees(end+1) = s; %#ok<AGROW>
        end
        lb_sessions.Items = items;
        lb_sessions.ItemsData = donnees;
    end

    function t = resume_fichiers(s)
        % 17 nov. : d'abord le fichier lu par le pipeline (même ordre que charger_couleurs :
        % colonne fichier_colorcell, puis GT, Maxwellnew, Maxwell), puis les fichiers présents
        d = sessions.dossier{s};
        court = @(nom) regexprep(regexprep(nom, '^colorcell_', ''), ...
            {'\.mat$', '^Maxwellnew$', '^Maxwell$', '^MaxwellV2$', '^MaxwellAuto$'}, {'', 'new', 'old', 'V2', 'Auto'});
        utilise = '';
        if ~isempty(sessions.fichier_colorcell{s})
            [~, nf] = fileparts(sessions.fichier_colorcell{s});
            utilise = [court(nf) ' (t)'];
        else
            for x = {'GT', 'Maxwellnew', 'Maxwell'}
                if isfile(fullfile(d, ['colorcell_' x{1} '.mat'])), utilise = court(x{1}); break; end
            end
        end
        if isempty(utilise), utilise = '-'; end
        t = '';
        for x = {'GT', 'Maxwellnew', 'Maxwell', 'MaxwellV2', 'MaxwellAuto'}
            if isfile(fullfile(d, ['colorcell_' x{1} '.mat'])), t = [t ' ' court(x{1})]; end %#ok<AGROW>
        end
        t = sprintf('%-9s|%s', utilise, t);
    end

    function ouvrir_session()
        s = lb_sessions.Value;
        if isempty(s) || iscell(s), return; end
        if isfield(S, 'modifie') && S.modifie
            r = uiconfirm(fig, 'The working copy of the current session is not saved. Discard it?', ...
                'Unsaved changes', 'Options', {'Discard', 'Cancel'}, 'DefaultOption', 2);
            if strcmp(r, 'Cancel'), return; end
        end
        d = uiprogressdlg(fig, 'Title', 'Loading', 'Indeterminate', 'on');
        try
            S = charger_session(sessions(s, :));
        catch err
            close(d);
            uialert(fig, err.message, sessions.nom{s});
            return
        end
        close(d);
        selection = 0;
        historique = {};
        dd_affiche.Items = {'(outlines only)'};     % nouvelle session : fichier par défaut (17 nov.)
        dd_affiche.Value = '(outlines only)';
        nouvelle_session = true;
        maj_listes();
        nouvelle_session = false;
        changer_affichage();
        etat_boutons();
    end

    function T = charger_session(ligne)
        T = struct();
        T.nom = ligne.nom{1};
        T.dossier = ligne.dossier{1};
        T.force = ligne.fichier_colorcell{1};
        F = load(fullfile(T.dossier, 'Fall.mat'), 'ops');
        % 18 nov. : masques corrigés (masques_corriges.mat) s'ils existent
        [F.stat, F.iscell, T.masques_corriges] = charger_masques(T.dossier);
        roi = find(F.iscell(:, 1) == 1);
        n = numel(roi);
        T.n = n;
        T.Ly = F.ops.Ly; T.Lx = F.ops.Lx;
        T.bx = cell(1, n); T.by = cell(1, n); T.centre = zeros(n, 2);
        T.cx = cell(1, n); T.cy = cell(1, n);                  % contours continus (17 nov.)
        T.label = zeros(T.Ly, T.Lx, 'uint16');
        for m = 1:n
            if iscell(F.stat), st = F.stat{roi(m)}; else, st = F.stat(roi(m)); end
            if T.masques_corriges
                % contour recalculé à partir des pixels du masque corrigé
                xq = double(st.xpix(:)) + 1; yq = double(st.ypix(:)) + 1;
                bwq = false(max(yq) - min(yq) + 3, max(xq) - min(xq) + 3);
                bwq(sub2ind(size(bwq), yq - min(yq) + 2, xq - min(xq) + 2)) = true;
                Bq = bwboundaries(bwq, 8, 'noholes'); [~, gq] = max(cellfun(@numel, Bq));
                T.bx{m} = Bq{gq}(:, 2) + min(xq) - 2; T.by{m} = Bq{gq}(:, 1) + min(yq) - 2;
            else
                T.bx{m} = double(st.xext(:)) + 1;
                T.by{m} = double(st.yext(:)) + 1;
            end
            [T.cx{m}, T.cy{m}] = contour_cellule(T.bx{m}, T.by{m});
            xp = double(st.xpix(:)) + 1; yp = double(st.ypix(:)) + 1;
            T.centre(m, :) = [mean(xp), mean(yp)];
            ok = xp >= 1 & xp <= T.Lx & yp >= 1 & yp <= T.Ly;
            T.label(sub2ind([T.Ly T.Lx], yp(ok), xp(ok))) = m;
        end
        % fichiers de couleur présents
        L = dir(fullfile(T.dossier, 'colorcell*.mat'));
        T.fichiers = struct('nom', {}, 'colorcell', {});
        for k = 1:numel(L)
            try
                C = load(fullfile(L(k).folder, L(k).name), 'colorcell');
                T.fichiers(end+1) = struct('nom', L(k).name, 'colorcell', double(C.colorcell(:))); %#ok<AGROW>
            catch
            end
        end
        T.travail = [];
        T.modifie = false;
        T.etapes = false(1, 10);
        T.P = struct();          % résultats du traitement
        T.canaux = {};
        T.fichier_image = trouver_composite_local(T.dossier);
        if ~isempty(T.fichier_image), T.canaux = lire_composite(T.fichier_image); end
        % images corrigées de la fuite de color_predictor_model (aligned_*_corrected.tif, 8 bits),
        % copiées du dossier d'analyse dans le dossier de session (17 nov.), sinon à côté du composite
        % 17 nov. : images corrigées comme MaxwellAuto (k et b de la fuite et fond de
        % couleurs_details_Auto.mat appliqués à compositecalcium_V2.tif)
        T.canaux_auto = {};
        fda = fullfile(T.dossier, 'couleurs_details_Auto.mat');
        if isfile(fda) && endsWith(T.fichier_image, 'compositecalcium_V2.tif') && numel(T.canaux) >= 4
            try
                Z = load(fda, 'params_fuite', 'fond_image');
                k = Z.params_fuite.k; b = Z.params_fuite.b;          % ordre vert rouge bleu
                fo = [0 0 0]; if isfield(Z, 'fond_image'), fo = Z.fond_image; end   % ordre R V B
                ca = T.canaux{4} / 65535;
                T.canaux_auto = {max(0, T.canaux{1} / 65535 - (k(2) * ca + b(2)) - fo(1)), ...
                    max(0, T.canaux{2} / 65535 - (k(1) * ca + b(1)) - fo(2)), ...
                    max(0, T.canaux{3} / 65535 - (k(3) * ca + b(3)) - fo(3)), T.canaux{4}};
            catch
                T.canaux_auto = {};
            end
        end
        T.canaux_corr = {};
        dossiers = {T.dossier};
        if ~isempty(T.fichier_image), dossiers{end+1} = fileparts(T.fichier_image); end
        for dd = dossiers
            fc = fullfile(dd{1}, {'aligned_red_corrected.tif', 'aligned_green_corrected.tif', ...
                'aligned_blue_corrected.tif'});
            if ~all(cellfun(@isfile, fc)), continue; end
            try
                chc = cellfun(@(f) double(imread(f)), fc, 'UniformOutput', false);
                if isequal(size(chc{1}), [T.Ly T.Lx]), T.canaux_corr = chc; break; end
            catch
            end
        end
    end

    function maj_listes()
        items = {S.fichiers.nom};
        tailles = arrayfun(@(f) numel(f.colorcell), S.fichiers);
        for k = 1:numel(items)
            if tailles(k) ~= S.n
                items{k} = sprintf('%s  [SIZE %d, %d cells]', items{k}, tailles(k), S.n);
            end
        end
        if ~isempty(S.travail), items{end+1} = TRAVAIL; end
        ancien = dd_affiche.Value;
        dd_affiche.Items = [{'(outlines only)'}, items];
        dd_compare.Items = [{'(none)'}, items];
        if ~nouvelle_session && any(strcmp(dd_affiche.Items, ancien))
            dd_affiche.Value = ancien;
        else
            % par défaut : le fichier que le pipeline lit
            dd_affiche.Value = choix_defaut(items);
        end
        if ~any(strcmp(dd_compare.Items, dd_compare.Value)), dd_compare.Value = '(none)'; end
        imgs = {'composite file'};
        if ~isempty(S.canaux_auto), imgs{end+1} = IMG_AUTO; end
        if ~isempty(S.canaux_corr), imgs{end+1} = IMG_CORR; end
        if isfield(S.P, 'aligned_red'), imgs{end+1} = 'aligned (step 2)'; end
        if isfield(S.P, 'aligned_red_corrected'), imgs{end+1} = 'corrected (step 5)'; end
        dd_image.Items = imgs;
    end

    function v = choix_defaut(items)
        v = '(outlines only)';
        if ~isempty(S.force)
            [~, nf, ext] = fileparts(S.force);
            k = find(startsWith(items, [nf ext]), 1);
            if ~isempty(k), v = items{k}; return; end
        end
        for x = {'colorcell_GT.mat', 'colorcell_Maxwellnew.mat', 'colorcell_Maxwell.mat'}
            k = find(startsWith(items, x{1}), 1);
            if ~isempty(k), v = items{k}; return; end
        end
        if ~isempty(items), v = items{1}; end
    end


%% ======================= affichage =======================
    function c = couleurs_de(nom)
        c = [];
        if strcmp(nom, TRAVAIL)
            c = S.travail;
        else
            nom = regexprep(nom, '\s+\[SIZE.*$', '');
            k = find(strcmp({S.fichiers.nom}, nom), 1);
            if ~isempty(k), c = S.fichiers(k).colorcell; end
        end
        if numel(c) ~= S.n, c = []; end
    end

    function montrer_figures()
        % triangle de Maxwell + proportions vs référence, pour le fichier affiché (17 nov.)
        if ~isfield(S, 'n'), uialert(fig, 'Choose a session first.', 'Triangle + proportions'); return; end
        c = couleurs_de(dd_affiche.Value);
        if isempty(c)
            uialert(fig, 'Choose a colour file (with the right number of cells) in "Colours of".', ...
                'Triangle + proportions');
            return
        end
        [D, source] = coordonnees_triangle();
        figure_couleurs(c, D, sprintf('%s - %s - triangle: %s', S.nom, dd_affiche.Value, source));
    end

    function [D, source] = coordonnees_triangle()
        % coordonnées du triangle, dans l'ordre :
        %   1. étape 7 faite dans cette session ;
        %   2. fichier de détails du dossier (couleurs_details_V2 / couleurs_details) :
        %      le triangle enregistré, sinon les valeurs corrigées des ROI ;
        %   3. sinon : moyenne de chaque ROI sur l'image affichée (menu « Image: »).
        D = []; source = 'none';
        if isfield(S.P, 'Donnees_Clustering') && size(S.P.Donnees_Clustering, 1) == S.n
            D = S.P.Donnees_Clustering; source = 'step 7';
            if isfield(S.P, 'r')                % cellules rejetées (r = g = b = 0) : hors du triangle
                D(S.P.r(:) + S.P.g(:) + S.P.b(:) == 0, :) = NaN;
            end
            return
        end
        details_fichiers = {'couleurs_details_V2.mat', 'couleurs_details.mat', 'couleurs_details_Auto.mat'};
        if contains(dd_affiche.Value, 'MaxwellAuto'), details_fichiers = details_fichiers([3 1 2]); end
        for fd = details_fichiers
            f = fullfile(S.dossier, fd{1});
            if ~isfile(f), continue; end
            presentes = who('-file', f);        % ne lire que ce qui sert (le fichier peut être gros)
            for champ = {'Donnees_Clustering', 'triangle'}
                if ~ismember(champ{1}, presentes), continue; end
                X = load(f, champ{1});
                if size(X.(champ{1}), 1) == S.n
                    D = X.(champ{1}); source = fd{1};
                    if all(ismember({'r', 'g', 'b'}, presentes))   % cellules rejetées : hors du triangle
                        Y = load(f, 'r', 'g', 'b');
                        if numel(Y.r) == S.n, D(Y.r(:) + Y.g(:) + Y.b(:) == 0, :) = NaN; end
                    end
                    return
                end
            end
            v3 = {'roi_red_corrected3', 'roi_green_corrected3', 'roi_blue_corrected3'};
            if all(ismember(v3, presentes))
                X = load(f, v3{:});
                if numel(X.roi_red_corrected3) == S.n
                    D = maxwell_proportions(X.roi_red_corrected3(:), X.roi_green_corrected3(:), X.roi_blue_corrected3(:));
                    source = [fd{1} ' (corrected ROI values)'];
                    return
                end
            end
            if ismember('roi_corrige', presentes)          % [rouge vert bleu] (color_predictor_model)
                X = load(f, 'roi_corrige');
                if size(X.roi_corrige, 1) == S.n
                    D = maxwell_proportions(X.roi_corrige(:, 1), X.roi_corrige(:, 2), X.roi_corrige(:, 3));
                    source = [fd{1} ' (corrected ROI values)'];
                    return
                end
            end
        end
        % 3. à partir des images : moyenne des pixels de chaque ROI (masques Fall.mat)
        nom_image = dd_image.Value;
        switch dd_image.Value
            case 'aligned (step 2)'
                ch = {S.P.aligned_red, S.P.aligned_green, S.P.aligned_blue};
            case 'corrected (step 5)'
                ch = {S.P.aligned_red_corrected, S.P.aligned_green_corrected, S.P.aligned_blue_corrected};
            otherwise
                % images corrigées de la fuite de color_predictor_model en priorité
                if ~isempty(S.canaux_corr)
                    ch = S.canaux_corr; nom_image = 'aligned_*_corrected.tif';
                else
                    ch = S.canaux;
                end
        end
        if numel(ch) < 3 || ~isequal(size(ch{1}), [S.Ly S.Lx]), return; end
        dedans = S.label > 0;
        lab = double(S.label(dedans));
        nb = accumarray(lab, 1, [S.n 1]);
        moy = zeros(S.n, 3);
        for k = 1:3
            v = double(ch{k});
            moy(:, k) = accumarray(lab, v(dedans), [S.n 1]) ./ max(nb, 1);
        end
        D = maxwell_proportions(moy(:, 1), moy(:, 2), moy(:, 3));
        source = sprintf('ROI means on %s', nom_image);
    end

    function changer_affichage()
        dessiner();
        comparer();
        afficher_cellule();
        dessiner_triangle();
    end

    function dessiner()
        if ~isempty(vue.cle) && ~isempty(ax.Children)
            vue.xlim = ax.XLim; vue.ylim = ax.YLim;     % zoom actuel, gardé si même image
        end
        cla(ax);
        if ~isfield(S, 'n'), return; end
        hold(ax, 'on');
        rgb = image_rgb();
        im = image(ax, rgb);
        im.ButtonDownFcn = @(~, e) clic(e);
        c1 = couleurs_de(dd_affiche.Value);
        c2 = couleurs_de(dd_compare.Value);
        montrer = true(S.n, 1);
        if strcmp(dd_masques.Value, 'disagreements only') && ~isempty(c1) && ~isempty(c2)
            montrer = c1 ~= c2;
        end
        if startsWith(dd_masques.Value, 'only ') && ~isempty(c1)
            % une seule couleur : rang du menu après 'no masks' = numéro de la classe (1 à 8)
            classe_vue = find(strcmp(dd_masques.Items, dd_masques.Value)) - 3;
            montrer = c1(:) == classe_vue;
        end
        if ~strcmp(dd_masques.Value, 'no masks')
            if isempty(c1)
                [X, Y] = joindre(S.cx, S.cy, find(montrer));
                plot(ax, X, Y, '-', 'Color', [1 1 1], 'LineWidth', 0.5, 'HitTest', 'off');
            else
                for c = 1:8
                    m = find(c1 == c & montrer);
                    if isempty(m), continue; end
                    [X, Y] = joindre(S.cx, S.cy, m);
                    plot(ax, X, Y, '-', 'Color', RGB_C(c, :), 'LineWidth', 1, 'HitTest', 'off');
                end
            end
            if ~isempty(c1) && ~isempty(c2)
                m = find(c1 ~= c2 & montrer);
                plot(ax, S.centre(m, 1), S.centre(m, 2), 'o', 'Color', [1 1 1], 'MarkerSize', 9, ...
                    'LineWidth', 1, 'HitTest', 'off');
            end
        end
        if cb_num.Value
            text(ax, S.centre(:, 1), S.centre(:, 2), arrayfun(@num2str, 1:S.n, 'UniformOutput', false), ...
                'Color', [1 1 0.6], 'FontSize', 7, 'HorizontalAlignment', 'center', 'Clipping', 'on', 'HitTest', 'off');
        end
        if selection >= 1 && selection <= S.n
            plot(ax, S.cx{selection}, S.cy{selection}, '-', 'Color', [1 1 0], 'LineWidth', 2.5, 'HitTest', 'off');
        end
        hold(ax, 'off');
        ax.DataAspectRatio = [1 1 1];
        ax.YDir = 'reverse';
        cle = sprintf('%s|%d|%d', S.nom, size(rgb, 1), size(rgb, 2));
        if strcmp(cle, vue.cle) && ~isempty(vue.xlim)
            ax.XLim = vue.xlim; ax.YLim = vue.ylim;
        else
            ax.XLim = [0.5, size(rgb, 2) + 0.5]; ax.YLim = [0.5, size(rgb, 1) + 0.5];
        end
        vue.cle = cle;
        ax.XTick = []; ax.YTick = [];
        t = sprintf('%s (%d cells)', S.nom, S.n);
        if strcmp(dd_image.Value, 'composite file')
            if isempty(S.fichier_image)
                t = [t ' - no composite image found: run step 1-2 or "Image file..."'];
            elseif size(rgb, 1) ~= S.Ly || size(rgb, 2) ~= S.Lx
                t = [t sprintf(' - WARNING: image %d x %d, masks %d x %d', size(rgb, 1), size(rgb, 2), S.Ly, S.Lx)];
            end
        end
        title(ax, t, 'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 10);
    end

    function rgb = image_rgb()
        % canaux : rouge, vert, bleu, calcium
        switch dd_image.Value
            case 'aligned (step 2)'
                ch = {S.P.aligned_red, S.P.aligned_green, S.P.aligned_blue, S.P.calcium_norm};
            case 'corrected (step 5)'
                ch = {S.P.aligned_red_corrected, S.P.aligned_green_corrected, S.P.aligned_blue_corrected, ...
                    S.P.calcium_norm};
            case IMG_AUTO
                ch = S.canaux_auto;
            case IMG_CORR
                ch = S.canaux_corr;
                if numel(S.canaux) >= 4 && isequal(size(S.canaux{4}), size(ch{1})), ch{4} = S.canaux{4}; end
            otherwise
                ch = S.canaux;
        end
        if isempty(ch)
            rgb = zeros(S.Ly, S.Lx, 3);
            return
        end
        [Ly, Lx] = size(ch{1});
        rgb = zeros(Ly, Lx, 3);
        actifs = [cb_R.Value, cb_V.Value, cb_B.Value, cb_Ca.Value];
        cols = [1 0 0; 0 1 0; 0 0 1; 1 1 1];
        for k = 1:min(4, numel(ch))
            if ~actifs(k), continue; end
            v = double(ch{k});
            bornes = [prctile(v(:), 1), prctile(v(:), 99.7)];
            v = (v - bornes(1)) / max(bornes(2) - bornes(1), eps) * sp_gain.Value;
            v = min(max(v, 0), 1);
            rgb = rgb + v .* reshape(cols(k, :), 1, 1, 3);
        end
        rgb = min(rgb, 1);
    end

    function choisir_image()
        if ~isfield(S, 'dossier'), return; end
        [f, p] = uigetfile({'*.tif;*.tiff', 'TIFF'}, 'Composite image (R, G, B, GCaMP)', S.dossier);
        figure(fig);
        if isequal(f, 0), return; end
        S.canaux = lire_composite(fullfile(p, f));
        S.fichier_image = fullfile(p, f);
        dd_image.Value = 'composite file';
        dessiner();
    end


%% ======================= comparaison =======================
    function comparer()
        L = {};
        if ~isfield(S, 'n'), ta_comp.Value = {''}; return; end
        c1 = couleurs_de(dd_affiche.Value);
        c2 = couleurs_de(dd_compare.Value);
        if ~isempty(c1)
            L{end+1} = sprintf('%s : %s', dd_affiche.Value, composition(c1));
            L{end+1} = sprintf('   same colour among the 3 nearest coloured neighbours / chance : %.2f', ...
                indice_voisins(c1));
        end
        if ~isempty(c2)
            L{end+1} = sprintf('%s : %s', dd_compare.Value, composition(c2));
            L{end+1} = sprintf('   same colour among the 3 nearest coloured neighbours / chance : %.2f', ...
                indice_voisins(c2));
        end
        if ~isempty(c1) && ~isempty(c2)
            col1 = c1 >= 1 & c1 <= 6; col2 = c2 >= 1 & c2 <= 6; deux = col1 & col2;
            p1 = histc(c1(deux), 1:6) / max(sum(deux), 1); p2 = histc(c2(deux), 1:6) / max(sum(deux), 1); %#ok<HISTC>
            L{end+1} = sprintf(['Agreement on cells coloured in both (%d): %.1f %% (chance %.1f %%) ; ' ...
                'different colour or class: %d cells (circled)'], sum(deux), 100 * mean(c1(deux) == c2(deux)), ...
                100 * sum(p1(:) .* p2(:)), sum(c1 ~= c2));
            L{end+1} = 'Rows: shown file ; columns: compared file';
            codes = unique([c1; c2])';
            L{end+1} = ['        ' strjoin(arrayfun(@(b) sprintf('%7s', nom_c(b)), codes, 'UniformOutput', false), '')];
            for a = codes
                L{end+1} = [sprintf('%7s ', nom_c(a)) sprintf('%7d', arrayfun(@(b) sum(c1 == a & c2 == b), codes))]; %#ok<AGROW>
            end
        end
        if isempty(L), L = {'Choose a colour file to show (and one to compare with).'}; end
        ta_comp.Value = L;
    end

    function t = composition(c)
        n_col = sum(c >= 1 & c <= 6);
        parts = arrayfun(@(k) sprintf('%s %d', NOMS{k}, sum(c == k)), 1:8, 'UniformOutput', false);
        t = sprintf('%d coloured / %d ; %s', n_col, numel(c), strjoin(parts, ', '));
    end

    function r = indice_voisins(c)
        r = NaN;
        idx = find(c >= 1 & c <= 6);
        if numel(idx) < 10, return; end
        P = S.centre(idx, :);
        D = sqrt((P(:, 1) - P(:, 1)').^2 + (P(:, 2) - P(:, 2)').^2);
        D(1:numel(idx) + 1:end) = Inf;
        [~, o] = sort(D, 2);
        nn = o(:, 1:3);
        cc = c(idx);
        meme = mean(mean(cc(nn) == cc));
        p = histc(cc, 1:6) / numel(cc); %#ok<HISTC>
        r = meme / sum(p .^ 2);
    end


%% ======================= cellule =======================
    function clic(e)
        if ~isfield(S, 'n'), return; end
        switch fig.SelectionType
            case 'alt'                                   % clic droit : déplacer l'image
                p = ax.CurrentPoint(1, 1:2);
                glisse = p;
                fig.WindowButtonMotionFcn = @(~, ~) deplacer();
                fig.WindowButtonUpFcn = @(~, ~) fin_deplacement();
                return
            case 'open'                                  % double-clic : image entière
                vue.xlim = []; vue.cle = '';
                dessiner();
                return
        end
        x = round(e.IntersectionPoint(1)); y = round(e.IntersectionPoint(2));
        m = 0;
        if x >= 1 && x <= S.Lx && y >= 1 && y <= S.Ly, m = double(S.label(y, x)); end
        if m == 0
            % pas sur un masque : la cellule la plus proche à moins de 8 pixels
            d = hypot(S.centre(:, 1) - x, S.centre(:, 2) - y);
            [dm, k] = min(d);
            if dm <= 8, m = k; end
        end
        selection = m;
        dessiner();
        afficher_cellule();
        dessiner_triangle();
    end

    function ph = position_hypothetique(m)
        % [x y r g b maxRGB] d'une cellule rejetée, avec la correction de fuite, le fond et le poids
        % du bleu de MaxwellAuto, mais sans le seuil de rejet. Vide si les valeurs manquent.
        ph = [];
        if isfield(S.P, 'roi_red') && isfield(S.P, 'params_fuite') && isfield(S.P, 'fond_image')
            V = [S.P.roi_red(m), S.P.roi_green(m), S.P.roi_blue(m)]; ca = S.P.roi_ca(m);
            pf = S.P.params_fuite;
            if isfield(pf, 'k_final'), k = pf.k_final; b = pf.b_final; else, k = pf.k; b = pf.b; end
            fo = S.P.fond_image;
        else
            f = fullfile(S.dossier, 'couleurs_details_Auto.mat');
            if ~isfile(f), return; end
            if ~isfield(tri, 'valeurs') || ~strcmp(tri.valeurs_cle, S.nom)
                Z = load(f, 'roi_red', 'roi_green', 'roi_blue', 'roi_ca', 'params_fuite', 'fond_image');
                if ~all(isfield(Z, {'roi_red', 'params_fuite', 'fond_image'})), return; end
                tri.valeurs = Z; tri.valeurs_cle = S.nom;
            end
            Z = tri.valeurs;
            if numel(Z.roi_red) ~= S.n, return; end
            V = [Z.roi_red(m), Z.roi_green(m), Z.roi_blue(m)]; ca = Z.roi_ca(m);
            k = Z.params_fuite.k; b = Z.params_fuite.b; fo = Z.fond_image;
        end
        % k et b dans l'ordre vert rouge bleu ; V dans l'ordre rouge vert bleu.
        % Demande de l'utilisateur (17 nov.) : les rapports de la cellule elle-même, sans retirer
        % le fond (qui écrasait les canaux faibles à 0) : valeurs corrigées de la fuite du GCaMP,
        % poids du bleu 0,7 ; si tout est à 0 après la fuite, valeurs brutes.
        net = max(0, V - [k(2) * ca + b(2), k(1) * ca + b(1), k(3) * ca + b(3)]) .* [1 1 0.7];
        if sum(net) <= 0, net = max(0, V) .* [1 1 0.7]; end
        if sum(net) <= 0, return; end
        p = net / sum(net);
        pb = max(0, V) / max(sum(max(0, V)), eps);         % rapports bruts (avant fuite), pour info
        ph = [p(1) + p(3) / 2, p(3) * sqrt(3) / 2, p, max(V), pb];
    end

    function dessiner_triangle()
        % petit triangle de Maxwell : toutes les cellules (couleur du fichier « Colours of »),
        % la cellule sélectionnée entourée (17 nov.)
        cla(axT); hold(axT, 'on');
        h = sqrt(3) / 2;
        plot(axT, [0 1 0.5 0], [0 0 h 0], 'k-', 'LineWidth', 1);
        text(axT, 1.02, -0.03, 'R', 'Color', [0.8 0 0], 'FontWeight', 'bold');
        text(axT, -0.07, -0.03, 'G', 'Color', [0 0.6 0], 'FontWeight', 'bold');
        text(axT, 0.47, h + 0.05, 'B', 'Color', [0 0 0.8], 'FontWeight', 'bold');
        axis(axT, 'equal'); axT.XLim = [-0.1 1.1]; axT.YLim = [-0.08 h + 0.1];
        lbl_tri.Text = '';
        if ~isfield(S, 'n'), hold(axT, 'off'); return; end
        cle = sprintf('%s|%s|%d|%s', S.nom, dd_affiche.Value, isfield(S.P, 'Donnees_Clustering'), dd_image.Value);
        if ~strcmp(cle, tri.cle)
            try
                [tri.D, tri.source] = coordonnees_triangle();
            catch
                tri.D = []; tri.source = 'none';
            end
            tri.cle = cle;
        end
        D = tri.D;
        if isempty(D) || size(D, 1) ~= S.n
            lbl_tri.Text = 'No triangle data for this session.';
            hold(axT, 'off'); return
        end
        c = couleurs_de(dd_affiche.Value);
        if isempty(c), c = 8 * ones(S.n, 1); end
        c = double(c(:)); c(~ismember(c, 1:8)) = 8;
        PAL = [1 0 0; 0 0.8 0; 0 0 1; 0.9 0.9 0; 1 0 1; 0 0.8 0.8; 0.85 0.85 0.85; 0.6 0.6 0.6];
        noire = c == 8;
        scatter(axT, D(noire, 1), D(noire, 2), 8, PAL(8, :), 'filled', 'HitTest', 'off');
        scatter(axT, D(~noire, 1), D(~noire, 2), 10, PAL(c(~noire), :), 'filled', 'HitTest', 'off');
        txt = sprintf('Triangle: %s', tri.source);
        if selection >= 1 && selection <= S.n
            m = selection;
            if any(isnan(D(m, :)))
                % cellule rejetée (trop faible) : où elle tomberait sans le rejet (17 nov.)
                ph = position_hypothetique(m);
                if isempty(ph)
                    txt = sprintf('Cell %d: rejected (too weak), not in the triangle (%s)\n%s', m, nom_c(c(m)), txt);
                else
                    plot(axT, ph(1), ph(2), 'o', 'MarkerSize', 14, 'LineWidth', 2, 'LineStyle', 'none', ...
                        'Color', [0.3 0.3 0.3], 'HitTest', 'off');
                    plot(axT, ph(1), ph(2), 'x', 'MarkerSize', 8, 'LineWidth', 1.5, 'Color', [0.3 0.3 0.3], 'HitTest', 'off');
                    txt = sprintf(['Cell %d: rejected (too weak, max R,G,B %.2f)\nwould be at r %.2f  g %.2f  b %.2f ' ...
                        '(grey cross, GCaMP leak removed)\nraw ratios r %.2f  g %.2f  b %.2f\n%s'], ...
                        m, ph(6), ph(3), ph(4), ph(5), ph(7), ph(8), ph(9), txt);
                end
            else
                plot(axT, D(m, 1), D(m, 2), 'o', 'MarkerSize', 13, 'LineWidth', 2.5, 'Color', 'k', 'HitTest', 'off');
                plot(axT, D(m, 1), D(m, 2), 'o', 'MarkerSize', 17, 'LineWidth', 1, 'Color', [1 0.6 0], 'HitTest', 'off');
                b = D(m, 2) / h; r = D(m, 1) - 0.5 * b; g = 1 - r - b;
                txt = sprintf('Cell %d: r %.2f  g %.2f  b %.2f (%s)\n%s', m, r, g, b, nom_c(c(m)), txt);
            end
        end
        lbl_tri.Text = txt;
        hold(axT, 'off');
    end

    function molette(e)
        % zoom autour du pointeur, sans passer par les outils
        if ~isfield(S, 'n') || isempty(ax.Children), return; end
        p = ax.CurrentPoint(1, 1:2);
        xl = ax.XLim; yl = ax.YLim;
        if p(1) < xl(1) || p(1) > xl(2) || p(2) < yl(1) || p(2) > yl(2), return; end
        f = 1.25 ^ e.VerticalScrollCount;               % molette vers soi : zoom arrière
        im = findobj(ax, 'Type', 'image');
        lx = [0.5, size(im(1).CData, 2) + 0.5]; ly = [0.5, size(im(1).CData, 1) + 0.5];
        nx = p(1) + (xl - p(1)) * f; ny = p(2) + (yl - p(2)) * f;
        if diff(nx) >= diff(lx) || diff(ny) >= diff(ly)
            nx = lx; ny = ly;
        else
            nx = nx - min(0, nx(1) - lx(1)) - max(0, nx(2) - lx(2));
            ny = ny - min(0, ny(1) - ly(1)) - max(0, ny(2) - ly(2));
        end
        ax.XLim = nx; ax.YLim = ny;
        vue.xlim = nx; vue.ylim = ny;
    end

    function deplacer()
        if isempty(glisse), return; end
        p = ax.CurrentPoint(1, 1:2);
        d = glisse - p;
        ax.XLim = ax.XLim + d(1); ax.YLim = ax.YLim + d(2);
    end

    function fin_deplacement()
        glisse = [];
        fig.WindowButtonMotionFcn = @(~, ~) rien();
        fig.WindowButtonUpFcn = '';
        vue.xlim = ax.XLim; vue.ylim = ax.YLim;
    end

    function afficher_cellule()
        if ~isfield(S, 'n') || selection < 1
            lbl_cellule.Text = 'none';
            return
        end
        m = selection;
        L = {sprintf('Cell %d (rank among iscell = 1 = row of colorcell)', m)};
        for k = 1:numel(S.fichiers)
            c = S.fichiers(k).colorcell;
            if numel(c) == S.n, L{end+1} = sprintf('  %s : %s', S.fichiers(k).nom, nom_c(c(m))); end %#ok<AGROW>
        end
        if ~isempty(S.travail), L{end+1} = sprintf('  working copy : %s', nom_c(S.travail(m))); end
        if isfield(S.P, 'roi_red')
            L{end+1} = sprintf('  raw R G B : %.3g  %.3g  %.3g', S.P.roi_red(m), S.P.roi_green(m), S.P.roi_blue(m));
        end
        if isfield(S.P, 'roi_red_corrected3')
            L{end+1} = sprintf('  corrected R G B : %.3g  %.3g  %.3g', S.P.roi_red_corrected3(m), ...
                S.P.roi_green_corrected3(m), S.P.roi_blue_corrected3(m));
        end
        if isfield(S.P, 'r')
            L{end+1} = sprintf('  proportions r g b : %.2f  %.2f  %.2f ; margin %.3f', S.P.r(m), S.P.g(m), ...
                S.P.b(m), S.P.marge(m));
        end
        if isfield(S.P, 'id_low_intensity') && numel(S.P.id_low_intensity) >= m && S.P.id_low_intensity(m)
            L{end+1} = '  below the low-intensity threshold';
        end
        lbl_cellule.Text = strjoin(L, newline);
    end

    function copie_travail()
        c = couleurs_de(dd_affiche.Value);
        if isempty(c)
            uialert(fig, 'Show a colour file with the right number of cells first.', 'Working copy');
            return
        end
        if ~isempty(S.travail) && S.modifie
            r = uiconfirm(fig, 'Replace the current working copy?', 'Working copy', ...
                'Options', {'Replace', 'Cancel'}, 'DefaultOption', 2);
            if strcmp(r, 'Cancel'), return; end
        end
        S.travail = c;
        S.source_travail = dd_affiche.Value;
        S.modifie = true;
        historique = {};
        maj_listes();
        dd_affiche.Value = TRAVAIL;
        changer_affichage();
        etat_boutons();
    end

    function mettre_couleur(c)
        if ~isfield(S, 'n') || selection < 1, return; end
        if isempty(S.travail)
            uialert(fig, ['No working copy: run step 7, or show a file and click ' ...
                '"Edit a copy of the displayed file".'], 'Working copy');
            return
        end
        historique{end+1} = [selection, S.travail(selection)];
        S.travail(selection) = c;
        S.modifie = true;
        changer_affichage();
    end

    function annuler()
        if isempty(historique), return; end
        h = historique{end}; historique(end) = [];
        S.travail(h(1)) = h(2);
        selection = h(1);
        changer_affichage();
    end


%% ======================= traitement =======================
    function lancer_1_3()
        for k = 1:3
            if ~lancer_etape(k), return; end
        end
    end

    function ok = lancer_etape(k)
        ok = false;
        if ~isfield(S, 'n'), uialert(fig, 'Choose a session first.', 'Processing'); return; end
        pathload = [regexprep(S.dossier, '[\\/]+$', '') filesep];
        d = [];
        if ~ismember(k, [4 5 6 8 9])
            d = uiprogressdlg(fig, 'Title', noms_etapes{k}, 'Indeterminate', 'on');
        end
        try
            P = S.P;
            switch k
                case 1
                    [red, red890, green, blue, calcium] = brainbow_image_read(pathload);
                    P.images_norm = brainbow_image_processing(pathload, red, red890, green, blue, calcium);
                case 2
                    [P.aligned_red, P.aligned_green, P.aligned_blue, P.calcium_norm, P.aligned_red890, ...
                        P.info_alignement] = image_alignment_multimodal(P.images_norm);
                case 3
                    [P.roi_red, P.roi_green, P.roi_blue, P.roi_ca, ~, ~, ~, ~, P.outline_gcampx, ...
                        P.outline_gcampy, P.extra_roi] = ROI_extraction2(pathload, P.aligned_red, P.aligned_green, ...
                        P.aligned_blue, P.images_norm.calcium);
                    if numel(P.roi_red) ~= S.n
                        error('ROI_extraction2 gives %d cells, Fall.mat %d (iscell = 1).', numel(P.roi_red), S.n);
                    end
                case 4
                    % proposition de départ (méthode MaxwellAuto, 17 nov.) : max(R, V, B) < 0,20 ;
                    % on peut ensuite cliquer des cellules ou bouger le curseur
                    id_prev = find(max([P.roi_red(:), P.roi_green(:), P.roi_blue(:)], [], 2) < 0.20);
                    [P.roi_red_corrected, P.roi_green_corrected, P.roi_blue_corrected, P.id_low_intensity] = ...
                        interactive_threshold_selection(P.aligned_red, P.aligned_green, P.aligned_blue, ...
                        P.calcium_norm, P.roi_red, P.roi_green, P.roi_blue, P.outline_gcampx, P.outline_gcampy, ...
                        ef_seuil.Value, id_prev);
                    id_low_intensity = P.id_low_intensity; %#ok<NASGU>
                    save(fullfile(pathload, 'id_low_intensity.mat'), 'id_low_intensity');
                case 5
                    P.limits_k = [0.4, 0.2, 0.1];       % comme color_predictor_model : [vert, rouge, bleu]
                    [P.roi_green_corrected2, P.roi_red_corrected2, P.roi_blue_corrected2, ...
                        P.aligned_green_corrected, P.aligned_red_corrected, P.aligned_blue_corrected, P.params_fuite] = ...
                        bleedthrough_correction3(P.roi_green_corrected, P.roi_red_corrected, P.roi_blue_corrected, ...
                        P.roi_ca, P.limits_k, P.aligned_green, P.aligned_red, P.aligned_blue, P.calcium_norm);
                case 6
                    % proposition de départ (méthode MaxwellAuto) : cellules de l'étape 4,
                    % plus max(R, V, B) corrigés < 0,15
                    P.id_low_intensity = union(P.id_low_intensity(:), find(max([P.roi_red_corrected2(:), ...
                        P.roi_green_corrected2(:), P.roi_blue_corrected2(:)], [], 2) < 0.15));
                    [P.roi_red_corrected3, P.roi_green_corrected3, P.roi_blue_corrected3, P.id_low_intensity] = ...
                        interactive_threshold_selection(P.aligned_red_corrected, P.aligned_green_corrected, ...
                        P.aligned_blue_corrected, P.calcium_norm, P.roi_red_corrected2, P.roi_green_corrected2, ...
                        P.roi_blue_corrected2, P.outline_gcampx, P.outline_gcampy, ef_seuil.Value, P.id_low_intensity);
                    id_low_intensity = P.id_low_intensity; %#ok<NASGU>
                    save(fullfile(pathload, 'id_low_intensity.mat'), 'id_low_intensity');
                case 7
                    % méthode MaxwellAuto (17 nov.) : échelle 'ancien', pas de classe blanche,
                    % noire si la couleur est proche d'une frontière (marge < marge_min)
                    % (calcul commun avec colorcell_maxwellauto : fond de l'image hors ROI,
                    % poids du bleu 0,7, frontières [0,4 0,3 0,4] ; ep n'est plus utilisé)
                    P.ep = 0; P.seuil_blanc = Inf; P.marge_min = ef_marge.Value;
                    rej = false(S.n, 1); rej(P.id_low_intensity) = true;
                    [P.Donnees_Clustering, P.r, P.g, P.b, P.classe_triangle, P.marge, P.info_maxwell] = ...
                        proportions_auto(P.roi_red_corrected3, P.roi_green_corrected3, P.roi_blue_corrected3, ...
                        {P.aligned_red_corrected, P.aligned_green_corrected, P.aligned_blue_corrected}, ...
                        S.dossier, rej | ~isfinite(P.roi_red_corrected3(:)));
                    P.fond_image = P.info_maxwell.fond;
                    P.colorcell_classification = P.classe_triangle;
                    P.colorcell_classification(isnan(P.marge) | P.marge < P.marge_min) = 8;
                    % 17 nov. : masques trop gros (> 2,5 x la médiane de la session) ou trop allongés (axes > 2,5) en noir
                    [gros, P.masques_trop_gros] = masques_trop_gros(S.dossier);
                    P.colorcell_classification(gros) = 8;
                    if ~isempty(S.travail) && S.modifie
                        r = uiconfirm(fig, 'Replace the working copy by the new classification?', 'Step 7', ...
                            'Options', {'Replace', 'Keep the working copy'}, 'DefaultOption', 1);
                    else
                        r = 'Replace';
                    end
                    if strcmp(r, 'Replace')
                        S.travail = P.colorcell_classification(:);
                        S.source_travail = 'step 7 (automatic classification)';
                        S.modifie = true;
                        historique = {};
                    end
                case 8
                    c8 = S.travail(:);
                    c8(~ismember(c8, 1:8)) = 8;         % 0 ou NaN (anciens fichiers) -> aucune, sinon palette(...) plante
                    triangle_maxwell(c8, P.Donnees_Clustering, 0, P.id_low_intensity);
                case 9
                    S.travail = interactive_color_correction(P.aligned_red, P.aligned_green, P.aligned_blue, ...
                        P.calcium_norm, S.travail, P.outline_gcampx, P.outline_gcampy);
                    S.travail = S.travail(:);
                    S.modifie = true;
                case 10
                    enregistrer(pathload);
            end
            S.P = P;
            S.etapes(k) = true;
            ok = true;
        catch err
            if ~isempty(d) && isvalid(d), close(d); end
            uialert(fig, sprintf('%s\n\n(%s)', err.message, err.identifier), noms_etapes{k});
            figure(fig);
            return
        end
        if ~isempty(d) && isvalid(d), close(d); end
        figure(fig);
        maj_listes();
        if k == 2, dd_image.Value = 'aligned (step 2)'; end
        if k == 5, dd_image.Value = 'corrected (step 5)'; end
        if ismember(k, [7 9]), dd_affiche.Value = TRAVAIL; end
        changer_affichage();
        etat_boutons();
    end

    function enregistrer(pathload)
        if isempty(S.travail) || numel(S.travail) ~= S.n
            error('No working copy with %d cells to save.', S.n);
        end
        % 17 nov. : enregistrement dans colorcell_MaxwellAuto.mat (fichier lu par le pipeline via la
        % colonne fichier_colorcell de la table) ; l'ancien et ses détails sont gardés avec la date
        f = fullfile(pathload, 'colorcell_MaxwellAuto.mat');
        fd = fullfile(pathload, 'couleurs_details_Auto.mat');
        date_f = datestr(now, 'yyyy_mm_dd_HHMMSS');
        details_avant = struct();
        if isfile(fd), details_avant = load(fd); end
        if isfile(f)
            r = uiconfirm(fig, sprintf(['colorcell_MaxwellAuto.mat already exists in this folder. It will be ' ...
                'renamed colorcell_MaxwellAuto_<date>.mat and replaced. Continue?']), 'Save', ...
                'Options', {'Save', 'Cancel'}, 'DefaultOption', 2);
            if strcmp(r, 'Cancel'), error('Saving cancelled.'); end
            movefile(f, fullfile(pathload, ['colorcell_MaxwellAuto_' date_f '.mat']));
        end
        if isfile(fd), movefile(fd, fullfile(pathload, ['couleurs_details_Auto_' date_f '.mat'])); end
        colorcell = S.travail(:); %#ok<NASGU>
        save(f, 'colorcell');
        % détails : tout ce qui a servi ; si les étapes n'ont pas été refaites (simple correction
        % à la main d'un fichier existant), on garde les détails d'avant (valeurs, triangle...)
        P = S.P;
        details = struct();
        if ~isfield(P, 'Donnees_Clustering'), details = details_avant; end
        details.date = datestr(now, 'yyyy-mm-dd HH:MM:SS');
        details.session = S.nom;
        details.colorcell = S.travail(:);
        details.source = '';
        if isfield(S, 'source_travail'), details.source = S.source_travail; end
        champs = {'colorcell_classification', 'id_low_intensity', 'roi_red', 'roi_green', 'roi_blue', 'roi_ca', ...
            'roi_red_corrected3', 'roi_green_corrected3', 'roi_blue_corrected3', 'r', 'g', 'b', 'marge', ...
            'Donnees_Clustering', 'info_maxwell', 'info_alignement', 'params_fuite', 'extra_roi', 'ep', 'seuil_blanc', ...
            'limits_k', 'marge_min', 'classe_triangle', 'fond_image', 'masques_trop_gros'};
        for z = 1:numel(champs)
            if isfield(P, champs{z}), details.(champs{z}) = P.(champs{z}); end
        end
        if isfield(P, 'colorcell_classification')
            details.cellules_corrigees = find(S.travail(:) ~= P.colorcell_classification(:));
        elseif isfield(details_avant, 'colorcell') && numel(details_avant.colorcell) == S.n
            details.cellules_corrigees_a_la_main = find(S.travail(:) ~= details_avant.colorcell(:));
        end
        save(fd, '-struct', 'details');
        S.modifie = false;
        % recharger la liste des fichiers
        Lf = dir(fullfile(S.dossier, 'colorcell*.mat'));
        S.fichiers = struct('nom', {}, 'colorcell', {});
        for k = 1:numel(Lf)
            try
                C = load(fullfile(Lf(k).folder, Lf(k).name), 'colorcell');
                S.fichiers(end+1) = struct('nom', Lf(k).name, 'colorcell', double(C.colorcell(:)));
            catch
            end
        end
        remplir_liste();
        lbl_etat.Text = sprintf('Saved: %s (%d cells).', f, S.n);
    end

    function etat_boutons()
        session = isfield(S, 'n');
        e = session & [true, false(1, 9)];
        if session
            P = S.P;
            e(2) = isfield(P, 'images_norm');
            e(3) = isfield(P, 'aligned_red');
            e(4) = isfield(P, 'roi_red');
            e(5) = isfield(P, 'roi_red_corrected');
            e(6) = isfield(P, 'roi_red_corrected2');
            e(7) = isfield(P, 'roi_red_corrected3');
            e(8) = isfield(P, 'Donnees_Clustering') && ~isempty(S.travail);
            e(9) = isfield(P, 'aligned_red') && ~isempty(S.travail);
            e(10) = ~isempty(S.travail);
        end
        for k = 1:10
            b_etape(k).Enable = onoff(e(k));
            if session && S.etapes(k), lbl_ok(k).Text = char(10003); else, lbl_ok(k).Text = ''; end
        end
        b_123.Enable = onoff(session);
        b_edit.Enable = onoff(session);
        b_fig.Enable = onoff(session);
        if session
            t = sprintf('Folder: %s', S.dossier);
            if ~isempty(S.travail)
                t = sprintf('%s\nWorking copy from: %s', t, S.source_travail);
                if S.modifie, t = [t ' (not saved)']; end
            end
            lbl_etat.Text = t;
        end
    end

    function t = nom_c(c)
        if c >= 1 && c <= 8, t = NOMS{c}; else, t = sprintf('%g', c); end
    end
end


%% ======================= fonctions =======================
function sessions = lire_sessions(fichier)
% Les sessions de la table à analyser (inclure = 1 depuis le 17 nov.), toutes conditions
opts = detectImportOptions(fichier);
texte = intersect({'id_animal', 'id_champ', 'id_session', 'id_condition', 'dossier_suite2p', 'fichier_colorcell'}, ...
    opts.VariableNames, 'stable');
if ~isempty(texte), opts = setvartype(opts, texte, 'char'); end
T = readtable(fichier, opts);
requis = {'id_animal', 'id_session', 'id_champ', 'dossier_suite2p'};
manque = setdiff(requis, T.Properties.VariableNames);
if ~isempty(manque), error('Session table: missing column(s) %s.', strjoin(manque, ', ')); end
if ~ismember('id_condition', T.Properties.VariableNames), T.id_condition = repmat({''}, height(T), 1); end
if ~ismember('fichier_colorcell', T.Properties.VariableNames), T.fichier_colorcell = repmat({''}, height(T), 1); end
T.fichier_colorcell(cellfun(@(x) ~ischar(x), T.fichier_colorcell)) = {''};
dossier = cell(height(T), 1);
for s = 1:height(T)
    d = strtrim(char(T.dossier_suite2p{s}));
    [parent, ~, ext] = fileparts(d);
    if strcmpi(ext, '.mat'), d = parent; end
    dossier{s} = d;
end
garder = ~cellfun(@isempty, dossier);
% 17 nov. (demandé par l'utilisateur) : seulement les sessions à analyser (inclure = 1)
if ismember('inclure', T.Properties.VariableNames)
    garder = garder & T.inclure == 1;
end
T = T(garder, :); dossier = dossier(garder);
nom = strcat(T.id_animal, '_', T.id_session, '_', T.id_champ, '_', T.id_condition);
nom = regexprep(nom, '_$', '');
sessions = table(nom, dossier, T.fichier_colorcell, 'VariableNames', {'nom', 'dossier', 'fichier_colorcell'});
sessions = sortrows(sessions, 'nom');
end


function f = trouver_composite_local(dossier)
% compositecalcium.tif (rouge, vert, bleu, calcium), sinon un fichier
% « composite », dans le dossier puis son parent (comme trouver_composite)
f = '';
d = regexprep(char(dossier), '[\\/]+$', '');
for niveau = 1:2
    % compositecalcium_V2.tif d'abord : recalage actuel sur le GCaMP (refaire_images_recalees, 17 nov.) ;
    % l'ancien compositecalcium.tif vient du recalage de 2023-2025
    for motif = {'compositecalcium_V2.tif', 'compositecalcium.tif', 'compositecalcium.tiff', '*omposite*.tif', '*omposite*.tiff'}
        L = dir(fullfile(d, motif{1}));
        if ~isempty(L)
            f = fullfile(L(1).folder, L(1).name);
            return
        end
    end
    d = fileparts(d);
end
end


function canaux = lire_composite(f)
% une page par canal (ImageJ), ou les plans d'une image RVB ; 4 au plus
info = imfinfo(f);
canaux = {};
for pg = 1:numel(info)
    I = double(imread(f, pg));
    for c = 1:size(I, 3), canaux{end+1} = I(:, :, c); end %#ok<AGROW>
end
canaux = canaux(1:min(4, numel(canaux)));
end


function rien()
end


function t = onoff(vrai)
if vrai, t = 'on'; else, t = 'off'; end
end


function v = ifelse(test, a, b)
if test, v = a; else, v = b; end
end


%% ---------------------------------------------------------------------
function [cx, cy] = contour_cellule(x, y)
% contour continu (fermé) d'une ROI à partir de ses pixels de bord (xext, yext de suite2p) :
% on remplit l'anneau de pixels puis on suit son bord avec bwboundaries. Sans la boîte à outils
% Image Processing : enveloppe convexe (comme interactive_color_correction).
cx = x(:); cy = y(:);
if numel(x) < 3, return; end
x0 = min(x) - 3; y0 = min(y) - 3;
M = false(max(y) - y0 + 3, max(x) - x0 + 3);
M(sub2ind(size(M), y - y0, x - x0)) = true;
try
    M = imfill(imclose(M, true(3)), 'holes');        % referme les trous d'un pixel dans l'anneau
    B = bwboundaries(M, 8, 'noholes');
    [~, k] = max(cellfun(@(b) size(b, 1), B));
    cx = B{k}(:, 2) + x0; cy = B{k}(:, 1) + y0;
catch
    k = convhull(x(:), y(:));
    cx = x(k); cx = cx(:); cy = y(k); cy = cy(:);
end
end


function [X, Y] = joindre(cx, cy, m)
% tous les contours m en une seule ligne, séparés par NaN (un seul objet graphique : rapide)
X = cellfun(@(v) [v(:); NaN], cx(m), 'UniformOutput', false);
Y = cellfun(@(v) [v(:); NaN], cy(m), 'UniformOutput', false);
X = vertcat(X{:}); Y = vertcat(Y{:});
end
