function fig = visionneuse_composite(dossier_s2p, nom_session, fichier_colorcell, fichier_image)
%VISIONNEUSE_COMPOSITE  Fenêtre pour vérifier les couleurs d'un champ sur
%   l'image Brainbow (10 nov., demande de l'utilisateur).
%   Image : compositecalcium.tif (11 nov. : prioritaire), sinon un fichier
%   dont le nom contient « composite », dans le dossier des données (sinon son
%   dossier parent) ; voir trouver_composite, ou fichier_image s'il est donné ; bouton « Other file... » sinon.
%   Image ImageJ à plusieurs canaux (une page par canal) ou image RVB.
%   Pour chaque canal : affiché ou non, couleur (rouge, vert, bleu, gris...),
%   gain. Contraste de chaque canal : 1er au 99,7e percentile, puis gain.
%   Masques GCaMP (contours suite2p, cellules iscell = 1) : affichés ou non,
%   en blanc ou dans la couleur attribuée (colorcell, lue brute), avec leur
%   numéro parmi les cellules iscell = 1 (= ligne de colorcell et du Raster ;
%   11 nov. : plus de numéro de ROI suite2p, à la demande de l'utilisateur).
%   « Highlight ROI » : entoure une cellule (même numéro) en jaune.
%   Zoom et déplacement : barre d'outils de l'image.
%   Ne modifie rien.
%
%   Ex. : visionneuse_composite('E:\Data\Aurelie\data\nocues\444175\221121_plane1\')
%   Appelé par le bouton « Composite viewer... » de interface_brainbow.

if nargin < 2 || isempty(nom_session), nom_session = ''; end
if nargin < 3, fichier_colorcell = ''; end
if nargin < 4, fichier_image = ''; end
dossier_s2p = avec_separateur(char(dossier_s2p));

% Géométrie (mêmes cellules et même ordre que le pipeline)
S2P = load([dossier_s2p 'Fall.mat'], 'stat', 'iscell');
roi = find(S2P.iscell(:, 1) == 1);
n = numel(roi);
bx = cell(1, n); by = cell(1, n); centre = zeros(n, 2);
for m = 1:n
    if iscell(S2P.stat), st = S2P.stat{roi(m)}; else, st = S2P.stat(roi(m)); end
    bx{m} = double(st.xext(:)) + 1;
    by{m} = double(st.yext(:)) + 1;
    centre(m, :) = [mean(double(st.xpix(:)) + 1), mean(double(st.ypix(:)) + 1)];
end
taille_masques = [max(cellfun(@max, by)), max(cellfun(@max, bx))];
try
    O = load([dossier_s2p 'Fall.mat'], 'ops');
    taille_masques = [O.ops.Ly, O.ops.Lx];
catch
end

% Couleurs attribuées (brutes)
colorcell = [];
try
    colorcell = charger_couleurs(dossier_s2p, fichier_colorcell);
    colorcell = colorcell(:);
    if numel(colorcell) ~= n, colorcell = []; end
catch
end
rgb_bb = [1 0.15 0.15; 0.15 1 0.15; 0.25 0.5 1; 1 1 0.1; 1 0.2 1; 0.1 1 1; 1 1 1; 0.6 0.6 0.6];

% Fenêtre
couleurs_noms = {'Red', 'Green', 'Blue', 'Grey', 'Cyan', 'Magenta', 'Yellow'};
couleurs_rgb = [1 0 0; 0 1 0; 0 0 1; 1 1 1; 0 1 1; 1 0 1; 1 1 0];
fig = uifigure('Name', ['Composite viewer - ' nom_session], 'Position', [80 60 1250 900]);
g = uigridlayout(fig, [1 2]);
g.ColumnWidth = {'1x', 300};
ax = uiaxes(g);
ax.Layout.Row = 1; ax.Layout.Column = 1;
p = uigridlayout(g, [16 3]);
p.Layout.Row = 1; p.Layout.Column = 2;
p.RowHeight = [{22, 40, 24, 18}, repmat({24}, 1, 4), {10, 22, 22, 22, 22, 24, 10, '1x'}];
p.ColumnWidth = {70, 110, '1x'};
uilabel(p, 'Text', 'Image', 'FontWeight', 'bold');
lbl_fichier = uilabel(p, 'Text', '', 'WordWrap', 'on', 'FontSize', 10);
lbl_fichier.Layout.Row = 2; lbl_fichier.Layout.Column = [1 3];
b_fichier = uibutton(p, 'Text', 'Other file...', 'ButtonPushedFcn', @(~, ~) choisir_fichier());
b_fichier.Layout.Row = 3; b_fichier.Layout.Column = [1 2];
e1 = uilabel(p, 'Text', 'Channel'); e1.Layout.Row = 4; e1.Layout.Column = 1;
e2 = uilabel(p, 'Text', 'Colour'); e2.Layout.Row = 4; e2.Layout.Column = 2;
e3 = uilabel(p, 'Text', 'Gain'); e3.Layout.Row = 4; e3.Layout.Column = 3;
cb_canal = gobjects(1, 4); dd_canal = gobjects(1, 4); sp_canal = gobjects(1, 4);
defaut_coul = {'Red', 'Green', 'Blue', 'Grey'};   % ordre de compositecalcium.tif (color_predictor_model) : rouge, vert, bleu, calcium
for k = 1:4
    cb_canal(k) = uicheckbox(p, 'Text', sprintf('Ch %d', k), 'Value', k < 4, ...
        'ValueChangedFcn', @(~, ~) dessiner());
    cb_canal(k).Layout.Row = 4 + k; cb_canal(k).Layout.Column = 1;
    dd_canal(k) = uidropdown(p, 'Items', couleurs_noms, 'Value', defaut_coul{k}, ...
        'ValueChangedFcn', @(~, ~) dessiner());
    dd_canal(k).Layout.Row = 4 + k; dd_canal(k).Layout.Column = 2;
    sp_canal(k) = uispinner(p, 'Value', 1, 'Limits', [0.1 10], 'Step', 0.1, ...
        'ValueChangedFcn', @(~, ~) dessiner());
    sp_canal(k).Layout.Row = 4 + k; sp_canal(k).Layout.Column = 3;
end
cb_masque = uicheckbox(p, 'Text', 'GCaMP masks', 'Value', true, 'ValueChangedFcn', @(~, ~) dessiner());
cb_masque.Layout.Row = 10; cb_masque.Layout.Column = [1 3];
cb_coul = uicheckbox(p, 'Text', 'Masks in assigned colour (colorcell)', 'Value', false, ...
    'ValueChangedFcn', @(~, ~) dessiner());
cb_coul.Layout.Row = 11; cb_coul.Layout.Column = [1 3];
if isempty(colorcell)
    cb_coul.Enable = 'off';
    cb_coul.Tooltip = 'colorcell not found or not the same number of cells as Fall.mat';
end
cb_num = uicheckbox(p, 'Text', 'Cell numbers (among iscell = 1, = row of colorcell)', 'Value', false, ...
    'ValueChangedFcn', @(~, ~) dessiner());
cb_num.Layout.Row = 12; cb_num.Layout.Column = [1 3];
l_h = uilabel(p, 'Text', 'Highlight ROI');
l_h.Layout.Row = 13; l_h.Layout.Column = 1;
ef_h = uieditfield(p, 'numeric', 'Value', 0, 'Limits', [0 n], 'RoundFractionalValues', 'on', ...
    'Tooltip', 'Cell number among iscell = 1 (0 = none)', 'ValueChangedFcn', @(~, ~) dessiner());
ef_h.Layout.Row = 13; ef_h.Layout.Column = 2;
lbl_info = uilabel(p, 'Text', '', 'WordWrap', 'on', 'FontSize', 10, 'VerticalAlignment', 'top');
lbl_info.Layout.Row = 16; lbl_info.Layout.Column = [1 3];

canaux = {};
plages = zeros(0, 2);
if isempty(fichier_image), fichier_image = trouver_composite(dossier_s2p); end
if isempty(fichier_image)
    lbl_fichier.Text = 'No file with "composite" in its name found: use "Other file...".';
else
    charger(fichier_image);
end
dessiner();


    function choisir_fichier()
        [f, pth] = uigetfile({'*.tif;*.tiff', 'TIFF'}, 'Composite image', dossier_s2p);
        figure(fig);
        if isequal(f, 0), return; end
        charger(fullfile(pth, f));
        dessiner();
    end

    function charger(f)
        % une page par canal (ImageJ), ou les plans d'une image RVB
        info = imfinfo(f);
        canaux = {};
        for pg = 1:numel(info)
            I = double(imread(f, pg));
            for c = 1:size(I, 3), canaux{end+1} = I(:, :, c); end %#ok<AGROW>
        end
        canaux = canaux(1:min(4, numel(canaux)));
        plages = zeros(numel(canaux), 2);
        for c = 1:numel(canaux)
            v = canaux{c}(:);
            plages(c, :) = [prctile(v, 1), prctile(v, 99.7)];
        end
        for c = 1:4
            cb_canal(c).Enable = onoff(c <= numel(canaux));
        end
        [~, nf, ext] = fileparts(f);
        lbl_fichier.Text = sprintf('%s%s : %d channel(s), %d x %d', nf, ext, numel(canaux), ...
            size(canaux{1}, 1), size(canaux{1}, 2));
    end

    function dessiner()
        cla(ax);
        hold(ax, 'on');
        Ly = taille_masques(1); Lx = taille_masques(2);
        txt = {};
        if ~isempty(canaux)
            [Ly, Lx] = size(canaux{1});
            rgb = zeros(Ly, Lx, 3);
            for c = 1:numel(canaux)
                if ~cb_canal(c).Value, continue; end
                v = (canaux{c} - plages(c, 1)) / max(plages(c, 2) - plages(c, 1), eps) * sp_canal(c).Value;
                v = min(max(v, 0), 1);
                col = couleurs_rgb(strcmp(couleurs_noms, dd_canal(c).Value), :);
                rgb = rgb + v .* reshape(col, 1, 1, 3);
            end
            image(ax, min(rgb, 1));
            if Ly ~= taille_masques(1) || Lx ~= taille_masques(2)
                txt{end+1} = sprintf(['WARNING: image %d x %d, suite2p masks %d x %d : masks are ' ...
                    'drawn in pixels of the masks, probably not aligned.'], Ly, Lx, taille_masques(1), ...
                    taille_masques(2));
            end
        end
        if cb_masque.Value
            if cb_coul.Value && ~isempty(colorcell)
                for c = 1:8
                    m = find(colorcell == c);
                    if isempty(m), continue; end
                    plot(ax, vertcat(bx{m}), vertcat(by{m}), '.', 'Color', rgb_bb(c, :), 'MarkerSize', 3);
                end
            else
                plot(ax, vertcat(bx{:}), vertcat(by{:}), '.', 'Color', [1 1 1], 'MarkerSize', 2);
            end
        end
        if cb_num.Value
            text(ax, centre(:, 1), centre(:, 2), arrayfun(@num2str, 1:n, 'UniformOutput', false), ...
                'Color', [1 1 0.6], 'FontSize', 7, 'HorizontalAlignment', 'center', 'Clipping', 'on');
        end
        h = ef_h.Value;
        if h >= 1 && h <= n
            plot(ax, bx{h}, by{h}, '.', 'Color', [1 1 0], 'MarkerSize', 8);
            if isempty(colorcell), tc = '?'; else, tc = nom_couleur(colorcell(h)); end
            txt{end+1} = sprintf('Cell %d: colorcell = %s', h, tc);
        end
        hold(ax, 'off');
        ax.DataAspectRatio = [1 1 1];
        ax.YDir = 'reverse';
        ax.XLim = [0.5, Lx + 0.5]; ax.YLim = [0.5, Ly + 0.5];
        ax.XTick = []; ax.YTick = [];
        title(ax, sprintf('%s (%d cells)', nom_session, n), 'Interpreter', 'none', 'FontWeight', 'normal');
        txt{end+1} = ['Numbers = rank among the cells iscell = 1 (order of Fall.mat) = row of colorcell ' ...
            'and of the Raster. Colorcell colours: R G B RG RB GB, white = all, grey = none.'];
        lbl_info.Text = strjoin(txt, newline);
    end
end


function t = nom_couleur(c)
noms = {'R', 'G', 'B', 'RG', 'RB', 'GB', 'white', 'none'};
if c >= 1 && c <= 8, t = noms{c}; else, t = sprintf('%g', c); end
end


function t = onoff(vrai)
if vrai, t = 'on'; else, t = 'off'; end
end
