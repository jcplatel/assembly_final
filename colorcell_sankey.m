function colorcell_sankey(files, varargin)
%COLORCELL_SANKEY  Sankey des transitions de couleur entre plusieurs fichiers
%   colorcell + carte "chemin de couleur par cellule".
%
%   Chaque fichier contient un vecteur colorcell : ligne i = code couleur de
%   la cellule i (même cellule à la même ligne dans tous les fichiers).
%
%   Codes : 1 red, 2 green, 3 blue, 4 yellow, 5 magenta, 6 cyan, 7 white, 8 black
%           0, NaN ou autre valeur = cellule absente à cette date
%           (seules les cellules présentes dans TOUS les fichiers sont gardées)
%
%   Formats : .mat (toutes versions, y compris v7.3), .txt / .csv (1 valeur par ligne)
%
%   Utilisation :
%     colorcell_sankey                 % fenêtre de sélection des fichiers
%     colorcell_sankey({'colorcell_220919.mat','colorcell_220922.mat','colorcell_220923.mat'})
%     colorcell_sankey(files, 'Labels', {'Ex vivo','220919','220922','220923'}, ...
%                      'Out', 'resultats/maxwell')
%
%   Options :
%     'Labels' : nom de chaque colonne (défaut : date à 6 ou 8 chiffres du nom de fichier)
%     'Var'    : nom de la variable dans le .mat (défaut 'colorcell')
%     'Out'    : préfixe des fichiers de sortie (défaut 'colorcell')
%     'Title'  : titre du Sankey
%
%   Sorties : <Out>_sankey.png/.svg, <Out>_paths.png,
%             <Out>_cells.csv, <Out>_transitions.csv

% ---------------------------------------------------------------- couleurs
NAMES = {'red','green','blue','yellow','magenta','cyan','white','black'};
RGB   = [1 0 0; 0 0.5 0; 0 0 1; 1 0.84 0; 1 0 1; 0 1 1; 1 1 1; 0 0 0];
ORDER = [1 4 2 6 3 5 7 8];   % ordre des nœuds de haut en bas

% ---------------------------------------------------------------- entrées
if nargin < 1 || isempty(files)
    [f, d] = uigetfile({'*.mat;*.txt;*.csv', 'Fichiers colorcell'}, ...
        'Choisir les fichiers colorcell (ordre chronologique)', 'MultiSelect', 'on');
    if isequal(f, 0), return; end
    files = fullfile(d, cellstr(f));
end
files = cellstr(files);

p = inputParser;
addParameter(p, 'Labels', {});
addParameter(p, 'Var', 'colorcell');
addParameter(p, 'Out', 'colorcell');
addParameter(p, 'Title', '');
parse(p, varargin{:});
opt = p.Results;

nd = numel(files);
if nd < 2, error('Il faut au moins 2 fichiers.'); end

labels = cellstr(opt.Labels);
if isempty(labels)
    labels = cellfun(@default_label, files, 'UniformOutput', false);
end
if numel(labels) ~= nd
    error('''Labels'' doit avoir autant d''éléments que de fichiers.');
end

% ---------------------------------------------------------------- lecture
V = cell(1, nd);
for j = 1:nd
    V{j} = load_colorcell(files{j}, opt.Var);
end
lens = cellfun(@numel, V);
n = min(lens);
if any(lens ~= n)
    warning('Longueurs différentes (%s) -> on garde les %d premières lignes.', ...
        num2str(lens), n);
end
M = zeros(n, nd);
for j = 1:nd
    M(:, j) = V{j}(1:n);
end

valid   = all(ismember(M, 1:8), 2);   % NaN, 0, etc. = absent
cellIds = find(valid);                % numéro de ligne d'origine
M       = M(valid, :);
fprintf('%d cellules au total, %d présentes dans les %d fichiers.\n', n, numel(cellIds), nd);
if isempty(cellIds), error('Aucune cellule commune.'); end

out = char(opt.Out);
outDir = fileparts(out);
if ~isempty(outDir) && ~exist(outDir, 'dir'), mkdir(outDir); end

ttl = opt.Title;
if isempty(ttl)
    ttl = {sprintf('Color transitions across %d dates', nd), ...
           sprintf('Cells present at all %d dates: n=%d', nd, numel(cellIds))};
end

plot_sankey(M, labels, out, ttl, RGB, NAMES, ORDER);
plot_paths(M, cellIds, labels, out, RGB, ORDER);
write_tables(M, cellIds, labels, out, NAMES, ORDER);

fprintf('Écrit : %s_sankey.png/.svg, %s_paths.png, %s_cells.csv, %s_transitions.csv\n', ...
    out, out, out, out);
end


% ========================================================================
function v = load_colorcell(f, varName)
[~, ~, ext] = fileparts(f);
if strcmpi(ext, '.mat')
    S = load(f);
    if isfield(S, varName)
        v = S.(varName);
    else
        fn = fieldnames(S);
        if numel(fn) ~= 1
            error('%s : variable ''%s'' absente (variables : %s)', f, varName, strjoin(fn, ', '));
        end
        v = S.(fn{1});
    end
else
    v = readmatrix(f);
end
v = double(v(:));
end


function lab = default_label(f)
[~, name] = fileparts(f);
lab = regexp(name, '\d{8}|\d{6}', 'match', 'once');
if isempty(lab), lab = name; end
end


% ========================================================================
function plot_sankey(M, labels, out, ttl, RGB, NAMES, ORDER)
[nC, nd] = size(M);
w = 0.035;        % largeur des nœuds
gap = 0.012;      % espace entre nœuds
alpha = 0.35;     % transparence des rubans
xs = linspace(0, 1, nd);

cnt = zeros(8, nd);
for j = 1:nd
    for c = 1:8
        cnt(c, j) = sum(M(:, j) == c);
    end
end
maxNodes = max(sum(cnt > 0, 1));
scale = (1 - gap * (maxNodes - 1)) / nC;   % hauteur d'une cellule

% position (haut, bas) de chaque nœud, y = 1 en haut
top = nan(8, nd); bot = nan(8, nd);
for j = 1:nd
    y = 1;
    for c = ORDER
        if cnt(c, j) > 0
            top(c, j) = y;
            bot(c, j) = y - cnt(c, j) * scale;
            y = bot(c, j) - gap;
        end
    end
end

wIn = 4.2 * nd + 2; hIn = 8;
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [0.5 0.5 wIn hIn]);
ax = axes(fig, 'Position', [0.02 0.02 0.96 0.86]);
hold(ax, 'on'); axis(ax, 'off');

t = linspace(0, 1, 60)';
bez = @(p0, p1, p2, p3) (1-t).^3*p0 + 3*(1-t).^2.*t*p1 + 3*(1-t).*t.^2*p2 + t.^3*p3;

% ---- rubans
for j = 1:nd-1
    F = zeros(8, 8);
    for a = 1:8
        for b = 1:8
            F(a, b) = sum(M(:, j) == a & M(:, j+1) == b);
        end
    end
    outCur = top(:, j);
    inCur  = top(:, j+1);
    inStart = nan(8, 8);
    for b = ORDER                     % entrées triées par ordre de la source
        for a = ORDER
            if F(a, b) > 0
                inStart(a, b) = inCur(b);
                inCur(b) = inCur(b) - F(a, b) * scale;
            end
        end
    end
    x0 = xs(j) + w/2; x1 = xs(j+1) - w/2; xm = (x0 + x1) / 2;
    for a = ORDER                     % sorties triées par ordre de la cible
        for b = ORDER
            if F(a, b) == 0, continue; end
            h  = F(a, b) * scale;
            y0 = outCur(a); outCur(a) = outCur(a) - h;
            y1 = inStart(a, b);
            col = RGB(a, :);
            if a == 7, col = [0.73 0.73 0.73]; end   % rubans "white" en gris
            xT = bez(x0, xm, xm, x1);
            yT = bez(y0, y0, y1, y1);
            patch(ax, [xT; flipud(xT)], [yT; flipud(yT) - h], col, ...
                'FaceAlpha', alpha, 'EdgeColor', 'none');
        end
    end
end

% ---- nœuds et étiquettes
for j = 1:nd
    last = (j == nd);
    for c = ORDER
        if cnt(c, j) == 0, continue; end
        rectangle(ax, 'Position', [xs(j) - w/2, bot(c, j), w, top(c, j) - bot(c, j)], ...
            'FaceColor', RGB(c, :), 'EdgeColor', [0.33 0.33 0.33], 'LineWidth', 0.6);
        if last
            xt = xs(j) - w/2 - 0.006; ha = 'right';
        else
            xt = xs(j) + w/2 + 0.006; ha = 'left';
        end
        text(ax, xt, (top(c, j) + bot(c, j)) / 2, sprintf('%s  n=%d', NAMES{c}, cnt(c, j)), ...
            'HorizontalAlignment', ha, 'VerticalAlignment', 'middle', 'FontSize', 9, ...
            'BackgroundColor', 'w', 'Margin', 0.5, 'Interpreter', 'none');
    end
    text(ax, xs(j), 1.02, labels{j}, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 11, 'FontWeight', 'bold', ...
        'Interpreter', 'none');
end

xlim(ax, [-0.04 1.04]);
ylim(ax, [min(bot(:)) - 0.02, 1.07]);
text(ax, 0, 1, ttl, 'Units', 'normalized', 'VerticalAlignment', 'bottom', ...
    'FontSize', 13, 'Interpreter', 'none');

save_fig(fig, [out '_sankey'], wIn, hIn, {'png', 'svg'});
end


% ========================================================================
function plot_paths(M, cellIds, labels, out, RGB, ORDER)
[nC, nd] = size(M);
[~, rk] = ismember(M, ORDER);
[~, ord] = sortrows(rk);              % tri : 1re date, puis 2e, ...
Ms = M(ord, :);

wIn = 1.4 * nd + 1.8; hIn = max(6, 0.09 * nC + 1.5);
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [0.5 0.5 wIn min(hIn, 10)]);
ax = axes(fig);
image(ax, Ms);                        % valeur k -> ligne k de la colormap
colormap(ax, RGB);
set(ax, 'XTick', 1:nd, 'XTickLabel', labels, 'YTick', 1:nC, ...
    'YTickLabel', compose('%d', cellIds(ord)), 'TickLength', [0 0], ...
    'TickLabelInterpreter', 'none');
ax.YAxis.FontSize = max(3, min(8, 600 / nC));
ax.XAxis.FontSize = 10;
xlabel(ax, 'Color assignment');
ylabel(ax, 'Cell ID');
title(ax, {'Color path per cell', sprintf('Cells present at all %d dates (n=%d)', nd, nC)}, ...
    'FontWeight', 'normal', 'FontSize', 10);

save_fig(fig, [out '_paths'], wIn, hIn, {'png'});
end


% ========================================================================
function save_fig(fig, base, wIn, hIn, fmts)
% taille d'impression indépendante de la taille de l'écran
set(fig, 'PaperUnits', 'inches', 'PaperPositionMode', 'manual', ...
    'PaperPosition', [0 0 wIn hIn], 'PaperSize', [wIn hIn]);
for k = 1:numel(fmts)
    switch fmts{k}
        case 'png'
            print(fig, [base '.png'], '-dpng', '-r200');
        case 'svg'
            set(fig, 'Renderer', 'painters');   % SVG vectoriel
            print(fig, [base '.svg'], '-dsvg');
    end
end
end


% ========================================================================
function write_tables(M, cellIds, labels, out, NAMES, ORDER)
[nC, nd] = size(M);

fid = fopen([out '_cells.csv'], 'w');
fprintf(fid, 'cell_id,%s\n', strjoin(labels, ','));
for i = 1:nC
    fprintf(fid, '%d,%s\n', cellIds(i), strjoin(NAMES(M(i, :)), ','));
end
fclose(fid);

fid = fopen([out '_transitions.csv'], 'w');
fprintf(fid, 'from_date,to_date,from_color,to_color,n\n');
for j = 1:nd-1
    for a = ORDER
        for b = ORDER
            k = sum(M(:, j) == a & M(:, j+1) == b);
            if k > 0
                fprintf(fid, '%s,%s,%s,%s,%d\n', labels{j}, labels{j+1}, NAMES{a}, NAMES{b}, k);
            end
        end
    end
end
fclose(fid);
end
