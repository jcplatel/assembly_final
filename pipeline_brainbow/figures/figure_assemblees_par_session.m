function fig = figure_assemblees_par_session(dossier_variante, fdr, alpha, nom_variante)
%FIGURE_ASSEMBLEES_PAR_SESSION  Toutes les sessions de tous les animaux :
%   nombre d'assemblées trouvées, testables, biaisées, et pour quelle couleur.
%
%   figure_assemblees_par_session('E:\...\agregation\champ_globale')
%   (lancée aussi par l'étape 4 pour chaque variante)
%
%   Lit resultats_biais_clonal.mat (table R de l'étape 4) dans
%   dossier_variante, et statut_sessions*.csv dans le dossier au-dessus
%   (sessions sans assemblée). Écrit assemblies_per_session.png.
%
%   A. une barre par session, rangées par animal :
%        couleur pleine   biaisée (q FDR < fdr), couleur de la couleur en excès
%        couleur claire   p_minp < alpha mais pas après correction FDR
%        gris             testable, sans biais
%        gris clair       non testable (trop peu de cellules colorées, ou
%                         critères de taille / d'activité de l'étape 4)
%      à droite : K, nombre d'assemblées, testables, biaisées (couleurs)
%   B. total par couleur (biaisées, puis p < alpha seulement)
%   Aucun test ici : seulement les résultats de l'étape 4, rangés.

if nargin < 2 || isempty(fdr), fdr = 0.10; end
if nargin < 3 || isempty(alpha), alpha = 0.05; end
if nargin < 4, nom_variante = ''; end
dossier_variante = regexprep(char(dossier_variante), '[\\/]+$', '');
S = load(fullfile(dossier_variante, 'resultats_biais_clonal.mat'), 'R');
R = S.R;
if isempty(nom_variante)
    [~, nom_variante] = fileparts(dossier_variante);
end

noms_coul = {'R', 'G', 'B', 'RG', 'RB', 'GB'};
cols = [0.86 0.16 0.16; 0.16 0.66 0.22; 0.20 0.38 0.92; ...
        0.93 0.78 0.10; 0.82 0.22 0.80; 0.10 0.74 0.80];
clair = 0.45 * cols + 0.55;
GRIS = [0.55 0.55 0.55]; GRIS_C = [0.85 0.85 0.85];

% Clé de session = nom des dossiers (animal_session_champ_condition)
R.cle = strcat(R.animal, '_', R.session, '_', R.champ, '_', R.condition);
biaisee = R.testable & R.q < fdr;
p05     = R.testable & R.p_minp < alpha & ~biaisee;
[~, ic_coul] = ismember(R.couleur_dominante, noms_coul);

% Sessions : celles de R, plus celles sans assemblée (statut_sessions*.csv)
cles = unique(R.cle, 'stable');
statut = repmat({''}, numel(cles), 1);
L = dir(fullfile(fileparts(dossier_variante), 'statut_sessions*.csv'));
if ~isempty(L)
    [~, i] = max([L.datenum]);
    try
        St = readtable(fullfile(L(i).folder, L(i).name), 'TextType', 'char', 'Delimiter', ',');
        autres = setdiff(St.session, cles, 'stable');
        cles = [cles; autres(:)];
        statut = [statut; St.statut(ismember(St.session, autres))];
    catch err
        warning('statut des sessions non lu : %s', err.message);
    end
end
% animal / session / champ depuis la clé, tri par animal puis session
morceaux = regexp(cles, '_', 'split');
animal = cellfun(@(m) m{1}, morceaux, 'UniformOutput', false);
etiquette = cellfun(@(m) strjoin(m(2:min(3, end)), ' '), morceaux, 'UniformOutput', false);
[~, ordre] = sort(strcat(animal, '|', etiquette));
cles = cles(ordre); statut = statut(ordre); animal = animal(ordre); etiquette = etiquette(ordre);
nS = numel(cles);

% Segments de chaque barre : 6 biaisées, 6 p < alpha, testables, non testables
M = zeros(nS, 14);
texte = cell(nS, 1);
for s = 1:nS
    m = strcmp(R.cle, cles{s});
    for c = 1:6
        M(s, c)     = sum(m & biaisee & ic_coul == c);
        M(s, 6 + c) = sum(m & p05 & ic_coul == c);
    end
    M(s, 13) = sum(m & R.testable) - sum(M(s, 1:12));
    M(s, 14) = sum(m & ~R.testable);
    if ~any(m)
        st = statut{s};
        if isempty(st), st = 'no assembly'; end
        texte{s} = strrep(st, 'exclue : ', 'excluded: ');
        continue
    end
    K = R.K(find(m, 1));
    if K == 0, tK = 'consensus'; else, tK = sprintf('K = %d', K); end
    coul_b = noms_coul(ic_coul(m & biaisee & ic_coul > 0));
    tb = sprintf('%d biased', numel(coul_b));
    if ~isempty(coul_b), tb = sprintf('%s (%s)', tb, strjoin(coul_b, ', ')); end
    texte{s} = sprintf('%s | %d asm, %d testable, %s', tK, sum(m), sum(m & R.testable), tb);
end

% Étiquettes : l'animal sur la première session de chaque groupe
nouvel_animal = [true; ~strcmp(animal(2:end), animal(1:end-1))];
lab = etiquette;
lab(nouvel_animal) = strcat(animal(nouvel_animal), {'  |  '}, etiquette(nouvel_animal));

%% Figure
hauteur = max(500, 140 + 19 * nS);
fig = figure('Color', 'w', 'Position', [40 30 1500 min(hauteur, 1900)]);
tl = tiledlayout(fig, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
n_b = sum(biaisee); n_t = sum(R.testable);
title(tl, {sprintf('Assemblies per session: %d sessions, %d animals, %d assemblies (%d testable), %d biased (FDR %g%%)', ...
    nS, numel(unique(animal)), height(R), n_t, n_b, 100 * fdr), texte_variante(nom_variante)}, ...
    'FontWeight', 'bold', 'FontSize', 12, 'Interpreter', 'none');

% A. une barre par session
ax = nexttile(tl, [1 3]); hold(ax, 'on');
y = (nS:-1:1)';
h = barh(ax, (1:nS)', flipud(M), 'stacked', 'EdgeColor', 'none', 'BarWidth', 0.75);
for c = 1:6
    h(c).FaceColor = cols(c, :);
    h(6 + c).FaceColor = clair(c, :);
end
h(13).FaceColor = GRIS; h(14).FaceColor = GRIS_C;
xmax = max(sum(M, 2));
if isempty(xmax) || xmax == 0, xmax = 1; end
for s = 1:nS
    if sum(M(s, :)) == 0, col_t = [0.6 0.2 0.2]; else, col_t = [0.25 0.25 0.25]; end
    text(ax, sum(M(s, :)) + 0.012 * xmax, y(s), texte{s}, 'FontSize', 8, 'Color', col_t, ...
        'VerticalAlignment', 'middle', 'Interpreter', 'none');
end
debuts = find(nouvel_animal);
for k = 2:numel(debuts)
    yline(ax, y(debuts(k)) + 0.5, '-', 'Color', [0.75 0.75 0.75]);
end
set(ax, 'YTick', 1:nS, 'YTickLabel', flipud(lab), 'YLim', [0.4 nS + 0.6], ...
    'XLim', [0 xmax * 1.75], 'TickLabelInterpreter', 'none', 'FontSize', 8, 'TickDir', 'out');
ax.XGrid = 'on'; box(ax, 'off');
xlabel(ax, 'number of assemblies');
title(ax, 'A. Assemblies of each session (grouped by animal)', 'FontWeight', 'bold', 'FontSize', 10);
leg = [h(1:6), h(7), h(13), h(14)];
noms_leg = [strcat({'biased, '}, noms_coul), {sprintf('p < %g only (light colour)', alpha), ...
    'testable, not biased', 'not testable'}];
legend(ax, leg, noms_leg, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
    'NumColumns', 5, 'Box', 'off', 'FontSize', 8);

% B. total par couleur
ax2 = nexttile(tl); hold(ax2, 'on');
tot_b = sum(M(:, 1:6), 1);
tot_p = sum(M(:, 7:12), 1);
% barre de bas en haut : GB ... R (R en haut)
hb = barh(ax2, (1:6)', [fliplr(tot_b); fliplr(tot_p)]', 'stacked', 'EdgeColor', 'none', 'BarWidth', 0.7);
hb(1).FaceColor = 'flat'; hb(2).FaceColor = 'flat';
hb(1).CData = flipud(cols); hb(2).CData = flipud(clair);
for c = 1:6
    text(ax2, tot_b(c) + tot_p(c) + 0.15, 7 - c, sprintf('%d + %d', tot_b(c), tot_p(c)), ...
        'FontSize', 8, 'VerticalAlignment', 'middle');
end
set(ax2, 'YTick', 1:6, 'YTickLabel', fliplr(noms_coul), 'YLim', [0.4 6.6], ...
    'XLim', [0 max(1, max(tot_b + tot_p)) * 1.4], 'FontSize', 9, 'TickDir', 'out');
ax2.XGrid = 'on'; box(ax2, 'off');
xlabel(ax2, 'assemblies');
title(ax2, {'B. Colour in excess', sprintf('(biased + p < %g only)', alpha)}, 'FontWeight', 'bold', 'FontSize', 10);

exportgraphics(fig, fullfile(dossier_variante, 'assemblies_per_session.png'), 'Resolution', 200);
end


function t = texte_variante(nom)
principal = contains(nom, 'PRINCIPAL');
base = strtrim(regexprep(nom, '\(.*\)', ''));
base = regexprep(base, '_PROVISOIRE$', '');
switch base
    case 'champ_globale',       t = 'field reference, colors shuffled across the field';
    case 'champ_blocs',         t = 'spatial control: colors shuffled only between nearby cells';
    case 'champ_territoire',    t = 'territory control: reference = colored cells within the assembly''s own area';
    case 'assemblees_blocs',    t = 'reference = colored assembly cells, shuffled within spatial blocks';
    otherwise,                  t = base;
end
if principal, t = ['MAIN TEST: ' t]; end
end
