function T = figure_screening_comptes(dossiers, reference, dossier_sortie)
%FIGURE_SCREENING_COMPTES  Screening : récapitulatif en nombres de tous les
%   réglages du prétraitement essayés (14 nov. ; 15 nov. : meilleure
%   combinaison de chaque réglage, plusieurs screenings, détail par animal).
%
%   Une ligne par réglage du prétraitement, sur toutes les sessions du
%   screening. Figure d'ensemble (figure_screening_comptes.png) :
%     A. nombre de transients (pics) détectés
%     B. nombre de SCE détectés (trait : SCE attendus par hasard)
%     C. nombre d'assemblées (clair : toutes ; foncé : testables)
%     D. nombre d'assemblées biaisées (trait : attendu par hasard) et %
%     E. les réglages d'assemblées utilisés pour C et D sur cette ligne
%   Figure par animal (figure_screening_comptes_animaux.png) : SCE,
%   assemblées testables, assemblées biaisées et excès, une barre par animal.
%   Tableau : screening_comptes.csv.
%
%   C et D dépendent des réglages des assemblées (2e argument) :
%     'meilleure' (défaut) : pour chaque réglage du prétraitement, SA meilleure
%                 combinaison classée (excès le plus haut). Attention : c'est
%                 un maximum choisi parmi des milliers, donc gonflé par le tri.
%     'premiere'  : pour tous, les réglages d'assemblées de la combinaison
%                 classée première du screening.
%     une ligne du résumé, copiée telle quelle : la même pour tous, ex.
%       'clusters rang 95 | FDR_Q=0.05 | v2 | consensus 0.5 K4-14 | >= 15 cellules, >= 5 SCE'
%   Les points noirs donnent toujours la moyenne de toutes les combinaisons
%   d'assemblées de la ligne (pas de tri, donc pas gonflée).
%
%   figure_screening_comptes                     screening de config_brainbow
%   figure_screening_comptes(dossier)            dossier ...\screening\resultats
%   figure_screening_comptes({dossier1, dossier2})   plusieurs screenings, les uns sous les autres
%   figure_screening_comptes(dossiers, reference)
%   figure_screening_comptes(dossiers, reference, dossier_sortie)   (défaut : le 1er dossier)
%
%   Ne recalcule rien : lit screening.mat et, pour A et B, pretraitement.mat
%   de chaque étape 1 du screening (les dossiers doivent encore exister).
%   Assemblées biaisées par animal : retrouvées à partir de l'excès et du
%   nombre de testables par animal (vérifié contre le total ; « ~ » dans le
%   tableau si le compte ne tombe pas juste).

if nargin < 1 || isempty(dossiers)
    config = config_brainbow();
    dossiers = fullfile(fileparts(regexprep(char(config.PathSave), '[\\/]+$', '')), 'screening', 'resultats');
end
if nargin < 2 || isempty(reference), reference = 'meilleure'; end
dossiers = cellstr(dossiers);
if nargin < 3 || isempty(dossier_sortie), dossier_sortie = dossiers{1}; end
reference = char(reference);
nD = numel(dossiers);

%% Lecture des screenings : une ligne par (screening, prétraitement)
L = struct('groupe', {}, 'reglage', {}, 'combinaison', {}, 'n_transients', {}, 'n_SCE', {}, 'n_SCE_hasard', {}, ...
    'n_asm', {}, 'n_testables', {}, 'n_biaisees', {}, 'attendu', {}, 'exces_points', {}, 'p_selection', {}, ...
    'moy_testables', {}, 'moy_biaisees', {}, 'animaux', {}, 'SCE_animal', {}, 'transients_animal', {}, ...
    'testables_animal', {}, 'biaisees_animal', {}, 'exces_animal', {}, 'biaisees_animal_exact', {});
noms_groupes = cell(1, nD);
titres = cell(1, nD);
for d = 1:nD
    Sf = load(fullfile(dossiers{d}, 'screening.mat'), 'S');
    S = Sf.S; C = S.C; G = S.G;
    [~, noms_groupes{d}] = fileparts(fileparts(fullfile(dossiers{d}, 'x')));
    if strcmpi(noms_groupes{d}, 'resultats'), [~, noms_groupes{d}] = fileparts(fileparts(dossiers{d})); end
    nS = numel(S.sessions);
    nA = numel(S.animaux);
    titres{d} = sprintf('%s: %s, %d sessions, base %s', noms_groupes{d}, strjoin(S.animaux, ', '), nS, ...
        texte_base(S));
    % animal de chaque session
    animal_s = zeros(1, nS);
    for s = 1:nS
        for a = 1:nA
            if strncmp(S.sessions{s}, [S.animaux{a} '_'], numel(S.animaux{a}) + 1), animal_s(s) = a; end
        end
    end
    % réglages d'assemblées communs à tout le screening (sauf mode 'meilleure')
    cle = C.exces;
    cle(~(C.eligible > 0) | isnan(cle)) = -Inf;
    ref = [];
    if strcmpi(reference, 'premiere')
        [~, k] = max(cle);
        ref = [C.r(k), C.m(k), C.v(k), C.c(k), C.imc(k), C.ims(k)];
    elseif ~strcmpi(reference, 'meilleure')
        ref = lire_reference(reference, G);
    end
    iP = find(ismember({G.pre.nom}, S.pre_complets));
    for i = 1:numel(iP)
        p = iP(i);
        x = struct();
        x.groupe = d;
        if strcmp(G.pre(p).facteur, 'base')
            x.reglage = 'base';
        else
            x.reglage = sprintf('%s: %s', G.pre(p).facteur, G.pre(p).niveau);
        end
        % A, B : transients et SCE, lus dans l'étape 1 de chaque session
        nt = zeros(1, nA); ns = zeros(1, nA); nh = 0; lues = 0;
        for s = 1:nS
            f = fullfile(G.pre(p).dossier, S.sessions{s}, 'pretraitement.mat');
            if ~isfile(f) || animal_s(s) == 0, continue; end
            try
                P = load(f, 'n_transients_total', 'diagnostic');
                nt(animal_s(s)) = nt(animal_s(s)) + double(P.n_transients_total);
                ns(animal_s(s)) = ns(animal_s(s)) + double(P.diagnostic.NRace);
                if isfield(P.diagnostic, 'n_SCE_hasard') && ~isempty(P.diagnostic.n_SCE_hasard)
                    nh = nh + double(P.diagnostic.n_SCE_hasard);
                else
                    nh = NaN;
                end
                lues = lues + 1;
            catch
            end
        end
        if lues == nS
            x.n_transients = sum(nt); x.n_SCE = sum(ns); x.n_SCE_hasard = nh;
            x.SCE_animal = ns; x.transients_animal = nt;
        else
            x.n_transients = NaN; x.n_SCE = NaN; x.n_SCE_hasard = NaN;
            x.SCE_animal = nan(1, nA); x.transients_animal = nan(1, nA);
        end
        % C, D : la combinaison d'assemblées de cette ligne
        du_p = C.ip == p;
        if isempty(ref)
            cand = find(du_p & C.eligible > 0 & ~isnan(C.exces));
            k = [];
            if ~isempty(cand)
                [~, o] = sortrows([-C.exces(cand), -C.n_test(cand)]);
                k = cand(o(1));
            end
        else
            k = find(du_p & C.r == ref(1) & C.m == ref(2) & C.v == ref(3) & C.c == ref(4) & C.imc == ref(5) & ...
                C.ims == ref(6), 1);
        end
        x.animaux = S.animaux;
        if isempty(k)
            x.combinaison = '(no ranked combination)';
            x.n_asm = NaN; x.n_testables = NaN; x.n_biaisees = NaN; x.attendu = NaN; x.exces_points = NaN;
            x.p_selection = NaN;
            x.testables_animal = nan(1, nA); x.biaisees_animal = nan(1, nA); x.exces_animal = nan(1, nA);
            x.biaisees_animal_exact = false;
        else
            x.combinaison = sprintf('%s %g | %s | %s | %s | >=%g cells, >=%g SCE', G.clusters{C.r(k), 1}, ...
                G.clusters{C.r(k), 2}, G.membres(C.m(k)).nom, G.versions{C.v(k)}, G.choix_K(C.c(k)).nom, ...
                G.min_cellules(C.imc(k)), G.min_SCE(C.ims(k)));
            x.n_asm = C.n_asm(k); x.n_testables = C.n_test(k); x.n_biaisees = C.n_sig(k);
            x.attendu = C.attendu(k); x.exces_points = 100 * C.exces(k);
            x.p_selection = NaN;
            if isfield(C, 'p_selection'), x.p_selection = C.p_selection(k); end
            x.testables_animal = C.nt_animal(k, :);
            x.exces_animal = 100 * C.exces_animal(k, :);
            % biaisées par animal = testables x (excès + part attendue), arrondi
            taux = C.attendu(k) / max(C.n_test(k), 1);
            b = round(C.nt_animal(k, :) .* (C.exces_animal(k, :) + taux));
            b(~(C.nt_animal(k, :) > 0)) = 0;
            b(isnan(b)) = 0;
            x.biaisees_animal = b;
            x.biaisees_animal_exact = sum(b) == C.n_sig(k);
        end
        if any(du_p)
            x.moy_testables = mean(C.n_test(du_p)); x.moy_biaisees = mean(C.n_sig(du_p));
        else
            x.moy_testables = NaN; x.moy_biaisees = NaN;
        end
        L(end+1) = x; %#ok<AGROW>
    end
end
nL = numel(L);
if nL == 0, error('figure_screening_comptes : aucun réglage classé dans %s.', strjoin(dossiers, ', ')); end

% animaux : la liste de tous les screenings, dans l'ordre
animaux = {};
for i = 1:nL
    animaux = [animaux, setdiff(L(i).animaux, animaux, 'stable')]; %#ok<AGROW>
end
nA = numel(animaux);
SCE_a = nan(nL, nA); test_a = nan(nL, nA); biais_a = nan(nL, nA); exces_a = nan(nL, nA);
for i = 1:nL
    [~, ou] = ismember(L(i).animaux, animaux);
    SCE_a(i, ou) = L(i).SCE_animal; test_a(i, ou) = L(i).testables_animal;
    biais_a(i, ou) = L(i).biaisees_animal; exces_a(i, ou) = L(i).exces_animal;
end
noms = cell(nL, 1);
for i = 1:nL
    if nD > 1
        noms{i} = sprintf('%s | %s', noms_groupes{L(i).groupe}, L(i).reglage);
    else
        noms{i} = L(i).reglage;
    end
end
switch lower(reference)
    case 'meilleure', texte_mode = 'best ranked assembly settings of each row (selected among thousands: inflated)';
    case 'premiere',  texte_mode = 'assembly settings of the top-ranked combination of each screening';
    otherwise,        texte_mode = ['same assembly settings for all rows: ' reference];
end

T = struct();
T.reglage = noms;
T.combinaison = {L.combinaison}';
champs = {'n_transients', 'n_SCE', 'n_SCE_hasard', 'n_asm', 'n_testables', 'n_biaisees', 'attendu', ...
    'exces_points', 'p_selection', 'moy_testables', 'moy_biaisees'};
for i = 1:numel(champs), T.(champs{i}) = [L.(champs{i})]'; end
T.animaux = animaux;
T.SCE_animal = SCE_a; T.testables_animal = test_a; T.biaisees_animal = biais_a; T.exces_animal = exces_a;
T.biaisees_animal_exact = [L.biaisees_animal_exact]';
T.mode = texte_mode;

%% Tableau (écran et csv)
fprintf('\n=== SCREENING : COMPTES PAR RÉGLAGE DU PRÉTRAITEMENT ===\n');
for d = 1:nD, fprintf('%s\n', titres{d}); end
fprintf('Assemblées : %s\n', texte_mode);
fprintf('%-46s %10s %6s %8s %5s %5s %12s %7s   %s\n', 'réglage', 'transients', 'SCE', '(hasard)', 'asm', ...
    'test.', 'biaisées', 'excès', ['par animal (biaisées/testables) : ' strjoin(animaux, ', ')]);
for i = 1:nL
    par_animal = '';
    for a = 1:nA
        if isnan(test_a(i, a)), par_animal = [par_animal '    -  ']; continue; end %#ok<AGROW>
        marque = ' ';
        if ~L(i).biaisees_animal_exact, marque = '~'; end
        par_animal = [par_animal sprintf(' %s%2d/%-2d ', marque, biais_a(i, a), test_a(i, a))]; %#ok<AGROW>
    end
    fprintf('%-46s %10.0f %6.0f %8.0f %5.0f %5.0f %5.0f (%2.0f %%) %+7.1f   %s\n', noms{i}, T.n_transients(i), ...
        T.n_SCE(i), T.n_SCE_hasard(i), T.n_asm(i), T.n_testables(i), T.n_biaisees(i), ...
        100 * T.n_biaisees(i) / max(T.n_testables(i), 1), T.exces_points(i), par_animal);
end
if any(isnan(T.n_SCE))
    fprintf(2, ['Transients et SCE absents pour %d réglage(s) : pretraitement.mat introuvable ou illisible pour ' ...
        'certaines sessions (dossiers d''étape 1 du screening déplacés ou supprimés ?).\n'], sum(isnan(T.n_SCE)));
end
try
    if ~isfolder(dossier_sortie), mkdir(dossier_sortie); end
    fid = fopen(fullfile(dossier_sortie, 'screening_comptes.csv'), 'w', 'n', 'UTF-8');
    entete = [{'reglage', 'reglages_assemblees'}, champs];
    for a = 1:nA
        entete = [entete, {['SCE_' animaux{a}], ['testables_' animaux{a}], ['biaisees_' animaux{a}], ...
            ['exces_points_' animaux{a}]}]; %#ok<AGROW>
    end
    entete{end+1} = 'biaisees_par_animal_exact';
    fprintf(fid, '%s\n', strjoin(entete, ','));
    for i = 1:nL
        v = cellfun(@(c) sprintf('%.6g', T.(c)(i)), champs, 'UniformOutput', false);
        for a = 1:nA
            v = [v, {sprintf('%.6g', SCE_a(i, a)), sprintf('%.6g', test_a(i, a)), sprintf('%.6g', biais_a(i, a)), ...
                sprintf('%.6g', exces_a(i, a))}]; %#ok<AGROW>
        end
        fprintf(fid, '"%s","%s",%s,%d\n', noms{i}, T.combinaison{i}, strjoin(v, ','), L(i).biaisees_animal_exact);
    end
    fclose(fid);
    fprintf('Tableau : %s\n', fullfile(dossier_sortie, 'screening_comptes.csv'));
catch err
    fprintf(2, 'screening_comptes.csv non écrit (%s)\n', err.message);
end

%% Figures
gris = [0.55 0.55 0.55]; bleu = [0.35 0.50 0.75]; bleu_clair = [0.75 0.82 0.92]; rouge = [0.80 0.15 0.15];
y = (1:nL)';
groupe = [L.groupe]';
limites = find(diff(groupe) ~= 0) + 0.5;           % traits entre deux screenings
hauteur = max(460, 19 * nL + 230);
en_tete = [{'Settings screening, counts per preprocessing setting'}, titres, {['Assemblies: ' texte_mode]}];

% --- Figure 1 : ensemble
try
    fig = figure('Color', 'w', 'Position', [30 40 1850 hauteur]);
    tl = tiledlayout(fig, 1, 6, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, en_tete, 'FontSize', 9, 'Interpreter', 'none');

    ax = nexttile(tl); hold(ax, 'on');
    barh(ax, y, T.n_transients, 0.7, 'FaceColor', gris, 'EdgeColor', 'none');
    mise_en_forme(ax, noms, true, limites, 'A. Transients (peaks)', {'transients, all sessions', ' '});

    ax = nexttile(tl); hold(ax, 'on');
    barh(ax, y, T.n_SCE, 0.7, 'FaceColor', gris, 'EdgeColor', 'none');
    plot(ax, [T.n_SCE_hasard, T.n_SCE_hasard]', [y - 0.35, y + 0.35]', '-', 'Color', 'k', 'LineWidth', 1.5);
    mise_en_forme(ax, noms, false, limites, 'B. SCE detected', {'SCE, all sessions', '(black line: expected by chance)'});

    ax = nexttile(tl); hold(ax, 'on');
    barh(ax, y, T.n_asm, 0.7, 'FaceColor', bleu_clair, 'EdgeColor', 'none');
    barh(ax, y, T.n_testables, 0.7, 'FaceColor', bleu, 'EdgeColor', 'none');
    plot(ax, T.moy_testables, y, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 4);
    mise_en_forme(ax, noms, false, limites, 'C. Assemblies (light: all, dark: testable)', ...
        {'assemblies, all sessions', '(dot: testable, mean of all assembly settings)'});

    ax = nexttile(tl); hold(ax, 'on');
    barh(ax, y, T.n_biaisees, 0.7, 'FaceColor', rouge, 'EdgeColor', 'none');
    plot(ax, [T.attendu, T.attendu]', [y - 0.35, y + 0.35]', '-', 'Color', 'k', 'LineWidth', 1.5);
    plot(ax, T.moy_biaisees, y, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 4);
    xmax = max([T.n_biaisees; T.moy_biaisees; 1], [], 'omitnan');
    for i = 1:nL
        if isnan(T.n_testables(i)), continue; end
        text(ax, max(T.n_biaisees(i), 0) + 0.03 * xmax, y(i), sprintf('%d/%d (%.0f%%)', T.n_biaisees(i), ...
            T.n_testables(i), 100 * T.n_biaisees(i) / max(T.n_testables(i), 1)), 'FontSize', 7, ...
            'VerticalAlignment', 'middle');
    end
    xlim(ax, [0, 1.5 * xmax]);
    mise_en_forme(ax, noms, false, limites, 'D. Biased assemblies (p_minp < 0.05)', ...
        {'biased assemblies, all sessions', '(line: expected by chance; dot: mean of all assembly settings)'});

    ax = nexttile(tl, [1 2]); hold(ax, 'on');
    for i = 1:nL
        text(ax, 0.01, y(i), L(i).combinaison, 'FontSize', 7, 'VerticalAlignment', 'middle', 'Interpreter', 'none');
    end
    xlim(ax, [0 1]);
    mise_en_forme(ax, noms, false, limites, 'E. Assembly settings used on this row (C, D)', {' ', ' '});
    set(ax, 'XTick', [], 'XColor', 'none');

    exportgraphics(fig, fullfile(dossier_sortie, 'figure_screening_comptes.png'), 'Resolution', 150);
    fprintf('Figure : %s\n', fullfile(dossier_sortie, 'figure_screening_comptes.png'));
catch err
    warning('figure_screening_comptes (ensemble) : %s', err.message);
end

% --- Figure 2 : par animal
try
    couleurs_a = [0.00 0.45 0.70; 0.90 0.60 0.00; 0.00 0.62 0.45; 0.80 0.47 0.65; 0.34 0.71 0.91; 0.84 0.37 0.00];
    couleurs_a = couleurs_a(mod((1:nA) - 1, size(couleurs_a, 1)) + 1, :);
    fig = figure('Color', 'w', 'Position', [60 40 1750 hauteur + 40]);
    tl = tiledlayout(fig, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, [{'Settings screening, counts per animal'}, titres, {['Assemblies: ' texte_mode]}], 'FontSize', 9, ...
        'Interpreter', 'none');
    donnees = {SCE_a, test_a, biais_a, exces_a};
    titres_p = {'A. SCE detected', 'B. Testable assemblies', 'C. Biased assemblies (p_minp < 0.05)', ...
        'D. Excess of biased assemblies'};
    etiquettes = {'SCE per animal', 'testable assemblies per animal', 'biased assemblies per animal', ...
        'points above chance, per animal'};
    for q = 1:4
        ax = nexttile(tl); hold(ax, 'on');
        M = donnees{q};
        if nA == 1, M = [M, nan(nL, 1)]; end       % une seule série : barh demande une matrice
        b = barh(ax, y, M, 'grouped', 'EdgeColor', 'none');
        for a = 1:nA, b(a).FaceColor = couleurs_a(a, :); end
        if q == 4, xline(ax, 0, '-', 'Color', [0.3 0.3 0.3]); end
        mise_en_forme(ax, noms, q == 1, limites, titres_p{q}, etiquettes{q});
        legend(ax, b(1:nA), animaux, 'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off', ...
            'Interpreter', 'none');
    end
    exportgraphics(fig, fullfile(dossier_sortie, 'figure_screening_comptes_animaux.png'), 'Resolution', 150);
    fprintf('Figure : %s\n', fullfile(dossier_sortie, 'figure_screening_comptes_animaux.png'));
catch err
    warning('figure_screening_comptes (par animal) : %s', err.message);
end
end


function mise_en_forme(ax, noms, avec_noms, limites, titre, etiquette)
n = numel(noms);
set(ax, 'YDir', 'reverse', 'YLim', [0.3, n + 0.7], 'YTick', 1:n, 'TickDir', 'out', 'Box', 'off', ...
    'FontSize', 7, 'TickLabelInterpreter', 'none');
if avec_noms
    set(ax, 'YTickLabel', noms);
else
    set(ax, 'YTickLabel', {});
end
for k = 1:numel(limites)
    yline(ax, limites(k), '-', 'Color', [0.6 0.6 0.6]);
end
title(ax, titre, 'FontSize', 9, 'Interpreter', 'none');
xlabel(ax, etiquette, 'Interpreter', 'none');
end


function t = texte_base(S)
% Seuil des SCE de la base du screening, en clair
t = '?';
try
    pre = S.config_base.pre;
    switch pre.sce_regle
        case 'fixe',       t = sprintf('fixed %g cells', pre.sce_n_cells_threshold);
        case 'fdr',        t = sprintf('fdr %g', pre.sce_fdr);
        case 'percentile', t = sprintf('percentile %g', pre.sce_percentile);
        otherwise,         t = char(pre.sce_regle);
    end
    if isfield(pre, 'population') && ~isempty(pre.population), t = [t ', population ' char(pre.population)]; end
catch
end
end


function ref = lire_reference(texte, G)
% Ligne du résumé -> indices [clusters, membres, version, choix_K, min_cellules, min_SCE]
parts = strtrim(strsplit(texte, '|'));
parts = parts(~cellfun(@isempty, parts));
if ~isempty(parts) && ~isempty(regexp(parts{1}, '^pr.traitement', 'once'))
    parts(1) = [];
end
if numel(parts) ~= 5
    error(['Référence non comprise : "%s". Attendu : ''meilleure'', ''premiere'', ou une ligne du résumé comme ' ...
        '''clusters rang 95 | FDR_Q=0.05 | v2 | consensus 0.5 K4-14 | >= 15 cellules, >= 5 SCE'''], texte);
end
noms_cl = cellfun(@(a, b) sprintf('clusters %s %g', a, b), G.clusters(:, 1), G.clusters(:, 2), 'UniformOutput', false);
r = find(strcmp(noms_cl, parts{1}), 1);
m = find(strcmp({G.membres.nom}, parts{2}), 1);
v = find(strcmp(G.versions, parts{3}), 1);
c = find(strcmp({G.choix_K.nom}, parts{4}), 1);
nb = regexp(parts{5}, '>=\s*([\d.]+)\s*cellules\s*,\s*>=\s*([\d.]+)\s*SCE', 'tokens', 'once');
imc = []; ims = [];
if ~isempty(nb)
    imc = find(G.min_cellules == str2double(nb{1}), 1);
    ims = find(G.min_SCE == str2double(nb{2}), 1);
end
noms = {'clusters', 'membres', 'version', 'choix du K', 'minimum de cellules', 'minimum de SCE'};
trouve = {r, m, v, c, imc, ims};
manque = cellfun(@isempty, trouve);
if any(manque)
    error('Référence : %s introuvable(s) dans la grille du screening (texte : "%s").', ...
        strjoin(noms(manque), ', '), texte);
end
ref = [r, m, v, c, imc, ims];
end
