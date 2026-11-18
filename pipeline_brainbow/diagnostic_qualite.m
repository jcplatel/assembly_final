function T = diagnostic_qualite(config)
%% =========================================================================
%  DIAGNOSTIC_QUALITE : qualité du signal de chaque session et de chaque
%  animal, SANS les couleurs (16 nov., demandé par l'utilisateur)
%  -------------------------------------------------------------------------
%  Sert à décider, sur l'activité seule, si un animal est à part (GCaMP
%  faible, peu d'assemblées) avant de le garder ou de le retirer.
%  Lit l'étape 1 déjà calculée (config.PathSave) : ne recalcule rien, ne lit
%  aucune couleur.
%
%  Lancer :   diagnostic_qualite            (config_brainbow.m)
%             diagnostic_qualite(config)
%
%  Par session :
%    n_cellules          cellules de la session
%    pct_actives         % de cellules avec au moins un transient au repos
%    repos_min           durée du repos (minutes)
%    transients_cell_min transients par cellule et par minute de repos (médiane des cellules actives)
%    bruit               bruit médian des cellules (unités de la trace)
%    amplitude_sur_bruit amplitude médiane des transients / bruit de la cellule (médiane des cellules)
%    n_SCE, SCE_hasard   SCE détectés, SCE attendus par hasard ; seuil_SCE
%    SCE_min             SCE par minute de repos
%    n_hyperactives      cellules hyperactives
%    asm_mediane, asm_max   nombre d'assemblées (version 1) : médiane et maximum sur les K
%  Par animal : médiane de ses sessions, et rapport à la médiane des autres
%  animaux (« x0.4 » = 40 % de la valeur des autres).
%  Écrit dans <PathSave>\diagnostic_qualite\ : qualite_sessions.csv,
%  qualite_animaux.csv, qualite.png.
%  =========================================================================

if nargin < 1 || isempty(config), config = config_brainbow(); end
addpath(genpath(fileparts(mfilename('fullpath'))));
config = dossiers_analyse(config);
sessions = lire_table_sessions(config);
n = height(sessions);

mesures = {'n_cellules', 'pct_actives', 'repos_min', 'transients_cell_min', 'bruit', 'amplitude_sur_bruit', ...
    'n_SCE', 'SCE_hasard', 'seuil_SCE', 'SCE_min', 'n_hyperactives', 'asm_mediane', 'asm_max'};
V = nan(n, numel(mesures));
fait = false(n, 1);
for s = 1:n
    d = char(sessions.dossier{s});
    if ismember('dossier_source', sessions.Properties.VariableNames) && ~isempty(sessions.dossier_source{s})
        d = char(sessions.dossier_source{s});
    end
    f = fullfile(d, 'pretraitement.mat');
    if ~isfile(f), continue; end
    try
        P = load(f, 'diagnostic', 'th_detection', 'PromALL', 'opts_utilisees', 'NCell');
        D = P.diagnostic;
        NCell = double(P.NCell);
        tpc = double(D.transients_percell(:));
        actives = tpc > 0;
        repos_min = D.duree_repos_s / 60;
        bruit = double(P.th_detection(:)) / P.opts_utilisees.threshold_peak;
        amp = nan(NCell, 1);
        for i = 1:NCell
            if ~isempty(P.PromALL{i}) && bruit(i) > 0
                amp(i) = median(double(P.PromALL{i})) / bruit(i);
            end
        end
        x = nan(1, numel(mesures));
        x(1) = NCell;
        x(2) = 100 * mean(actives);
        x(3) = repos_min;
        x(4) = median(tpc(actives)) / max(repos_min, eps);
        x(5) = median(bruit(bruit > 0), 'omitnan');
        x(6) = median(amp, 'omitnan');
        x(7) = D.NRace;
        x(8) = D.n_SCE_hasard;
        x(9) = D.seuil_SCE;
        x(10) = D.NRace / max(repos_min, eps);
        x(11) = numel(D.cellules_hyperactives_idx);
        % assemblées (version 1) à chaque K
        nb = [];
        for K = config.liste_K
            fk = fullfile(d, sprintf('k%d', K), 'assemblees.mat');
            if ~isfile(fk), continue; end
            A = load(fk, 'assemblyortho');
            nb(end+1) = numel(A.assemblyortho); %#ok<AGROW>
        end
        if ~isempty(nb), x(12) = median(nb); x(13) = max(nb); end
        V(s, :) = x;
        fait(s) = true;
    catch err
        fprintf(2, '%s : %s\n', sessions.nom{s}, err.message);
    end
end
if ~any(fait)
    error('diagnostic_qualite : aucune étape 1 trouvée dans %s.', config.PathEtape1);
end

%% Par session
T = struct();
T.session = sessions.nom;
T.animal = sessions.id_animal;
for k = 1:numel(mesures), T.(mesures{k}) = V(:, k); end
fprintf('\n=== QUALITÉ DU SIGNAL PAR SESSION (sans les couleurs) : %s ===\n', config.PathEtape1);
fprintf('%-42s %5s %6s %6s %9s %7s %8s %5s %7s %5s %7s %5s %5s\n', 'session', 'cell.', '%act.', 'repos', ...
    'tr/c/min', 'bruit', 'amp/br.', 'SCE', 'hasard', 'seuil', 'SCE/min', 'asm', 'max');
for s = 1:n
    fprintf('%-42s %5.0f %6.0f %6.1f %9.2f %7.3f %8.1f %5.0f %7.1f %5.0f %7.2f %5.1f %5.0f\n', sessions.nom{s}, ...
        V(s, 1), V(s, 2), V(s, 3), V(s, 4), V(s, 5), V(s, 6), V(s, 7), V(s, 8), V(s, 9), V(s, 10), V(s, 12), V(s, 13));
end

%% Par animal : médiane des sessions, et rapport aux autres animaux
[animaux, ~, ia] = unique(sessions.id_animal, 'stable');
nA = numel(animaux);
M = nan(nA, numel(mesures));
for a = 1:nA
    M(a, :) = median(V(ia == a & fait, :), 1, 'omitnan');
end
R = nan(nA, numel(mesures));
for a = 1:nA
    autres = setdiff(1:nA, a);
    R(a, :) = M(a, :) ./ median(M(autres, :), 1, 'omitnan');
end
a_montrer = [4 6 10 12];      % transients/cellule/min, amplitude/bruit, SCE/min, assemblées
fprintf('\n=== PAR ANIMAL (médiane de ses sessions ; entre parenthèses : rapport à la médiane des autres animaux) ===\n');
fprintf('%-10s %4s %22s %22s %22s %22s\n', 'animal', 'ses.', 'transients/cell./min', 'amplitude/bruit', ...
    'SCE/min', 'assemblées (méd. des K)');
for a = 1:nA
    fprintf('%-10s %4d', animaux{a}, sum(ia == a & fait));
    for k = a_montrer
        fprintf(' %12.2f (x%4.2f)  ', M(a, k), R(a, k));
    end
    fprintf('\n');
end
fprintf(['Lecture : un animal nettement sous les autres sur plusieurs mesures (par ex. x0.5 ou moins) est à part.\n' ...
    'Le seuil pour retirer un animal est à fixer AVANT de regarder ses couleurs.\n']);

%% Fichiers
dossier = fullfile(config.PathEtape1, 'diagnostic_qualite');
try
    if ~isfolder(dossier), mkdir(dossier); end
    Ts = array2table(V, 'VariableNames', mesures);
    Ts = [table(sessions.nom, sessions.id_animal, 'VariableNames', {'session', 'animal'}), Ts];
    writetable(Ts, fullfile(dossier, 'qualite_sessions.csv'));
    Ta = array2table([M, R], 'VariableNames', [mesures, strcat('rapport_', mesures)]);
    Ta = [table(animaux, 'VariableNames', {'animal'}), Ta];
    writetable(Ta, fullfile(dossier, 'qualite_animaux.csv'));
    fprintf('Tableaux : %s\n', dossier);
catch err
    fprintf(2, 'Tableaux non écrits (%s)\n', err.message);
end

%% Figure : une ligne par session, couleur par animal
try
    palette = [0.00 0.45 0.70; 0.90 0.60 0.00; 0.00 0.62 0.45; 0.80 0.47 0.65; 0.34 0.71 0.91; 0.84 0.37 0.00; ...
        0.95 0.90 0.25; 0.40 0.40 0.40; 0.60 0.30 0.10];
    coul = palette(mod(ia - 1, size(palette, 1)) + 1, :);
    titres = {'Transients per cell per min of rest', 'Transient amplitude / noise', 'SCE per min of rest', ...
        'SCE detected (line: expected by chance)', 'Assemblies (median over K)'};
    colonnes = [4 6 10 7 12];
    fig = figure('Color', 'w', 'Position', [40 60 1750 max(420, 20 * n + 170)]);
    tl = tiledlayout(fig, 1, numel(colonnes), 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, 'Signal quality per session (no colour information); one colour per animal', 'FontSize', 10);
    y = (1:n)';
    for q = 1:numel(colonnes)
        ax = nexttile(tl); hold(ax, 'on');
        b = barh(ax, y, V(:, colonnes(q)), 0.7, 'FaceColor', 'flat', 'EdgeColor', 'none');
        b.CData = coul;
        if colonnes(q) == 7
            plot(ax, [V(:, 8), V(:, 8)]', [y - 0.35, y + 0.35]', '-', 'Color', 'k', 'LineWidth', 1.5);
        end
        set(ax, 'YDir', 'reverse', 'YLim', [0.3, n + 0.7], 'YTick', 1:n, 'TickDir', 'out', 'Box', 'off', ...
            'FontSize', 7, 'TickLabelInterpreter', 'none');
        if q == 1
            set(ax, 'YTickLabel', sessions.nom);
        else
            set(ax, 'YTickLabel', {});
        end
        title(ax, titres{q}, 'FontSize', 9);
    end
    exportgraphics(fig, fullfile(dossier, 'qualite.png'), 'Resolution', 150);
    fprintf('Figure : %s\n', fullfile(dossier, 'qualite.png'));
catch err
    warning('diagnostic_qualite (figure) : %s', err.message);
end
end
