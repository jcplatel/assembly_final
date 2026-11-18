function fig = figure_diagnostic_SCE(dossier)
%FIGURE_DIAGNOSTIC_SCE  Figure de diagnostic_SCE_assemblees, sans rien
%   recalculer : lit diagnostic_SCE_sessions.csv et diagnostic_SCE_assemblees.csv
%   dans dossier (8 nov.).
%     A. chaque session : SCE attendus par hasard / observés, et part des
%        activations d'assemblées que le hasard pourrait imiter ;
%     B. activations imitables selon la taille de l'assemblée ;
%     C. chaque assemblée : membres présents dans ses activations à petits
%        SCE, contre le 95e percentile des petits SCE du hasard ;
%     D. chaque assemblée : participation des membres selon la taille.
%   Écrit diagnostic_SCE.png dans dossier.
%   Ex. : figure_diagnostic_SCE('E:\...\resultats\diagnostic_SCE')

dossier = regexprep(char(dossier), '[\\/]+$', '');
Ts = readtable(fullfile(dossier, 'diagnostic_SCE_sessions.csv'));
Ta = readtable(fullfile(dossier, 'diagnostic_SCE_assemblees.csv'));
fig = figure('Color', 'w', 'Position', [40 60 1500 820]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Small SCE and assembly activations (%d sessions, %d assemblies)', height(Ts), height(Ta)), ...
    'FontWeight', 'bold', 'FontSize', 12);

% A. hasard / observés -> activations imitables
ax = nexttile(tl); hold(ax, 'on');
ratio = Ts.n_SCE_hasard ./ Ts.n_SCE;
imit = 100 * Ts.n_activations_imitables ./ max(Ts.n_activations, 1);
scatter(ax, 100 * ratio, imit, 40, [0.3 0.3 0.3], 'filled');
xline(ax, 50, ':', 'Color', [0.5 0.5 0.5]);
xline(ax, 100, '--', 'Color', [0.8 0.1 0.1]);
xlabel(ax, 'SCE expected by chance / observed (%)');
ylabel(ax, '% of activations chance could mimic');
title(ax, 'A. Each session', 'FontSize', 10);
box(ax, 'off'); set(ax, 'FontSize', 9, 'TickDir', 'out');

% B. selon la taille de l'assemblée
ax = nexttile(tl); hold(ax, 'on');
bornes = [0 10 15 20 30 Inf];
etiq = {'< 10', '10-14', '15-19', '20-29', '>= 30'};
v = nan(1, 5); nb = zeros(1, 5);
for k = 1:5
    m = Ta.n_membres >= bornes(k) & Ta.n_membres < bornes(k + 1);
    nb(k) = sum(m);
    v(k) = 100 * sum(Ta.activations_imitables(m)) / max(sum(Ta.n_activations(m)), 1);
end
bar(ax, 1:5, v, 0.6, 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', 'none');
for k = 1:5
    text(ax, k, v(k), sprintf('n = %d', nb(k)), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 8);
end
set(ax, 'XTick', 1:5, 'XTickLabel', etiq, 'FontSize', 9, 'TickDir', 'out');
xlabel(ax, 'assembly size (cells)');
ylabel(ax, '% of activations chance could mimic');
title(ax, 'B. By assembly size', 'FontSize', 10);
box(ax, 'off');

% C. membres dans les petites activations contre le hasard
ax = nexttile(tl); hold(ax, 'on');
m = isfinite(Ta.membres_dans_petites_activations) & isfinite(Ta.membres_petits_SCE_hasard_95e);
x = Ta.membres_petits_SCE_hasard_95e(m) + 0.15 * randn(sum(m), 1);
y = Ta.membres_dans_petites_activations(m) + 0.15 * randn(sum(m), 1);
scatter(ax, x, y, 25, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.6);
lim = [0, max([x; y; 1]) + 1];
plot(ax, lim, lim, '--', 'Color', [0.8 0.1 0.1]);
xlim(ax, lim); ylim(ax, lim);
xlabel(ax, 'members in small chance SCE (95th percentile)');
ylabel(ax, 'members in small activations (median)');
title(ax, 'C. Small SCE: real activations vs chance', 'FontSize', 10);
box(ax, 'off'); set(ax, 'FontSize', 9, 'TickDir', 'out');

% D. participation selon la taille
ax = nexttile(tl);
scatter(ax, Ta.n_membres, 100 * Ta.participation_mediane, 25, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.6);
set(ax, 'XScale', 'log', 'FontSize', 9, 'TickDir', 'out');
xlabel(ax, 'assembly size (cells, log)');
ylabel(ax, 'members present per activation (%)');
title(ax, 'D. Participation of members', 'FontSize', 10);
box(ax, 'off');
exportgraphics(fig, fullfile(dossier, 'diagnostic_SCE.png'), 'Resolution', 200);
end
