function fig = figure_sequence_session(dossier)
%FIGURE_SEQUENCE_SESSION  Séquence pendant la course d'UNE session (étape 5),
%   sans rien recalculer : lit <session>\sequences\sequences.mat.
%     A. moitié des courses qui range : activité moyenne autour du début de
%        course des cellules de la séquence, dans leur ordre ;
%     B. l'autre moitié des courses, MÊME ordre (une vraie séquence garde
%        sa diagonale) ;
%     C. décalage de chaque cellule de la séquence sur les deux moitiés
%        (rho de Spearman, unilatéral).
%   Écrit sequence_session.png dans le même dossier.

dossier = regexprep(char(dossier), '[\\/]+$', '');
R = load(fullfile(dossier, 'sequences.mat'));
fig = figure('Color', 'w', 'Position', [60 60 1500 520]);
if ~isfield(R, 'ordre_app')
    annotation(fig, 'textbox', [0.1 0.4 0.8 0.2], 'String', {R.nom, R.statut}, 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'Interpreter', 'none', 'FontSize', 12);
    exportgraphics(fig, fullfile(dossier, 'sequence_session.png'), 'Resolution', 150);
    return
end
fs = R.sampling_rate;
t = (-R.opts.psth_avant:R.opts.psth_apres) / fs;
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s: %d runs, %.0f/%d cells in the sequence', R.nom, R.n_courses, R.n_sequence, R.NCell), ...
    'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 11);
ordre = R.ordre_app;
panneaux = {R.psth_app, R.n_psth_app, 'A. Half of the runs (sorting)'; ...
    R.psth_test, R.n_psth_test, 'B. Other half, same order'};
for k = 1:2
    ax = nexttile(tl);
    M = double(panneaux{k, 1});
    if isempty(ordre) || all(isnan(M(:)))
        text(ax, 0.5, 0.5, 'No run far enough from the edges', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center');
        axis(ax, 'off');
        continue
    end
    M = smoothdata(M, 2, 'gaussian', 10);
    imagesc(ax, t, 1:numel(ordre), M(ordre, :));
    colormap(ax, 'hot');
    caxis(ax, [0.1 2]);
    xline(ax, 0, 'w', 'LineWidth', 2);
    xlabel(ax, 'time from run onset (s)');
    ylabel(ax, 'sequence cells (sorted on the first half)');
    title(ax, {panneaux{k, 3}, sprintf('%d runs', panneaux{k, 2})}, 'FontSize', 10);
    set(ax, 'TickDir', 'out', 'FontSize', 9);
end
ax = nexttile(tl); hold(ax, 'on');
seq = R.sequence_app;
x = -R.decalage_app(seq) / fs;
y = -R.decalage_test(seq) / fs;
scatter(ax, x, y, 18, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.6);
lim = [min([x, y, 0]), max([x, y, 1])];
plot(ax, lim, lim, '--', 'Color', [0.6 0.6 0.6]);
set(ax, 'XLim', lim, 'YLim', lim, 'TickDir', 'out', 'FontSize', 9);
axis(ax, 'square');
xlabel(ax, 'delay after run onset, first half (s)');
ylabel(ax, 'delay after run onset, other half (s)');
if isfield(R, 'kfold')
    % 19 nov. : the scatter is one 50/50 split; the session result is the cross-validation
    title(ax, {'C. Same order on the other half? (one split)', ...
        sprintf('this split: rho = %.2f, %s', R.rho_moitie, texte_p(R.p_rho_moitie)), ...
        sprintf('%s: mean rho = %.2f, %s', R.validation, R.rho, texte_p(R.p_rho))}, 'FontSize', 10);
else
    title(ax, {'C. Same order on the other half?', sprintf('rho = %.2f, %s', R.rho, texte_p(R.p_rho))}, ...
        'FontSize', 10);
end
box(ax, 'off');
exportgraphics(fig, fullfile(dossier, 'sequence_session.png'), 'Resolution', 150);
end


function t = texte_p(p)
if isnan(p)
    t = 'p = NA';
elseif p < 0.001
    t = 'p < 0.001';
else
    t = sprintf('p = %.3f', p);
end
end
