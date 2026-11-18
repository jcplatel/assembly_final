function fig = figure_clonal_tous_K(dossier)
%FIGURE_CLONAL_TOUS_K  Biais clonal des assemblées pour chaque K
%   (main_clonal_tous_K), sans rien recalculer : lit clonal_tous_K.mat.
%     A. part des assemblées testables à p < 0,05 : p initial (non corrigé,
%        gris) et p_minp (corrigé, noir) ; pointillés = attendu par hasard ;
%     B. assemblées et assemblées testables ;
%     C. assemblées significatives après FDR (toutes les sessions, à ce K) ;
%     D. chaque animal : % à p_minp < 0,05 moins l'attendu.
%   K = 0 (consensus) est montré à gauche (« cons. »).
%   Écrit clonal_tous_K.png.

dossier = regexprep(char(dossier), '[\\/]+$', '');
S = load(fullfile(dossier, 'clonal_tous_K.mat'));
P = S.P; T = S.T;
Ks = P.K';
x = Ks;
x(Ks == 0) = min(Ks(Ks > 0)) - 2;                 % consensus à gauche
etiq = arrayfun(@num2str, Ks, 'UniformOutput', false);
etiq(Ks == 0) = {'cons.'};
fig = figure('Color', 'w', 'Position', [40 60 1500 760]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Clonal bias of assemblies at each K (exploratory; %s_%s)', S.variante{1}, S.variante{2}), ...
    'FontWeight', 'bold', 'FontSize', 12, 'Interpreter', 'none');
nt = max(P.n_testables', 1);

ax = nexttile(tl); hold(ax, 'on');
h1 = plot(ax, x, 100 * P.n_p_initial' ./ nt, '-o', 'Color', [0.6 0.6 0.6], 'MarkerFaceColor', [0.6 0.6 0.6], ...
    'LineWidth', 1.5);
plot(ax, x, 100 * P.attendu_p_initial' ./ nt, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
h2 = plot(ax, x, 100 * P.n_p_minp' ./ nt, '-o', 'Color', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 2);
plot(ax, x, 100 * P.attendu_p_minp' ./ nt, '--', 'Color', 'k', 'LineWidth', 1);
legend(ax, [h1 h2], {'initial p (uncorrected)', 'p_{minp} (corrected)'}, 'Location', 'best', 'Box', 'off');
reglages(ax, x, etiq);
ylabel(ax, sprintf('%% of testable assemblies at p < %g', S.seuil));
title(ax, 'A. Biased assemblies (dashed: chance)', 'FontSize', 10);

ax = nexttile(tl); hold(ax, 'on');
plot(ax, x, P.n_assemblees', '-o', 'Color', [0.6 0.6 0.6], 'MarkerFaceColor', [0.6 0.6 0.6], 'LineWidth', 1.5);
plot(ax, x, P.n_testables', '-o', 'Color', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 2);
legend(ax, {'assemblies', 'testable'}, 'Location', 'best', 'Box', 'off');
reglages(ax, x, etiq);
ylabel(ax, 'number (all sessions)');
title(ax, 'B. Assemblies', 'FontSize', 10);

ax = nexttile(tl); hold(ax, 'on');
bar(ax, x, P.n_fdr', 0.6, 'FaceColor', [0.80 0.10 0.10], 'EdgeColor', 'none');
for r = 1:height(P)
    text(ax, x(r), P.n_fdr(r), sprintf('%d/%d', P.n_animaux_exces(r), P.n_animaux(r)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
end
reglages(ax, x, etiq);
ylabel(ax, sprintf('assemblies significant (FDR %g%%)', 100 * S.q));
title(ax, 'C. After FDR (text: animals above chance)', 'FontSize', 10);

ax = nexttile(tl); hold(ax, 'on');
t = T(T.testable, :);
animaux = unique(t.animal);
Y = nan(numel(animaux), numel(Ks));
for i = 1:numel(animaux)
    for k = 1:numel(Ks)
        m = strcmp(t.animal, animaux{i}) & t.K == Ks(k);
        if any(m), Y(i, k) = 100 * (sum(t.p_minp(m) < S.seuil) - sum(t.taux_nul(m))) / sum(m); end
    end
    plot(ax, x, Y(i, :), '-', 'Color', [0.75 0.75 0.75]);
end
plot(ax, x, mean(Y, 1, 'omitnan'), '-o', 'Color', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 2);
yline(ax, 0, ':', 'Color', [0.5 0.5 0.5]);
reglages(ax, x, etiq);
ylabel(ax, '% at p_{minp} < 0.05 minus chance');
title(ax, 'D. Each animal (grey) and mean', 'FontSize', 10);
exportgraphics(fig, fullfile(dossier, 'clonal_tous_K.png'), 'Resolution', 200);
end


function reglages(ax, x, etiq)
set(ax, 'XTick', x, 'XTickLabel', etiq, 'XLim', [min(x) - 1, max(x) + 1], 'FontSize', 9, 'TickDir', 'out');
xlabel(ax, 'K');
box(ax, 'off');
end
