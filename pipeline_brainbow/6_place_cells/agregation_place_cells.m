function T = agregation_place_cells(dossier)
%AGREGATION_PLACE_CELLS  (18 nov. 2026) Regroupe les place_cells.mat de
%   <dossier>\<session>\ : une ligne par session, combinaison des p par
%   animal (Fisher) puis sur les animaux (Stouffer), figure de synthèse.
%   Écrit place_cells_par_session.csv, place_cells_par_animal.csv,
%   place_cells_biais_couleur.png dans <dossier>.

l = dir(fullfile(dossier, '*', 'place_cells.mat'));
lignes = {};
champs = {'NCell', 'n_tours', 'n_pc', 'pct_pc', 'n_colorees', 'n_pc_colorees', 'n_paires_meme', ...
    'p_proportion', 'dist_meme', 'dist_diff', 'p_distance', 'p_distance_blocs', ...
    'corr_meme', 'corr_diff', 'p_correlation', 'p_correlation_blocs', 'rho_anatomie_champ'};
for k = 1:numel(l)
    R = load(fullfile(l(k).folder, l(k).name));
    r = struct('session', R.nom, 'animal', R.animal, 'statut', R.statut);
    for c = champs
        if isfield(R, c{1}) && ~isempty(R.(c{1})), r.(c{1}) = double(R.(c{1})); else, r.(c{1}) = NaN; end
    end
    for c = 1:6
        if isfield(R, 'n_par_couleur'), r.(sprintf('pc_c%d', c)) = R.pc_par_couleur(c); r.(sprintf('n_c%d', c)) = R.n_par_couleur(c);
        else, r.(sprintf('pc_c%d', c)) = NaN; r.(sprintf('n_c%d', c)) = NaN; end
    end
    lignes{end+1} = r; %#ok<AGROW>
end
T = struct2table([lignes{:}]);
writetable(T, fullfile(dossier, 'place_cells_par_session.csv'));

tests = {'p_proportion', 'p_distance', 'p_distance_blocs', 'p_correlation', 'p_correlation_blocs'};
animaux = unique(T.animal);
A = table(animaux, 'VariableNames', {'animal'});
for t = tests
    pa = nan(numel(animaux), 1); na = zeros(numel(animaux), 1);
    for a = 1:numel(animaux)
        p = T.(t{1})(strcmp(T.animal, animaux{a})); p = p(~isnan(p));
        na(a) = numel(p);
        if ~isempty(p), pa(a) = 1 - chi2cdf(-2 * sum(log(p)), 2 * numel(p)); end   % Fisher
    end
    A.(t{1}) = pa; A.(['n_' t{1}]) = na;
end
A.n_pc = arrayfun(@(a) sum(T.n_pc(strcmp(T.animal, a)), 'omitnan'), animaux);
writetable(A, fullfile(dossier, 'place_cells_par_animal.csv'));

% Stouffer sur les animaux
glob = struct();
for t = tests
    p = A.(t{1}); p = p(~isnan(p));
    if isempty(p), glob.(t{1}) = NaN; else, glob.(t{1}) = 1 - normcdf(sum(norminv(1 - p)) / sqrt(numel(p))); end
end

fig = figure('Color', 'w', 'Position', [50 50 1500 450], 'Visible', 'off');
tl = tiledlayout(fig, 1, numel(tests) + 1, 'TileSpacing', 'compact');
ax = nexttile(tl); histogram(ax, T.pct_pc, 0:2:40); xlabel(ax, '% de place cells'); ylabel(ax, 'sessions');
title(ax, sprintf('%d sessions, %d place cells', sum(~isnan(T.n_pc)), sum(T.n_pc, 'omitnan')));
for k = 1:numel(tests)
    ax = nexttile(tl); p = T.(tests{k}); p = p(~isnan(p));
    histogram(ax, p, 0:0.1:1); yline(ax, numel(p) / 10, 'r--');
    xlabel(ax, 'p par session'); title(ax, {strrep(tests{k}, '_', ' '), sprintf('%d sessions ; animaux combinés p = %.3g', numel(p), glob.(tests{k}))});
end
title(tl, 'Place cells : biais de couleur (p uniformes = pas d''effet ; trait rouge = attendu)', 'FontWeight', 'bold');
exportgraphics(fig, fullfile(dossier, 'place_cells_biais_couleur.png'), 'Resolution', 150);
close(fig);
save(fullfile(dossier, 'place_cells_agregation.mat'), 'T', 'A', 'glob');
end
