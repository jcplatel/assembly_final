function [Tsce, Tses] = etendue_spatiale_SCE(cfg, dossier_analyse, dossier_sortie, nTirages)
%ETENDUE_SPATIALE_SCE  (18 nov., demandée par l'utilisateur) Les SCE sont-ils
%   locaux (cellules actives regroupées dans le champ) ou étalés sur tout le
%   champ ? AUCUNE couleur.
%
%   Pour chaque SCE : distance moyenne entre les cellules actives (pixels),
%   comparée à nTirages tirages du même nombre de cellules parmi les
%   cellules qui participent à au moins un SCE, avec une probabilité
%   proportionnelle à leur nombre de SCE (garde le fait que certaines zones
%   du champ sont plus actives que d'autres).
%     rapport = distance observée / distance moyenne des tirages (< 1 : regroupé)
%     p       = part des tirages aussi regroupés ou plus (unilatéral)
%   Assemblées : K retenu à l'étape 2 (dossier_analyse\best_k_par_session.mat),
%   assemblees_v2.mat, mêmes critères que l'agrégation (>= 15 cellules, >= 5 SCE).
%   Écrit etendue_SCE_par_SCE.csv, etendue_SCE_par_session.csv,
%   etendue_spatiale_SCE.png, exemples_SCE_<session>.png (une session).

if nargin < 4, nTirages = 500; end
if ~isfolder(dossier_sortie), mkdir(dossier_sortie); end
cfg = dossiers_analyse(cfg);
S = lire_table_sessions(cfg);
B = load(fullfile(dossier_analyse, 'best_k_par_session.mat'), 'best_k'); B = B.best_k;
lignes = {}; ses = {};
for s = 1:height(S)
    nom = S.nom{s};
    fprintf('%d/%d %s\n', s, height(S), nom);
    try
        P = load(fullfile(S.dossier_source{s}, 'pretraitement.mat'), 'Race');
        Race = P.Race > 0;
        [~, ~, cen] = charger_geometrie([S.dossier_suite2p{s} filesep]);
        [nC, nS] = size(Race);
        if size(cen, 1) ~= nC, error('positions (%d) et Race (%d) : nombres de cellules différents', size(cen, 1), nC); end
        D = squareform(pdist(cen));
        % assemblées
        nAsm = zeros(1, nS); K = NaN;
        i = find(strcmp(B.session, nom), 1);
        if ~isempty(i) && ~isnan(B.K(i))
            K = B.K(i);
            A = load(fullfile(S.dossier_source{s}, sprintf('k%d', K), 'assemblees_v2.mat'), 'CellCl_final', 'PCl_final');
            C = A.CellCl_final > 0; Sa = A.PCl_final > 0;
            garde = sum(C, 2) >= 15 & sum(Sa, 2) >= 5;
            nAsm = sum(Sa(garde, :), 1);
        end
        part = sum(Race, 2); elig = find(part > 0); w = part(elig) / sum(part(elig));
        rng(s);
        rapport = nan(nS, 1); p = nan(nS, 1); dobs = nan(nS, 1); dnul = nan(nS, 1); n = sum(Race, 1)';
        for j = 1:nS
            a = find(Race(:, j)); k = numel(a);
            if k < 3, continue; end
            sub = D(a, a); dobs(j) = sum(sub(:)) / (k * (k - 1));
            nul = zeros(nTirages, 1);
            for t = 1:nTirages
                b = elig(tirage_sans_remise(w, k));
                sb = D(b, b); nul(t) = sum(sb(:)) / (k * (k - 1));
            end
            dnul(j) = mean(nul); rapport(j) = dobs(j) / dnul(j);
            p(j) = (1 + sum(nul <= dobs(j))) / (nTirages + 1);
        end
        Tj = table(repmat({nom}, nS, 1), (1:nS)', n, nAsm', dobs, dnul, rapport, p, ...
            'VariableNames', {'session', 'SCE', 'n_cellules', 'n_assemblees', 'dist_obs_px', 'dist_hasard_px', 'rapport', 'p'});
        lignes{end+1} = Tj; %#ok<AGROW>
        ok = ~isnan(p);
        r = struct('session', nom, 'animal', S.id_animal{s}, 'K', K, 'n_cellules', nC, 'n_SCE', nS, ...
            'pct_SCE_regroupes', 100 * mean(p(ok) < 0.05), ...
            'rapport_median', median(rapport(ok)), ...
            'pct_regroupes_sans_asm', 100 * mean(p(ok & nAsm' == 0) < 0.05), ...
            'pct_regroupes_avec_asm', 100 * mean(p(ok & nAsm' > 0) < 0.05), ...
            'rapport_sans_asm', median(rapport(ok & nAsm' == 0)), ...
            'rapport_avec_asm', median(rapport(ok & nAsm' > 0)), ...
            'rho_taille_rapport', corr(n(ok), rapport(ok), 'Type', 'Spearman'));
        ses{end+1} = r; %#ok<AGROW>
    catch err
        fprintf('  ERREUR : %s\n', err.message);
    end
end
Tsce = vertcat(lignes{:});
Tses = struct2table([ses{:}]);
writetable(Tsce, fullfile(dossier_sortie, 'etendue_SCE_par_SCE.csv'));
writetable(Tses, fullfile(dossier_sortie, 'etendue_SCE_par_session.csv'));

%% figure de synthèse
fig = figure('Color', 'w', 'Position', [50 50 1500 420]);
tl = tiledlayout(fig, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl); histogram(ax, Tsce.p, 0:0.05:1, 'Normalization', 'probability', 'FaceColor', [0.3 0.3 0.3]);
hold(ax, 'on'); yline(ax, 0.05, 'r--', 'hasard');
xlabel(ax, 'p (regroupement spatial)'); ylabel(ax, 'part des SCE');
title(ax, sprintf('A. Tous les SCE (n = %d) : %.0f %% à p < 0,05', sum(~isnan(Tsce.p)), 100 * mean(Tsce.p(~isnan(Tsce.p)) < 0.05)));
ax = nexttile(tl); hold(ax, 'on');
[~, o] = sort(Tses.pct_SCE_regroupes);
bar(ax, Tses.pct_SCE_regroupes(o), 'FaceColor', [0.4 0.55 0.8]); yline(ax, 5, 'r--', '5 % attendu');
set(ax, 'XTick', 1:height(Tses), 'XTickLabel', Tses.session(o), 'TickLabelInterpreter', 'none', 'FontSize', 6); xtickangle(ax, 70);
ylabel(ax, '% de SCE regroupés (p < 0,05)'); title(ax, 'B. Par session');
ax = nexttile(tl); hold(ax, 'on');
g = Tsce.n_assemblees > 0; ok = ~isnan(Tsce.rapport);
histogram(ax, Tsce.rapport(ok & ~g), 0.5:0.025:1.3, 'Normalization', 'probability', 'DisplayName', 'SCE sans assemblée');
histogram(ax, Tsce.rapport(ok & g), 0.5:0.025:1.3, 'Normalization', 'probability', 'DisplayName', 'SCE avec assemblée');
xline(ax, 1, 'k--'); legend(ax, 'Location', 'northwest', 'Box', 'off');
xlabel(ax, 'distance observée / hasard (< 1 : regroupé)'); ylabel(ax, 'part des SCE');
title(ax, sprintf('C. sans assemblée : %.0f %% regroupés ; avec : %.0f %%', ...
    100 * mean(Tsce.p(ok & ~g) < 0.05), 100 * mean(Tsce.p(ok & g) < 0.05)));
ax = nexttile(tl); hold(ax, 'on');
bords = [3 8 12 16 20 30 50 1000];
[~, ~, bin] = histcounts(Tsce.n_cellules, bords);
pc = arrayfun(@(b) 100 * mean(Tsce.p(ok & bin == b) < 0.05), 1:numel(bords) - 1);
nb = arrayfun(@(b) sum(ok & bin == b), 1:numel(bords) - 1);
bar(ax, pc, 'FaceColor', [0.55 0.4 0.7]); yline(ax, 5, 'r--');
lab = arrayfun(@(b) sprintf('%d-%d', bords(b), bords(b + 1) - 1), 1:numel(bords) - 1, 'UniformOutput', false);
lab{end} = sprintf('>=%d', bords(end - 1));
set(ax, 'XTick', 1:numel(pc), 'XTickLabel', lab, 'FontSize', 8);
for b = 1:numel(pc), text(ax, b, pc(b) + 0.5, sprintf('n=%d', nb(b)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7); end
xlabel(ax, 'cellules dans le SCE'); ylabel(ax, '% regroupés (p < 0,05)'); title(ax, 'D. Selon la taille du SCE');
title(tl, 'Étendue spatiale des SCE : cellules actives plus proches que le hasard (tirages pondérés par l''activité) ?', 'FontWeight', 'bold');
exportgraphics(fig, fullfile(dossier_sortie, 'etendue_spatiale_SCE.png'), 'Resolution', 170);
end


function idx = tirage_sans_remise(w, k)
% k indices tirés sans remise, probabilité proportionnelle à w
% (méthode d'Efraimidis-Spirakis : les k plus grandes clés u^(1/w))
cle = rand(numel(w), 1) .^ (1 ./ w);
[~, o] = maxk(cle, k); idx = o;
end
