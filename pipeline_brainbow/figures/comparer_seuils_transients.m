function T = comparer_seuils_transients(variantes, sessions, dossier_sortie)
%COMPARER_SEUILS_TRANSIENTS  (18 nov., demandée par l'utilisateur) Compare
%   l'activité et la qualité des assemblées entre plusieurs réglages de
%   l'étape 1 (ex. threshold_peak 3,09 / 3,5 / 4 / 5). AUCUNE couleur.
%
%   variantes : struct array, champs nom, etape1 (dossier de l'étape 1),
%               analyse (dossier des étapes 2 à 4 : best_k_par_session.mat)
%   sessions  : cellstr des noms de sessions
%
%   Pour chaque session et réglage, au K retenu à l'étape 2 :
%     transients / cellule / min de repos, SCE / min de repos, seuil des SCE
%     assemblées gardées (>= 15 cellules et >= 5 SCE, comme l'agrégation)
%     % de cellules dans une assemblée (Malvache : la plupart)
%     % de SCE à 0, 1, >= 2 assemblées (Malvache : 33 / 63 / 3,4 %)
%     contraste : participation des membres dans les SCE de leur assemblée
%       (dedans) contre dans les autres SCE (dehors), rapport dedans/dehors
%     pureté : part des cellules actives d'un SCE de l'assemblée qui en
%       sont membres
%   Écrit comparaison_seuils_transients.csv et un raster cellules x SCE
%   (toutes les cellules, triées par assemblée) par session :
%   rasters_<session>.png

if nargin < 3, dossier_sortie = pwd; end
min_cell = 15; min_sce = 5;
lignes = {};
for s = 1:numel(sessions)
    fig = figure('Color', 'w', 'Position', [40 40 380 * numel(variantes) 520], 'Visible', 'off');
    tl = tiledlayout(fig, 1, numel(variantes), 'TileSpacing', 'compact', 'Padding', 'compact');
    for v = 1:numel(variantes)
        V = variantes(v);
        r = struct('session', sessions{s}, 'reglage', V.nom);
        try
            B = load(fullfile(V.analyse, 'best_k_par_session.mat'), 'best_k');
            i = find(strcmp(B.best_k.session, sessions{s}), 1);
            K = B.best_k.K(i);
            d1 = fullfile(V.etape1, sessions{s});
            P = load(fullfile(d1, 'pretraitement.mat'), 'Raster', 'Race', 'WinRest', 'sampling_rate', 'diagnostic', 'NCell');
            A = load(fullfile(d1, sprintf('k%d', K), 'assemblees_v2.mat'), 'CellCl_final', 'PCl_final');
            fs = P.sampling_rate; rep = P.WinRest(:)'; min_rep = numel(rep) / fs / 60;
            Race = P.Race > 0; [nC, nS] = size(Race);
            r.K = K; r.n_cellules = nC; r.repos_min = min_rep;
            r.transients_cell_min = sum(sum(P.Raster(:, rep))) / nC / min_rep;
            r.seuil_SCE = P.diagnostic.seuil_SCE; r.n_SCE = nS; r.SCE_min = nS / min_rep;
            r.pct_SCE_hasard = 100 * P.diagnostic.proportion_SCE_hasard;
            C = A.CellCl_final > 0; S = A.PCl_final > 0;
            garde = sum(C, 2) >= min_cell & sum(S, 2) >= min_sce;
            C = C(garde, :); S = S(garde, :); nA = sum(garde);
            r.n_asm = nA;
            r.taille_mediane = median(sum(C, 2));
            r.pct_cellules_asm = 100 * mean(any(C, 1));
            nAsm = sum(S, 1);
            r.pct_SCE_0 = 100 * mean(nAsm == 0); r.pct_SCE_1 = 100 * mean(nAsm == 1); r.pct_SCE_2plus = 100 * mean(nAsm >= 2);
            dedans = nan(nA, 1); dehors = nan(nA, 1); purete = nan(nA, 1);
            for a = 1:nA
                m = C(a, :); e = S(a, :);
                dedans(a) = mean(mean(Race(m, e)));
                if any(~e), dehors(a) = mean(mean(Race(m, ~e))); end
                act = sum(Race(:, e), 1); purete(a) = mean(sum(Race(m, e), 1) ./ max(act, 1));
            end
            r.particip_dedans = 100 * median(dedans); r.particip_dehors = 100 * median(dehors);
            r.contraste = median(dedans ./ dehors); r.purete = 100 * median(purete);
            % raster : cellules triées par assemblée (hors assemblée en haut), SCE triés par assemblée
            [~, ac] = max([~any(C, 1); C], [], 1);            % 1 = hors assemblée
            [~, as] = max([nAsm == 0; S], [], 1);
            [~, oc] = sort(ac); [~, os] = sort(as);
            img = ones(nC, nS, 3);
            cols = lines(max(nA, 1));
            R = Race(oc, os); acs = ac(oc); ass = as(os);
            for c = 1:nC
                for j = find(R(c, :))
                    if acs(c) > 1 && acs(c) == ass(j), img(c, j, :) = cols(acs(c) - 1, :);
                    else, img(c, j, :) = [0.35 0.35 0.35];
                    end
                end
            end
            ax = nexttile(tl); image(ax, img); axis(ax, 'tight');
            title(ax, {sprintf('%s : K %d, %d assemblées', V.nom, K, nA), ...
                sprintf('%.1f SCE/min (seuil %d), %.0f %% des SCE à 1 assemblée', r.SCE_min, r.seuil_SCE, r.pct_SCE_1), ...
                sprintf('contraste %.1f, pureté %.0f %%, %.0f %% des cellules en assemblée', r.contraste, r.purete, r.pct_cellules_asm)}, ...
                'FontSize', 8, 'Interpreter', 'none');
            xlabel(ax, 'SCE (triés)'); ylabel(ax, 'cellules (triées)');
        catch err
            r.erreur = err.message;
            ax = nexttile(tl); axis(ax, 'off'); title(ax, [V.nom ' : ' err.message], 'Interpreter', 'none', 'FontSize', 8);
        end
        lignes{end+1} = r; %#ok<AGROW>
    end
    title(tl, sessions{s}, 'Interpreter', 'none');
    exportgraphics(fig, fullfile(dossier_sortie, ['rasters_' sessions{s} '.png']), 'Resolution', 150);
    close(fig);
end
noms = {};
for k = 1:numel(lignes), noms = union(noms, fieldnames(lignes{k}), 'stable'); end
for k = 1:numel(lignes)
    for f = setdiff(noms, fieldnames(lignes{k}))'
        if strcmp(f{1}, 'erreur'), lignes{k}.(f{1}) = ''; else, lignes{k}.(f{1}) = NaN; end
    end
    lignes{k} = orderfields(lignes{k}, noms);
end
T = struct2table([lignes{:}]);
writetable(T, fullfile(dossier_sortie, 'comparaison_seuils_transients.csv'));
end
