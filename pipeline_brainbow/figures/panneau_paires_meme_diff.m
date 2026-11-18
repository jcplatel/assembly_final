function panneau_paires_meme_diff(ax, S, mesure, permutation, pd, prefixe, type)
%PANNEAU_PAIRES_MEME_DIFF  Le panneau « même couleur contre couleurs
%   différentes » de l'étape 3b, dans les axes ax (utilisé par
%   figure_paires_resultats et figure_paires_meme_diff).
%   S : contenu de paires_agregation.mat ; mesure : 'sttc', 'dff', ... ;
%   permutation : 'couplage', ... ; pd : [dist_min dist_max] (disque : [0 R]) ;
%   prefixe : début du titre (ex. 'A. ').
%   type : 'residu' (défaut) ou 'brut' (5 nov. : mesure telle quelle, sans
%          retirer la distance ni le niveau de chaque cellule ; Wilcoxon sur
%          les animaux, pas de permutation).
%   Gris : chaque animal (moyenne de ses sessions) ; noir : moyenne des
%   animaux +/- sem. Valeurs = similarité corrigée de la distance et du
%   niveau de chaque cellule (résidus). Dans le titre : différence même -
%   différente couleur, Wilcoxon unilatéral sur les animaux et p de
%   permutation des sessions combinés.

if nargin < 6, prefixe = ''; end
if nargin < 7 || isempty(type), type = 'residu'; end
brut = strcmp(type, 'brut');
GRIS = [0.7 0.7 0.7];
if brut
    if ~isfield(S, 'AB') || isempty(S.AB)
        text(ax, 0.5, 0.5, {'No raw values yet:', 'run agregation\_paires again'}, 'Units', 'normalized', ...
            'HorizontalAlignment', 'center');
        axis(ax, 'off');
        return
    end
    selb = @(T) T(strcmp(T.mesure, mesure) & T.dist_min_um == pd(1) & T.dist_max_um == pd(2), :);
    T0 = selb(S.AniB);
    r0 = selb(S.AB);
    if ~isempty(T0)
        T0.moyenne_meme = T0.brut_meme;
        T0.moyenne_diff = T0.brut_diff;
    end
else
    Ani = S.Ani; A = S.A;
    sel = @(T) T(strcmp(T.mesure, mesure) & strcmp(T.permutation, permutation) & strcmp(T.couleur, 'toutes') & ...
        T.dist_min_um == pd(1) & T.dist_max_um == pd(2), :);
    T0 = sel(Ani);
    r0 = sel(A);
end
if S.cumulatif
    texte_d = sprintf('\\leq %g \\mum', pd(2));
else
    texte_d = sprintf('%g-%g \\mum', pd(1), pd(2));
end
texte_d = strrep(texte_d, '\leq Inf \mum', 'all pairs');
hold(ax, 'on');
if isempty(T0) || isempty(r0) || ~ismember('moyenne_meme', T0.Properties.VariableNames) || ...
        ~any(isfinite(T0.moyenne_meme))
    text(ax, 0.5, 0.5, {sprintf('No same / different means for %s, %s', mesure, texte_d), ...
        'run step 3b again (main\_paires)'}, 'Units', 'normalized', 'HorizontalAlignment', 'center');
    axis(ax, 'off');
    return
end
y = [T0.moyenne_meme, T0.moyenne_diff];
for a = 1:size(y, 1)
    plot(ax, [1 2], y(a, :), '-', 'Color', GRIS);
end
n_a = sum(all(isfinite(y), 2));
errorbar(ax, [1 2], mean(y, 1, 'omitnan'), std(y, 0, 1, 'omitnan') / sqrt(max(n_a, 1)), '-o', ...
    'Color', 'k', 'LineWidth', 2, 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'CapSize', 6);
set(ax, 'XTick', [1 2], 'XTickLabel', {'same colour', 'different colour'}, 'XLim', [0.6 2.4], ...
    'FontSize', 10, 'TickDir', 'out');
ax.YAxis.Exponent = 0;
if brut
    ylabel(ax, sprintf('%s (raw)', nom_mesure(mesure)));
    title(ax, {sprintf('%s%s raw, %s', prefixe, nom_mesure(mesure), texte_d), ...
        sprintf('%s, p = %s', texte_effet(r0.delta_pct(1), r0.delta_moyen(1)), texte_p(r0.p_wilcoxon(1)))}, ...
        'FontSize', 10);
else
    ylabel(ax, sprintf('%s (residual)', nom_mesure(mesure)));
    title(ax, {sprintf('%s%s residual, %s', prefixe, nom_mesure(mesure), texte_d), ...
        sprintf('%s, p = %s', texte_effet(r0.delta_pct(1), r0.delta_moyen(1)), texte_p(r0.p_wilcoxon(1)))}, ...
        'FontSize', 10);
end
box(ax, 'off');
end


function t = texte_effet(pct, delta)
% en % si la mesure a une échelle positive, sinon la différence brute (phi des SCE)
if isfinite(pct), t = sprintf('%+.1f%%', pct); else, t = sprintf('\\Delta = %+.4f', delta); end
end


function t = nom_mesure(m)
switch m
    case 'sttc', t = 'STTC';
    case 'dff',  t = 'Pearson dF/F';
    case 'sce',  t = 'SCE co-participation';
    otherwise
        t = strrep(strrep(m, 'sttc_', 'STTC '), 'ms', ' ms');
end
end


function t = texte_p(p)
if isnan(p)
    t = 'NA';
elseif p < 0.001
    t = '< 0.001';
else
    t = sprintf('%.3f', p);
end
end
