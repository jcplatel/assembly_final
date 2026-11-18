function lignes = compacite_assemblees(R, sessions, config, dossier_agreg)
%COMPACITE_ASSEMBLEES  Les assemblées sont-elles plus compactes dans l'espace
%   que des groupes de cellules tirés au hasard ? (étape 4, aucune couleur)
%
%   Pour chaque assemblée, l'étape 3 calcule déjà (stats_clonales_champ3) :
%     étalement = distance moyenne entre ses cellules (actives) ;
%     1000 groupes de même taille tirés au hasard parmi les cellules actives
%     du champ (au moins un SCE) -> rapport = étalement / étalement moyen du
%     hasard (< 1 : plus compacte) et compacite_pct = % des tirages aussi
%     compacts ou plus (p unilatéral x 100).
%   Ce n'est pas un test hypergéométrique (il compte des catégories, comme
%   les couleurs) : c'est un test par tirages au hasard (Monte-Carlo) sur une
%   distance.
%   Agrégation (l'animal est l'unité) : écart = rapport - 1 (0 sous
%   l'hypothèse nulle ; pas le log, qui est biaisé vers « compact » par
%   l'inégalité de Jensen), moyenné par session puis par animal ; Wilcoxon
%   unilatéral sur les animaux (plus compactes) ; LMM écart ~ 1 + (1|animal)
%   + (1|session), ddl de Satterthwaite.
%
%   R : table des assemblées de la variante principale (toutes, testables ou
%       non : la compacité ne dépend pas des couleurs)
%   Écrit <dossier_agreg>\controles_champ\compacite_assemblees.csv / .mat et
%   assembly_compactness.png (figure_compacite). lignes : texte du résumé.

lignes = {};
besoin = {'etalement_px', 'ratio_etalement', 'compacite_pct', 'animal', 'session', 'champ', 'condition'};
if ~all(ismember(besoin, R.Properties.VariableNames))
    lignes = {'Compacité des assemblées : colonnes absentes des résultats de l''étape 3 (la refaire).'};
    return
end
dossier = fullfile(dossier_agreg, 'controles_champ');
if ~isfolder(dossier), mkdir(dossier); end

% µm par pixel de chaque session (table des sessions, sinon valeur par défaut)
cle = strcat(R.animal, '_', R.session, '_', R.champ, '_', R.condition);
um = nan(height(R), 1);
if ismember('um_par_px', sessions.Properties.VariableNames)
    [ok_s, is] = ismember(cle, sessions.nom);
    um(ok_s) = sessions.um_par_px(is(ok_s));
end
um(isnan(um)) = config.um_par_px_defaut;

n_cel = nan(height(R), 1);
if ismember('n_total', R.Properties.VariableNames), n_cel = R.n_total; end
num = (1:height(R))';
if ismember('assemblee', R.Properties.VariableNames), num = R.assemblee; end
T = table(R.animal, cle, num, n_cel, R.etalement_px .* um, R.etalement_px .* um ./ R.ratio_etalement, ...
    R.ratio_etalement, R.compacite_pct / 100, 'VariableNames', ...
    {'animal', 'session', 'assemblee', 'n_cellules', 'etalement_um', 'etalement_hasard_um', 'rapport', 'p'});
T = T(isfinite(T.rapport) & T.rapport > 0, :);
if height(T) < 3
    lignes = {'Compacité des assemblées : trop peu d''assemblées.'};
    return
end
T.ecart = T.rapport - 1;
writetable(T, fullfile(dossier, 'compacite_assemblees.csv'));

% Session, puis animal
[gs, sess_u, anim_s] = findgroups(T.session, T.animal);
ecart_s = splitapply(@mean, T.ecart, gs);
[ga, animaux] = findgroups(anim_s);
ecart_a = splitapply(@mean, ecart_s, ga);
p_wil = NaN;
if numel(ecart_a) >= 2, p_wil = signrank(ecart_a, 0, 'tail', 'left'); end

% Modèle mixte
est = NaN; ic = [NaN NaN]; p_lmm = NaN;
try
    T.animal_c = categorical(T.animal);
    T.session_c = categorical(T.session);
    lme = fitlme(T, 'ecart ~ 1 + (1|animal_c) + (1|session_c)');
    [~, ~, c] = fixedEffects(lme, 'DFMethod', 'satterthwaite');
    est = c.Estimate(1); ic = [c.Lower(1), c.Upper(1)]; p_lmm = c.pValue(1);
catch err
    warning('compacite_assemblees : LMM impossible (%s)', err.message);
end
n_compactes = sum(T.p <= 0.05);

C = struct('T', T, 'ecart_a', ecart_a, 'animaux', {animaux}, 'ecart_s', ecart_s, 'sessions', {sess_u}, ...
    'p_wilcoxon', p_wil, 'rapport_lmm', 1 + est, 'rapport_ic', 1 + ic, 'p_lmm', p_lmm, ...
    'n_compactes', n_compactes); %#ok<NASGU>
save(fullfile(dossier, 'compacite_assemblees.mat'), '-struct', 'C');

% Texte du résumé
lignes{end+1} = sprintf('Compacité des assemblées (sans couleur ; %d assemblées, %d sessions, %d animaux) :', ...
    height(T), numel(sess_u), numel(animaux));
lignes{end+1} = sprintf(['   « Les cellules d''une assemblée sont en moyenne à %s µm les unes des autres, contre %s µm ' ...
    'pour des groupes de même taille tirés au hasard parmi les cellules actives : rapport %s ' ...
    '[IC 95 %% %s à %s] (LMM, animal et session aléatoires, ddl de Satterthwaite, p = %s) ; Wilcoxon ' ...
    'unilatéral sur les animaux (plus compactes) : p = %s, %d/%d animaux sous 1. »'], ...
    nb_c(mean(T.etalement_um), '%.0f'), nb_c(mean(T.etalement_hasard_um), '%.0f'), nb_c(1 + est, '%.2f'), ...
    nb_c(1 + ic(1), '%.2f'), nb_c(1 + ic(2), '%.2f'), pv_c(p_lmm), pv_c(p_wil), sum(ecart_a < 0), numel(ecart_a));
lignes{end+1} = sprintf(['   %d / %d assemblées plus compactes que 95 %% des tirages (%s attendues par hasard). ' ...
    'Test par tirages au hasard (1000 groupes de même taille par assemblée), pas hypergéométrique.'], ...
    n_compactes, height(T), nb_c(0.05 * height(T), '%.1f'));
fprintf('%s\n', lignes{:});

try
    figure_compacite(dossier);
catch err
    warning('figure_compacite non exécutée : %s', err.message);
end
end


function t = nb_c(x, fmt)
if isnan(x), t = 'NA'; else, t = strrep(sprintf(fmt, x), '.', ','); end
end


function t = pv_c(p)
if isnan(p)
    t = 'NA';
elseif p < 0.001
    t = '< 0,001';
else
    t = strrep(sprintf('%.3f', p), '.', ',');
end
end
