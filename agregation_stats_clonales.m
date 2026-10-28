%% =========================================================================
%  BIAIS CLONAL DES ASSEMBLÉES — regroupement de TOUS les champs et animaux
%  -------------------------------------------------------------------------
%  Lit tous les fichiers stats_clonales.mat produits par
%  stats_clonales_champ.m, puis :
%    A. TEST PRINCIPAL : la proportion d'assemblées biaisées dépasse-t-elle
%       les 5 % attendus par hasard ? (unité d'indépendance = l'animal)
%    B. Test secondaire : test global par champ (z_champ), moyenné par animal
%    C. CLASSEMENT de chaque assemblée : biaisée / non biaisée /
%       indéterminée / non testable (FDR de Benjamini-Hochberg)
%    D. CARACTÉRISATION : en quoi les assemblées biaisées diffèrent des
%       autres (séquentialité, condition), à taille égale
%    E. Robustesse : mesure continue -log10(p), autres seuils de FDR
%    F. Figure p_minp en fonction de la taille de l'assemblée
%
%  Pourquoi pas de correction entre assemblées pour le test A ?
%    Comparer la proportion d'assemblées à p < 5 % aux 5 % attendus tient
%    déjà compte du nombre d'assemblées : sur 40 assemblées sans effet, on
%    en attend ~2 significatives. Corriger chaque p avant de compter serait
%    une double correction. La correction (FDR) ne sert qu'au classement C.
%  =========================================================================


%% ---------- PARAMÈTRES (à fixer AVANT de regarder les résultats) ---------
dossier_racine   = 'E:\Data\Aurelie\analysis\Final\assembly\nocues\';  % dossier contenant les sessions
% Ne garder que les stats_clonales.mat du meilleur K de chaque session
% (liste produite par analyse_assembly_results5parfor_taille.m).
% Mettre '' pour prendre tous les fichiers trouvés.
fichier_best_k   = fullfile(dossier_racine, 'best_k_par_session.mat');
alpha            = 0.05;   % seuil par assemblée (p_minp) pour le test A
fdr_seuil        = 0.10;   % q < 0.10 -> assemblée « biaisée » (classement C)
seuil_non_biaise = 0.5;    % p_minp > 0.5 (et testable) -> « non biaisée »


%% ---------- Chargement --------------------------------------------------
fichiers = dir(fullfile(dossier_racine, '**', 'stats_clonales.mat'));
if isempty(fichiers)
    error('Aucun fichier stats_clonales.mat trouvé sous %s', dossier_racine);
end
% Filtre : uniquement les dossiers du meilleur K. Sinon, des analyses faites
% sur plusieurs K d'une même session seraient comptées plusieurs fois.
if ~isempty(fichier_best_k)
    bk = load(fichier_best_k, 'best_k');
    garde = ismember({fichiers.folder}, bk.best_k.dossier);
    fprintf('%d fichiers stats_clonales.mat trouvés, %d gardés (meilleur K de chaque session)\n', ...
        numel(fichiers), sum(garde));
    fichiers = fichiers(garde);
    if isempty(fichiers)
        error(['Aucun stats_clonales.mat dans les dossiers du meilleur K. ' ...
            'Vérifier que namefull pointait bien vers ces dossiers.']);
    end
end

R = table();
cles = cell(numel(fichiers), 1);
for f = 1:numel(fichiers)
    s = load(fullfile(fichiers(f).folder, fichiers(f).name), 'res');
    cles{f} = strjoin([s.res.animal(1), s.res.champ(1), s.res.session(1)], ' / ');
    R = [R; s.res]; %#ok<AGROW>
end
% Contrôle : une session ne doit apparaître qu'une fois
[cles_u, ~, ic] = unique(cles);
doublons = accumarray(ic, 1) > 1;
if any(doublons)
    error(['Ces sessions apparaissent dans plusieurs fichiers (plusieurs K, ' ...
        'ou identifiants non mis à jour) :\n%s'], strjoin(cles_u(doublons), '\n'));
end
% Identifiant unique du champ (un même champ est revu sur plusieurs jours)
R.champ_u    = strcat(R.animal, '_', R.champ);
R.session_u  = strcat(R.champ_u, '_', R.session);
R.animal_c   = categorical(R.animal);
R.champ_c    = categorical(R.champ_u);
R.condition_c = categorical(R.condition);
R.log_taille = log(max(R.n_colorees, 1));

fprintf('\n%d fichiers, %d animaux, %d champs, %d assemblées (%d testables)\n', ...
    numel(fichiers), numel(unique(R.animal)), numel(unique(R.champ_u)), ...
    height(R), sum(R.testable));

Rt = R(R.testable, :);          % seules les assemblées testables entrent dans les tests
Rt.sig = double(Rt.p_minp < alpha);


%% ---------- A. TEST PRINCIPAL : proportion d'assemblées biaisées --------
% A1. Une proportion par animal, puis test sur les animaux (simple, robuste)
[G, animaux] = findgroups(Rt.animal);
prop_animal  = splitapply(@mean,  Rt.sig, G);
n_animal     = splitapply(@numel, Rt.sig, G);
T_animal = table(animaux, n_animal, prop_animal, ...
    'VariableNames', {'animal', 'n_assemblees_testables', 'proportion_biaisees'});
disp(T_animal)
p_wilcoxon = signrank(prop_animal, alpha, 'tail', 'right');
fprintf(['A1. Proportion moyenne d''assemblées biaisées : %.1f %% (attendu par hasard : %.0f %%)\n' ...
         '    Wilcoxon sur les animaux (unilatéral) : p = %.4f  (n = %d animaux)\n'], ...
    100 * mean(prop_animal), 100 * alpha, p_wilcoxon, numel(prop_animal));
% NB : avec n animaux, la plus petite p possible du Wilcoxon est 1/2^n
% (ex. 0.03 avec 5 animaux). Avec peu d'animaux, A2 est plus puissant.

% A2. GLMM binomial avec offset : l'intercept teste « proportion = 5 % »
%     (l'offset place le zéro de l'intercept sur 5 % au lieu de 50 %)
try
    offs = repmat(log(alpha / (1 - alpha)), height(Rt), 1);
    glme0 = fitglme(Rt, 'sig ~ 1 + (1|animal_c) + (1|champ_c)', ...
        'Distribution', 'Binomial', 'Offset', offs);
    b0 = glme0.Coefficients.Estimate(1);
    p0 = glme0.Coefficients.pValue(1);          % bilatéral
    prop_est = 1 / (1 + exp(-(log(alpha / (1 - alpha)) + b0)));
    fprintf(['A2. GLMM : proportion estimée = %.1f %%  (p bilatéral = %.4g ; ' ...
             'unilatéral = %.4g si l''estimation dépasse 5 %%)\n'], ...
        100 * prop_est, p0, p0 / 2);
catch err
    fprintf('A2. GLMM impossible (%s). Garder A1.\n', err.message);
end


%% ---------- B. Test secondaire : test global par champ ------------------
% Une valeur par enregistrement, moyennée par animal
[Gs, ~] = findgroups(R.session_u);
z_sess   = splitapply(@(z) z(1), R.z_champ, Gs);
anim_ses = splitapply(@(x) x(1), R.animal, Gs);
[Ga, ~]  = findgroups(anim_ses);
z_animal = splitapply(@mean, z_sess, Ga);
fprintf('B.  z moyen du test par champ : %.2f ; Wilcoxon > 0 sur les animaux : p = %.4f\n', ...
    mean(z_animal), signrank(z_animal, 0, 'tail', 'right'));


%% ---------- C. CLASSEMENT de chaque assemblée ---------------------------
% FDR de Benjamini-Hochberg sur toutes les assemblées testables ensemble
R.q = nan(height(R), 1);
R.q(R.testable) = bh_qvalues(R.p_minp(R.testable));

R.categorie = repmat({'non_testable'}, height(R), 1);
R.categorie(R.testable & R.q <  fdr_seuil)                                  = {'biaisee'};
R.categorie(R.testable & R.q >= fdr_seuil & R.p_minp >  seuil_non_biaise)   = {'non_biaisee'};
R.categorie(R.testable & R.q >= fdr_seuil & R.p_minp <= seuil_non_biaise)   = {'indeterminee'};

fprintf('\nC.  Classement (FDR %.0f %%) :\n', 100 * fdr_seuil);
disp(groupcounts(R, 'categorie'))
disp(groupcounts(R, {'animal', 'categorie'}))


%% ---------- D. CARACTÉRISATION : biaisées contre non biaisées -----------
% Les « indéterminées » et « non testables » sont laissées de côté.
% La taille est TOUJOURS dans le modèle : la puissance du test augmente
% avec la taille, sinon on risque de « découvrir » un simple effet de taille.
% Ne jamais comparer les groupes sur la composition en couleurs (circulaire).
C = R(ismember(R.categorie, {'biaisee', 'non_biaisee'}), :);
C.biaisee = double(strcmp(C.categorie, 'biaisee'));

a_score = any(~isnan(C.score_sequence));
a_cond  = numel(unique(C.condition)) > 1;
termes  = {'log_taille'};
if a_score && a_cond
    termes{end+1} = 'score_sequence*condition_c';   % inclut l'interaction type x condition
elseif a_score
    termes{end+1} = 'score_sequence';
elseif a_cond
    termes{end+1} = 'condition_c';
end
formule_D = ['biaisee ~ 1 + ' strjoin(termes, ' + ') ' + (1|animal_c) + (1|champ_c)'];
try
    glme1 = fitglme(C, formule_D, 'Distribution', 'Binomial');
    fprintf('\nD.  %s\n', formule_D);
    disp(glme1.Coefficients)
catch err
    fprintf('D.  GLMM impossible (%s).\n', err.message);
end
% Les tests annoncés à l'avance (ex. score_sequence, condition, interaction)
% peuvent être corrigés entre eux par Holm. Le reste = exploratoire.


%% ---------- E. ROBUSTESSE ------------------------------------------------
% E1. Sans classement : -log10(p_minp) comme mesure continue du biais
Rt.mlogp = -log10(Rt.p_minp);
formule_E = ['mlogp ~ 1 + ' strjoin(termes, ' + ') ' + (1|animal_c) + (1|champ_c)'];
try
    lme2 = fitlme(Rt, formule_E);
    fprintf('\nE1. %s\n', formule_E);
    disp(lme2.Coefficients)
catch err
    fprintf('E1. LMM impossible (%s).\n', err.message);
end
% E2. Le nombre d'assemblées biaisées change-t-il beaucoup avec le seuil ?
for s = [0.05 0.10 0.20]
    fprintf('E2. FDR %.0f %% : %d assemblées biaisées\n', 100 * s, sum(R.q < s));
end


%% ---------- F. FIGURE : p_minp en fonction de la taille ------------------
% Gris : la meilleure p possible pour chaque assemblée. Une assemblée sous
% la ligne pointillée est biaisée ; une assemblée dont le point gris est
% au-dessus de la ligne ne pouvait pas l'être (non testable).
fig = figure('Color', 'w');
semilogy(R.n_colorees, R.p_min_atteignable, '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 10); hold on
cats = {'biaisee', 'indeterminee', 'non_biaisee', 'non_testable'};
cols = [0.85 0.1 0.1; 0.95 0.6 0.1; 0.2 0.4 0.8; 0.5 0.5 0.5];
for k = 1:numel(cats)
    ix = strcmp(R.categorie, cats{k});
    semilogy(R.n_colorees(ix), R.p_minp(ix), 'o', 'MarkerFaceColor', cols(k, :), ...
        'MarkerEdgeColor', 'none', 'MarkerSize', 6);
end
yline(alpha, '--k');
hold off
xlabel('Cellules colorées dans l''assemblée'); ylabel('p_{min-p}');
legend([{'p minimale atteignable'}, strrep(cats, '_', ' ')], 'Location', 'southwest');
title('Biais clonal des assemblées');
exportgraphics(fig, fullfile(dossier_racine, 'biais_clonal_p_vs_taille.png'), 'Resolution', 200);


%% ---------- Export -------------------------------------------------------
writetable(R(:, ~ismember(R.Properties.VariableNames, ...
    {'animal_c', 'champ_c', 'condition_c'})), ...
    fullfile(dossier_racine, 'resultats_biais_clonal.csv'));
writetable(T_animal, fullfile(dossier_racine, 'resultats_par_animal.csv'));
save(fullfile(dossier_racine, 'resultats_biais_clonal.mat'), 'R', 'T_animal');


%% ---------- Fonction locale ----------------------------------------------
function q = bh_qvalues(p)
% q-values de Benjamini-Hochberg
p = p(:);
m = numel(p);
[ps, ordre] = sort(p);
qs = ps .* m ./ (1:m)';
qs = min(1, cummin(qs, 'reverse'));
q = nan(m, 1);
q(ordre) = qs;
end
