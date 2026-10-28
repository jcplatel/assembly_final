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
% Même enregistrement présent plusieurs fois (analysé deux fois) :
%   'plus_recent' : garde l'analyse la plus récente, retire les autres (avec un avertissement)
%   'erreur'      : s'arrête et liste les dossiers (à utiliser pour l'analyse définitive)
gestion_doublons = 'plus_recent';
figure_non_testables = false;   % true : montrer aussi les assemblées non testables sur la figure F
% Critères minimaux de validité d'une assemblée (définis sur l'ACTIVITÉ, jamais
% sur la couleur). Mettre 0 pour désactiver (analyse de sensibilité).
min_cellules_asm = 10;   % nombre minimal de cellules de l'assemblée (toutes, colorées ou non)
min_SCE_asm      = 20;   % nombre minimal d'activations significatives (n_SCE_signif, RACEORTHO)
fdr_seuil_figure     = 0.10;    % seuil FDR utilisé pour colorer la figure F (affichage seulement ;
                                % le classement C et la caractérisation D utilisent fdr_seuil)


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

tous_res = cell(numel(fichiers), 1);
cles = cell(numel(fichiers), 1);
for f = 1:numel(fichiers)
    s = load(fullfile(fichiers(f).folder, fichiers(f).name), 'res');
    tous_res{f} = s.res;
    % La condition fait partie de l'identifiant : un même plan enregistré le
    % même jour avec et sans indices = deux sessions distinctes
    cles{f} = strjoin([s.res.animal(1), s.res.champ(1), s.res.session(1), s.res.condition(1)], ' / ');
end

% Contrôle : une session ne doit apparaître qu'une fois
[cles_u, ~, ic] = unique(cles);
doublons = find(accumarray(ic, 1) > 1);
garde = true(numel(fichiers), 1);
if ~isempty(doublons) && strcmp(gestion_doublons, 'erreur')
    msg = '';
    for d = doublons'
        msg = [msg sprintf('\n%s, trouvée dans :\n   %s', cles_u{d}, ...
            strjoin({fichiers(ic == d).folder}, '\n   '))]; %#ok<AGROW>
    end
    error(['Ces sessions apparaissent dans plusieurs fichiers (même enregistrement ' ...
        'analysé deux fois, ou identifiants identiques) :%s'], msg);
end
for d = doublons'
    idx = find(ic == d);
    dossiers_d = {fichiers(idx).folder};
    % Le nom du dossier de session se termine par la date de l'analyse
    % (aa_mm_jj_hh_mm_ss) : l'ordre alphabétique est l'ordre chronologique.
    % On garde la plus récente. Règle indépendante des résultats clonaux.
    [~, ordre] = sort(dossiers_d);
    garde(idx(ordre(1:end-1))) = false;
    warning('Doublon %s : gardé\n   %s\n retiré :\n   %s', cles_u{d}, ...
        dossiers_d{ordre(end)}, strjoin(dossiers_d(ordre(1:end-1)), '\n   '));
end
fichiers = fichiers(garde);
tous_res = tous_res(garde);
% Les fichiers produits par une ancienne version n'ont pas toutes les
% colonnes : on complète les colonnes manquantes par NaN.
listes   = cellfun(@(t) t.Properties.VariableNames, tous_res, 'UniformOutput', false);
noms_var = unique([listes{:}], 'stable');
for f = 1:numel(tous_res)
    manq = setdiff(noms_var, tous_res{f}.Properties.VariableNames, 'stable');
    for v = 1:numel(manq)
        % même type que dans les fichiers qui ont cette colonne
        modele = tous_res{find(cellfun(@(t) ismember(manq{v}, t.Properties.VariableNames), tous_res), 1)};
        if iscell(modele.(manq{v}))
            tous_res{f}.(manq{v}) = repmat({''}, height(tous_res{f}), 1);
        else
            tous_res{f}.(manq{v}) = nan(height(tous_res{f}), size(modele.(manq{v}), 2));
        end
    end
    tous_res{f} = tous_res{f}(:, noms_var);
end
R = vertcat(tous_res{:});
% Identifiant unique du champ (un même champ est revu sur plusieurs jours)
R.champ_u    = strcat(R.animal, '_', R.champ);
R.session_u  = strcat(R.champ_u, '_', R.session, '_', R.condition);
R.animal_c   = categorical(R.animal);
R.champ_c    = categorical(R.champ_u);
R.condition_c = categorical(R.condition);
R.log_taille = log(max(R.n_colorees, 1));

% Critères de validité des assemblées : appliqués AVANT tous les tests et le
% FDR. Les assemblées exclues sont traitées comme non testables.
R.exclue_criteres = R.n_total < min_cellules_asm;
if ismember('n_SCE_signif', R.Properties.VariableNames)
    R.exclue_criteres = R.exclue_criteres | R.n_SCE_signif < min_SCE_asm;
elseif min_SCE_asm > 0
    warning('n_SCE_signif absent : le critère sur le nombre de SCE n''est pas appliqué.');
end
fprintf('\nCritères : >= %d cellules et >= %d SCE significatifs -> %d assemblée(s) exclue(s) sur %d testables\n', ...
    min_cellules_asm, min_SCE_asm, sum(R.exclue_criteres & R.testable), sum(R.testable));
R.testable = R.testable & ~R.exclue_criteres;

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
    % Intervalle de confiance à 95 % de l'intercept, ramené en proportion
    ic_b0   = [glme0.Coefficients.Lower(1), glme0.Coefficients.Upper(1)];
    ic_prop = 1 ./ (1 + exp(-(log(alpha / (1 - alpha)) + ic_b0)));
    fprintf('    IC 95 %% (GLMM) : %.1f %% à %.1f %% des assemblées\n', 100 * ic_prop);
catch err
    fprintf('A2. GLMM impossible (%s). Garder A1.\n', err.message);
end

% A3. Intervalle de confiance par rééchantillonnage des ANIMAUX (bootstrap
%     par grappes) : on tire 7 animaux avec remise, on recalcule la
%     proportion, 10 000 fois. Ne suppose aucun modèle ; plus prudent que
%     l'IC du GLMM quand l'effet varie beaucoup d'un animal à l'autre.
rng(1);
anim_u = unique(Rt.animal);
sig_par_animal = cellfun(@(a) Rt.sig(strcmp(Rt.animal, a)), anim_u, 'UniformOutput', false);
nB = 10000;
prop_boot = zeros(nB, 1);
for b = 1:nB
    tir = randi(numel(anim_u), numel(anim_u), 1);
    prop_boot(b) = mean(vertcat(sig_par_animal{tir}));
end
ic_boot = prctile(prop_boot, [2.5 97.5]);
fprintf(['A3. Proportion d''assemblées à p < %.2f : %.1f %% (%d / %d) ; ' ...
         'IC 95 %% (rééchantillonnage des animaux) : %.1f %% à %.1f %%\n'], ...
    alpha, 100 * mean(Rt.sig), sum(Rt.sig), height(Rt), 100 * ic_boot);

% A4. L'excès tient-il sans chaque animal ? (un résultat porté par un seul
%     animal doit être présenté comme tel)
fprintf('A4. Proportion sans chaque animal (attendu par hasard : %.0f %%) :\n', 100 * alpha);
for i_a = 1:numel(anim_u)
    garde_a = ~strcmp(Rt.animal, anim_u{i_a});
    fprintf('    sans %-8s : %2d / %3d = %4.1f %%\n', anim_u{i_a}, ...
        sum(Rt.sig(garde_a)), sum(garde_a), 100 * mean(Rt.sig(garde_a)));
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
R.categorie(R.exclue_criteres) = {'exclue_criteres'};

fprintf('\nC.  Classement (FDR %.0f %%) :\n', 100 * fdr_seuil);
disp(groupcounts(R, 'categorie'))
disp(groupcounts(R, {'animal', 'categorie'}))


%% ---------- C bis. FDR sur toutes les paires assemblée x couleur --------
% Description : QUELLE couleur est enrichie dans QUELLE assemblée.
% Utilise les p hypergéométriques brutes (une par couleur), sans
% permutation. La correction pour les 6 couleurs vient du fait que le FDR
% porte sur ~6 fois plus de tests que le classement C.
% Ne remplace ni le test principal (A) ni le classement des assemblées (C).
P = R.p_brut;                    % une ligne par assemblée : R G B RG RB GB
P(~R.testable, :) = NaN;         % mêmes assemblées que pour le reste
ok = ~isnan(P);                  % NaN = couleur exclue du champ (trop rare)
q_paire = nan(size(P));
q_paire(ok) = bh_qvalues(P(ok));

[ia, ic] = find(q_paire < fdr_seuil);
lin = sub2ind(size(P), ia, ic);
noms_couleur = {'R', 'G', 'B', 'RG', 'RB', 'GB'};
paires = table(R.animal(ia), R.champ(ia), R.session(ia), R.assemblee(ia), ...
    reshape(noms_couleur(ic), [], 1), R.n_colorees(ia), P(lin), q_paire(lin), ...
    'VariableNames', {'animal', 'champ', 'session', 'assemblee', 'couleur', ...
    'n_colorees', 'p_brut', 'q'});
paires = sortrows(paires, 'q');
fprintf('\nC bis. FDR sur %d paires assemblée x couleur : %d paire(s) avec q < %.2f\n', ...
    sum(ok(:)), height(paires), fdr_seuil);
disp(paires)
writetable(paires, fullfile(dossier_racine, 'paires_assemblee_couleur_FDR.csv'));


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
% E3. (Exploratoire) Biais selon la compacité spatiale de l'assemblée.
%     Si le biais n'apparaît que dans les assemblées compactes, la
%     proximité peut l'expliquer autant que la lignée.
if ismember('compacite_pct', Rt.Properties.VariableNames) && any(~isnan(Rt.compacite_pct))
    compacte = Rt.compacite_pct <= 5;
    fprintf('E3. Assemblées très compactes (percentile <= 5) : %d / %d à p < %.2f (%.1f %%)\n', ...
        sum(Rt.sig(compacte)), sum(compacte), alpha, 100 * mean(Rt.sig(compacte)));
    fprintf('    Autres assemblées                         : %d / %d à p < %.2f (%.1f %%)\n', ...
        sum(Rt.sig(~compacte & ~isnan(Rt.compacite_pct))), sum(~compacte & ~isnan(Rt.compacite_pct)), ...
        alpha, 100 * mean(Rt.sig(~compacte & ~isnan(Rt.compacite_pct))));
    % Le percentile de compacité dépend de la taille (les grandes assemblées
    % atteignent plus facilement des percentiles extrêmes) et la taille joue
    % sur la puissance du test : on ajuste donc sur la taille.
    Rt.compacte = double(compacte);
    try
        glme3 = fitglme(Rt(~isnan(Rt.compacite_pct), :), ...
            'sig ~ 1 + compacte + log_taille + (1|animal_c) + (1|champ_c)', 'Distribution', 'Binomial');
        fprintf('    Ajusté sur la taille : sig ~ compacte + log_taille + (1|animal) + (1|champ)\n');
        disp(glme3.Coefficients)
    catch err
        fprintf('    GLMM compacité impossible (%s).\n', err.message);
    end
end
% E3b. Même question avec la taille d'effet continue (rapport d'étalement :
%      1 = dispersée comme au hasard, < 1 = regroupée). Un coefficient
%      négatif = plus l'assemblée est regroupée, plus elle est biaisée.
if ismember('ratio_etalement', Rt.Properties.VariableNames) && any(~isnan(Rt.ratio_etalement))
    try
        glme3b = fitglme(Rt(~isnan(Rt.ratio_etalement), :), ...
            'sig ~ 1 + ratio_etalement + log_taille + (1|animal_c) + (1|champ_c)', 'Distribution', 'Binomial');
        fprintf('E3b. sig ~ ratio_etalement + log_taille + (1|animal) + (1|champ)\n');
        disp(glme3b.Coefficients)
    catch err
        fprintf('E3b. GLMM impossible (%s).\n', err.message);
    end
end


%% ---------- F. FIGURE : p_minp en fonction de la taille ------------------
% Gris : la meilleure p possible pour chaque assemblée. Une assemblée sous
% la ligne pointillée est biaisée ; une assemblée dont le point gris est
% au-dessus de la ligne ne pouvait pas l'être (non testable).
% Catégories de la figure, au seuil fdr_seuil_figure (mêmes q que le classement C)
cat_fig = repmat({'non_testable'}, height(R), 1);
cat_fig(R.testable & R.q <  fdr_seuil_figure)                                = {'biaisee'};
cat_fig(R.testable & R.q >= fdr_seuil_figure & R.p_minp >  seuil_non_biaise) = {'non_biaisee'};
cat_fig(R.testable & R.q >= fdr_seuil_figure & R.p_minp <= seuil_non_biaise) = {'indeterminee'};
n_biais_fig = sum(strcmp(cat_fig, 'biaisee'));
fprintf('\nF.  Figure (FDR %.0f %%) : %d assemblées biaisées sur %d testables\n', ...
    100 * fdr_seuil_figure, n_biais_fig, sum(R.testable));
disp(R(strcmp(cat_fig, 'biaisee'), {'animal', 'champ', 'session', 'assemblee', ...
    'couleur_dominante', 'n_colorees', 'p_minp', 'q'}))

sel = R.testable | figure_non_testables;
Rf = R(sel, :);
cf = cat_fig(sel);
fig = figure('Color', 'w');
h_leg = gobjects(0); noms_leg = {};
h_leg(end+1) = semilogy(Rf.n_colorees, Rf.p_min_atteignable, '.', ...
    'Color', [0.7 0.7 0.7], 'MarkerSize', 10); hold on
noms_leg{end+1} = 'p minimale atteignable';
cats = {'biaisee', 'indeterminee', 'non_biaisee', 'non_testable'};
noms_cats = {sprintf('biaisée (q FDR < %.2f)', fdr_seuil_figure), ...
             sprintf('indéterminée (q \\geq %.2f, p \\leq %.1f)', fdr_seuil_figure, seuil_non_biaise), ...
             sprintf('non biaisée (p > %.1f)', seuil_non_biaise), 'non testable'};
cols = [0.85 0.1 0.1; 0.95 0.6 0.1; 0.2 0.4 0.8; 0.5 0.5 0.5];
for k = 1:numel(cats)
    ix = strcmp(cf, cats{k});
    if ~any(ix), continue; end            % pas d'entrée de légende pour une catégorie vide
    h_leg(end+1) = semilogy(Rf.n_colorees(ix), Rf.p_minp(ix), 'o', ...
        'MarkerFaceColor', cols(k, :), 'MarkerEdgeColor', 'none', 'MarkerSize', 6); %#ok<SAGROW>
    noms_leg{end+1} = noms_cats{k}; %#ok<SAGROW>
end
h_leg(end+1) = yline(alpha, '--k');
noms_leg{end+1} = sprintf('seuil p = %.2f', alpha);
hold off
xlabel('Cellules colorées dans l''assemblée'); ylabel('p_{min-p}');
legend(h_leg, noms_leg, 'Location', 'southwest');
title(sprintf('Biais clonal des assemblées : %d biaisées sur %d testables (FDR %.0f %%)', ...
    n_biais_fig, sum(R.testable), 100 * fdr_seuil_figure));
exportgraphics(fig, fullfile(dossier_racine, 'biais_clonal_p_vs_taille.png'), 'Resolution', 200);


%% ---------- Export -------------------------------------------------------
writetable(R(:, ~ismember(R.Properties.VariableNames, ...
    {'animal_c', 'champ_c', 'condition_c'})), ...
    fullfile(dossier_racine, 'resultats_biais_clonal.csv'));
writetable(T_animal, fullfile(dossier_racine, 'resultats_par_animal.csv'));
save(fullfile(dossier_racine, 'resultats_biais_clonal.mat'), 'R', 'T_animal');

%% ---------- G. Figure des analyses complémentaires -----------------------
% (script séparé : figures_biais_clonal.m, dans le même dossier que celui-ci)
if exist('figures_biais_clonal', 'file')
    figures_biais_clonal
else
    fprintf('figures_biais_clonal.m introuvable : figure complémentaire non tracée.\n');
end


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
