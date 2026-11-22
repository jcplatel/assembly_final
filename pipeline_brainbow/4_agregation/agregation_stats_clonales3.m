function [fichier_resume, lignes_resume] = agregation_stats_clonales3(sessions, config, id_execution)
%% =========================================================================
%  ÉTAPE 4 : BIAIS CLONAL DES ASSEMBLÉES — toutes les sessions, tous les animaux
%  -------------------------------------------------------------------------
%  = agregation_stats_clonales2 (tests A à F inchangés dans leur logique),
%  avec :
%   - la liste des sessions vient de la table : plus de recherche de
%     fichiers, plus de doublons possibles
%   - vérification que TOUTES les sessions sont finies, par la même
%     exécution de l'étape 3 ; sinon arrêt avec la liste des sessions
%     manquantes (config.agreg.autoriser_partiel = true : on agrège quand
%     même, et les dossiers de sortie se terminent par _PROVISOIRE)
%   - les sessions exclues pour une raison prévue (pas d'assemblée, ...)
%     sont comptées et listées
%   - une assemblée sans n_SCE_signif est exclue (avant, NaN passait le critère)
%   - le hasard attendu est le taux_nul de chaque assemblée (<= 5 %) au lieu
%     de 5 % pour toutes : test plus juste, et plus puissant
%   - l'analyse est faite pour CHAQUE variante (référence / permutation) ;
%     la première de config.clonal.variantes est le test principal
%   - (3 nov.) RÉSUMÉ final : les deux phrases que les résultats permettent
%     d'écrire, avec les tests et les chiffres, pour le test principal et le
%     contrôle spatial -> agregation\resume_resultats.txt
%       1. biais clonal en général : Wilcoxon signé sur les animaux + GLMM
%          binomial (animal et session en effets aléatoires) + robustesse
%       2. lesquelles : Benjamini-Hochberg sur toutes les assemblées testables
%   - vérifie que l'étape 1 a les mêmes réglages pour toutes les sessions
%     (seuil des SCE, réglages des assemblées) : sinon, avertissement
%   - test du champ (B) et analyses exploratoires (D, E) seulement sur
%     demande : config.agreg.test_champ / config.agreg.analyses_exploratoires
%  Sorties dans <PathSave>\agregation\<variante>\
%  =========================================================================

fichier_resume = '';      % phrases récapitulatives (affichées à la fin de main_brainbow)
lignes_resume  = {};
dossier_agreg = fullfile(config.PathSave, 'agregation');
if ~isfolder(dossier_agreg), mkdir(dossier_agreg); end
nSess = height(sessions);

%% ---------- 1. Toutes les sessions sont-elles finies ? -------------------
statuts  = repmat({''}, nSess, 1);
ids      = repmat({''}, nSess, 1);
contenus = cell(nSess, 1);
manquantes = {};
for s = 1:nSess
    f = fullfile(sessions.dossier{s}, 'clonal', 'stats_clonales.mat');
    if ~isfile(f)
        manquantes{end+1} = sprintf('%s : pas de stats_clonales.mat', sessions.nom{s}); %#ok<AGROW>
        continue
    end
    C = load(f);
    statuts{s} = C.statut;
    ids{s}     = C.id_execution;
    contenus{s} = C;
end
% Une seule exécution de l'étape 3 pour toutes les sessions
if isempty(id_execution)
    id_presents = unique(ids(~cellfun(@isempty, ids)));
    if numel(id_presents) > 1
        manquantes{end+1} = sprintf(['les résultats viennent de %d exécutions différentes de ' ...
            'l''étape 3 (%s) : relancer l''étape 3 pour toutes les sessions'], ...
            numel(id_presents), strjoin(id_presents, ', '));
    end
else
    for s = 1:nSess
        if ~isempty(ids{s}) && ~strcmp(ids{s}, id_execution)
            manquantes{end+1} = sprintf('%s : résultat d''une autre exécution (%s)', ...
                sessions.nom{s}, ids{s}); %#ok<AGROW>
        end
    end
end
suffixe = '';
if ~isempty(manquantes)
    message = sprintf('%d problème(s) :\n   %s', numel(manquantes), strjoin(manquantes, '\n   '));
    if config.agreg.autoriser_partiel
        warning('AGRÉGATION PROVISOIRE. %s', message);
        suffixe = '_PROVISOIRE';
    else
        error(['Agrégation refusée : toutes les sessions doivent être finies.\n%s\n' ...
            'Relancer l''étape 3 pour ces sessions (ou config.agreg.autoriser_partiel = true ' ...
            'pour un aperçu provisoire).'], message);
    end
end

% Sessions exclues pour une raison prévue
est_ok     = strcmp(statuts, 'ok');
est_exclue = ~cellfun(@isempty, statuts) & ~est_ok;
T_sessions = table(sessions.nom, statuts, 'VariableNames', {'session', 'statut'});
writetable(T_sessions, fullfile(dossier_agreg, ['statut_sessions' suffixe '.csv']));
fprintf('\n%d sessions : %d analysées, %d exclues (raison prévue), %d manquantes\n', ...
    nSess, sum(est_ok), sum(est_exclue), nSess - sum(est_ok) - sum(est_exclue));
if any(est_exclue)
    disp(T_sessions(est_exclue, :))
end
if ~any(est_ok)
    error('Aucune session analysée : rien à agréger.');
end

% Biais de détection, une valeur par session
p_detection_sessions = cellfun(@(C) C.p_detection, contenus(est_ok));
p_F0_sessions        = cellfun(@(C) C.p_F0, contenus(est_ok));
p_compo_sessions     = cellfun(@(C) C.p_composition_membres, contenus(est_ok));
fprintf(['Biais de détection selon la couleur, sessions à p < 0.05 :\n' ...
         '   bruit de détection : %d / %d\n   F0                 : %d / %d\n' ...
         '   couleurs membres vs non-membres : %d / %d  (biais de détection OU effet clonal)\n'], ...
    sum(p_detection_sessions < 0.05), numel(p_detection_sessions), ...
    sum(p_F0_sessions < 0.05), numel(p_F0_sessions), ...
    sum(p_compo_sessions < 0.05), numel(p_compo_sessions));
writetable(table(sessions.nom(est_ok), p_detection_sessions, p_F0_sessions, p_compo_sessions, ...
    'VariableNames', {'session', 'p_detection', 'p_F0', 'p_composition_membres'}), ...
    fullfile(dossier_agreg, ['biais_detection' suffixe '.csv']));

%% ---------- 1b. analysis_brainbow_ALL.xlsx : une ligne par session -----
% (ancien export_data_brainbow, réécrit entièrement : plus de doublons)
contenus_ok = contenus(est_ok);
lignes = cellfun(@(C) cell2table(struct2cell(C.resume_session)', ...
    'VariableNames', fieldnames(C.resume_session)'), contenus_ok, 'UniformOutput', false);
path_global = fullfile(dossier_agreg, ['analysis_brainbow_ALL' suffixe '.xlsx']);
if isfile(path_global), delete(path_global); end
writetable(vertcat(lignes{:}), path_global);

%% ---------- 1c. Contrôles du champ, sans les assemblées (5 nov.) --------
% Couleurs regroupées dans l'espace ? La couleur change-t-elle la détection
% ou l'activité ? -> agregation\controles_champ\ (CSV, figures) et résumé
try
    lignes_controles = controles_champ(contenus_ok, sessions(est_ok, :), config, dossier_agreg);
catch err
    lignes_controles = {sprintf('Contrôles du champ non faits : %s', err.message)};
    warning('controles_champ : %s', err.message);
end

%% ---------- 2. Une analyse par variante ---------------------------------
noms_variantes = contenus_ok{1}.noms_variantes;
bilans = cell(numel(noms_variantes), 1);
lignes_compacite = {};
for v = 1:numel(noms_variantes)
    liste_res = cellfun(@(C) avec_composition(C, v), contenus_ok, 'UniformOutput', false);
    R = vertcat(liste_res{:});
    if v == 1
        % (5 nov.) compacité spatiale des assemblées, sans couleur (une fois)
        try
            lignes_compacite = compacite_assemblees(R, sessions(est_ok, :), config, dossier_agreg);
        catch err
            lignes_compacite = {sprintf('Compacité des assemblées non faite : %s', err.message)};
            warning('compacite_assemblees : %s', err.message);
        end
    end
    if v == 1
        nom_affiche = [noms_variantes{v} ' (TEST PRINCIPAL)'];
    else
        nom_affiche = [noms_variantes{v} ' (contrôle)'];
    end
    bilans{v} = analyser_variante(R, config, fullfile(dossier_agreg, [noms_variantes{v} suffixe]), nom_affiche);
end

% (22 nov., demandé par l'utilisateur) profil des assemblées biaisées et
% taille des événements CASCADE -> agregation\profil_assemblees\
lignes_profil = {};
try
    lignes_profil = profil_assemblees(fullfile(dossier_agreg, [noms_variantes{1} suffixe], ...
        'resultats_biais_clonal.csv'), sessions(est_ok, :), config, dossier_agreg);
catch err
    lignes_profil = {sprintf('Profil des assemblées non fait : %s', err.message)};
    warning('profil_assemblees : %s', err.message);
end

% Réglages de cette étape 4 (bouton « Analysis parameters » de l'interface)
agreg = config.agreg;
agreg.seuilstat = config.clonal.seuilstat;
date_etape4 = datestr(now, 'yyyy-mm-dd HH:MM'); %#ok<NASGU>
save(fullfile(dossier_agreg, 'parametres_etape4.mat'), 'agreg', 'date_etape4');

%% ---------- 3. Résumé : ce que les résultats permettent de dire ---------
lignes_etape1 = verifier_etape1(sessions(est_ok, :));
fichier_resume = fullfile(dossier_agreg, ['resume_resultats' suffixe '.txt']);
if ~isempty(lignes_compacite)
    lignes_controles = [lignes_controles(:)', lignes_compacite(:)'];
end
if ~isempty(lignes_profil)
    lignes_controles = [lignes_controles(:)', lignes_profil(:)'];
end
lignes_resume = ecrire_resume(bilans, noms_variantes, config, fichier_resume, lignes_etape1, ...
    lignes_controles);
end


%% =========================================================================
function bilan = analyser_variante(R, config, dossier_racine, nom_variante)
% Tests A à F de agregation_stats_clonales2 pour une variante.
% bilan : les chiffres du résumé final (phrases 1 et 2)
bilan = struct('n_testables', 0, 'n_sig', 0, 'n_attendu', NaN, 'n_animaux', 0, ...
    'p_wilcoxon', NaN, 'p_wilcoxon_min', NaN, 'OR', NaN, 'OR_ic', [NaN NaN], 'p_glmm', NaN, ...
    'exces_ic', [NaN NaN], 'exces_sans_animal', [NaN NaN], 'n_biaisees', 0, ...
    'n_sessions_biaisees', 0, 'n_animaux_biaisees', 0, 'couleurs_biaisees', '', 'fdr', NaN, ...
    'score_moyen', NaN, 'p_wil_score', NaN, 'lmm_score', [NaN NaN NaN], 'p_lmm_score', NaN, ...
    'pi1', NaN, 'pi1_ic', [NaN NaN], 'n_score', 0);
bilan.ids_biaisees = {};
faire_test_champ   = isfield(config.agreg, 'test_champ') && config.agreg.test_champ;
faire_exploratoire = isfield(config.agreg, 'analyses_exploratoires') && config.agreg.analyses_exploratoires;

%% ---------- PARAMÈTRES (fixés AVANT de regarder les résultats) ---------
alpha            = config.clonal.seuilstat;   % le même seuil que taux_nul
fdr_seuil        = config.agreg.fdr_seuil;
seuil_non_biaise = config.agreg.seuil_non_biaise;
min_cellules_asm = config.agreg.min_cellules_asm;
min_SCE_asm      = config.agreg.min_SCE_asm;
fdr_seuil_figure = config.agreg.fdr_seuil_figure;
figure_non_testables = config.agreg.figure_non_testables;
nPerm            = config.clonal.nPerm;
if ~isfolder(dossier_racine), mkdir(dossier_racine); end
fprintf('\n################ VARIANTE %s ################\n', nom_variante);

% Identifiants
R.session_u   = strcat(R.animal, '_', R.champ, '_', R.session, '_', R.condition);
% (5 nov.) Un même plan (plane0...) revu un autre jour n'est en général PAS
% à la même profondeur : ce sont d'autres cellules. Chaque session est donc
% son propre champ de vue ; l'effet aléatoire « champ » des modèles mixtes
% (GLMM A2, LMM A5, figures) est la session. (Avant : animal + id_champ,
% qui regroupait les plane0 de jours différents.)
R.champ_u     = R.session_u;
R.animal_c    = categorical(R.animal);
R.champ_c     = categorical(R.champ_u);
R.condition_c = categorical(R.condition);
R.log_taille  = log(max(R.n_colorees, 1));

% Critères de validité (ACTIVITÉ seulement), appliqués AVANT tous les tests
R.exclue_criteres = R.n_total < min_cellules_asm;
if min_SCE_asm > 0
    R.exclue_criteres = R.exclue_criteres | ~(R.n_SCE_signif >= min_SCE_asm);   % NaN -> exclue
end
fprintf('Critères : >= %d cellules et >= %d SCE significatifs -> %d assemblée(s) exclue(s) sur %d testables\n', ...
    min_cellules_asm, min_SCE_asm, sum(R.exclue_criteres & R.testable), sum(R.testable));
R.testable = R.testable & ~R.exclue_criteres;
fprintf('%d animaux, %d sessions (chacune son champ de vue), %d assemblées (%d testables)\n', ...
    numel(unique(R.animal)), numel(unique(R.session_u)), height(R), sum(R.testable));

Rt = R(R.testable, :);
if isempty(Rt)
    warning('Aucune assemblée testable pour cette variante.');
    return
end
Rt.sig = double(Rt.p_minp < alpha);

%% ---------- A. TEST PRINCIPAL : plus d'assemblées biaisées que le hasard ?
% A1. Par animal : proportion observée - proportion attendue (taux_nul)
[G, animaux] = findgroups(Rt.animal);
prop_animal    = splitapply(@mean,  Rt.sig, G);
attendu_animal = splitapply(@mean,  Rt.taux_nul, G);
n_animal       = splitapply(@numel, Rt.sig, G);
T_animal = table(animaux, n_animal, prop_animal, attendu_animal, ...
    'VariableNames', {'animal', 'n_assemblees_testables', 'proportion_biaisees', 'proportion_attendue'});
disp(T_animal)
p_wilcoxon = signrank(prop_animal - attendu_animal, 0, 'tail', 'right');
fprintf(['A1. Assemblées biaisées : %.1f %% observé, %.1f %% attendu par hasard\n' ...
         '    Wilcoxon sur les animaux (observé - attendu > 0) : p = %.4f (n = %d animaux)\n'], ...
    100 * mean(prop_animal), 100 * mean(attendu_animal), p_wilcoxon, numel(prop_animal));
% NB : avec n animaux, la plus petite p possible du Wilcoxon est 1/2^n
bilan.n_testables    = height(Rt);
bilan.n_sig          = sum(Rt.sig);
bilan.n_attendu      = sum(Rt.taux_nul);
bilan.n_animaux      = numel(prop_animal);
bilan.p_wilcoxon     = p_wilcoxon;
bilan.p_wilcoxon_min = 1 / 2^numel(prop_animal);

% A2. GLMM binomial, offset = logit(taux_nul) de chaque assemblée :
%     exp(intercept) = rapport des cotes observé / hasard (1 = pas d'effet)
try
    taux = min(max(Rt.taux_nul, 1 / nPerm), 1 - 1 / nPerm);
    offs = log(taux ./ (1 - taux));
    glme0 = fitglme(Rt, 'sig ~ 1 + (1|animal_c) + (1|champ_c)', ...
        'Distribution', 'Binomial', 'Offset', offs);
    b0 = glme0.Coefficients.Estimate(1);
    ic_b0 = [glme0.Coefficients.Lower(1), glme0.Coefficients.Upper(1)];
    fprintf(['A2. GLMM : rapport des cotes observé / hasard = %.2f (IC 95 %% %.2f à %.2f), ' ...
             'p bilatéral = %.4g\n'], exp(b0), exp(ic_b0), glme0.Coefficients.pValue(1));
    bilan.OR     = exp(b0);
    bilan.OR_ic  = exp(ic_b0);
    bilan.p_glmm = glme0.Coefficients.pValue(1);
catch err
    fprintf('A2. GLMM impossible (%s). Garder A1.\n', err.message);
end

% A3. IC par rééchantillonnage des ANIMAUX (bootstrap par grappes)
rng(1);
anim_u = unique(Rt.animal);
sig_par_animal = cellfun(@(a) Rt.sig(strcmp(Rt.animal, a)), anim_u, 'UniformOutput', false);
att_par_animal = cellfun(@(a) Rt.taux_nul(strcmp(Rt.animal, a)), anim_u, 'UniformOutput', false);
nB = 10000;
exces_boot = zeros(nB, 1);
for b = 1:nB
    tir = randi(numel(anim_u), numel(anim_u), 1);
    exces_boot(b) = mean(vertcat(sig_par_animal{tir})) - mean(vertcat(att_par_animal{tir}));
end
ic_boot = prctile(exces_boot, [2.5 97.5]);
bilan.exces_ic = 100 * ic_boot;
fprintf(['A3. Assemblées à p < %.2f : %d / %d (%.1f %%), attendu %.1f ; excès %.1f points, ' ...
         'IC 95 %% (animaux rééchantillonnés) : %.1f à %.1f points\n'], ...
    alpha, sum(Rt.sig), height(Rt), 100 * mean(Rt.sig), sum(Rt.taux_nul), ...
    100 * (mean(Rt.sig) - mean(Rt.taux_nul)), 100 * ic_boot);

% A4. L'excès tient-il sans chaque animal ?
fprintf('A4. Sans chaque animal (observé / attendu) :\n');
exces_sans = nan(numel(anim_u), 1);
for i_a = 1:numel(anim_u)
    garde_a = ~strcmp(Rt.animal, anim_u{i_a});
    exces_sans(i_a) = 100 * (mean(Rt.sig(garde_a)) - mean(Rt.taux_nul(garde_a)));
    fprintf('    sans %-8s : %2d / %3d = %4.1f %% (attendu %4.1f %%)\n', anim_u{i_a}, ...
        sum(Rt.sig(garde_a)), sum(garde_a), 100 * mean(Rt.sig(garde_a)), 100 * mean(Rt.taux_nul(garde_a)));
end
bilan.exces_sans_animal = [min(exces_sans), max(exces_sans)];

%% ---------- A5. Toutes les p, pas seulement p < 0,05 (score continu) ------
% score_minp (étape 3) : p_minp en score centré-réduit sur les permutations
% de chaque assemblée (0 sans effet). Moyenne par animal, Wilcoxon sur les
% animaux ; LMM sur les assemblées (animal et session aléatoires).
% Simulation (8 animaux, 16 champs, 128 assemblées ; 16 études par cas) :
%   aucun effet                                : 0 % de faux positifs (A1 : 0 %)
%   15 % des assemblées, effet fort            : A5 81 %, A1 100 % des études
%   60 % des assemblées, effet faible          : A5 100 %, A1 88 %
% A1 est meilleur quand peu d'assemblées ont un effet net ; A5 quand beaucoup
% ont un petit effet. Les deux sont donnés.
%% ---------- A6. Proportion d'assemblées clonales --------------------------
% pi1 = 1 - (assemblées à p_minp > 0,5) / (nombre attendu sans effet, somme
% des p0_sup). Sans effet : -3,5 % +/- 5,9 en simulation ; 30 % d'assemblées
% clonales : 24 % estimés. C'est un minimum : une assemblée faiblement
% clonale dont p_minp > 0,5 compte comme non clonale.
if all(ismember({'score_minp', 'p0_sup'}, Rt.Properties.VariableNames)) && any(~isnan(Rt.score_minp))
    [Gs5, ~] = findgroups(Rt.animal);
    score_animal = splitapply(@(x) mean(x, 'omitnan'), Rt.score_minp, Gs5);
    bilan.n_score     = sum(~isnan(Rt.score_minp));
    bilan.score_moyen = mean(Rt.score_minp, 'omitnan');
    bilan.p_wil_score = signrank(score_animal, 0, 'tail', 'right');
    fprintf(['A5. Toutes les p (score continu, 0 sans effet) : score moyen %.2f ; ' ...
        'Wilcoxon sur les animaux : p = %.4f\n'], bilan.score_moyen, bilan.p_wil_score);
    try
        lme5 = fitlme(Rt(~isnan(Rt.score_minp), :), 'score_minp ~ 1 + (1|animal_c) + (1|champ_c)');
        % (5 nov.) ddl de Satterthwaite : avec 8 animaux, le test par défaut
        % (ddl des résidus) est trop optimiste dès que les animaux diffèrent
        [~, ~, c5] = fixedEffects(lme5, 'DFMethod', 'satterthwaite');
        bilan.lmm_score   = [c5.Estimate(1), c5.Lower(1), c5.Upper(1)];
        bilan.p_lmm_score = c5.pValue(1);
        fprintf('    LMM : score %.2f (IC 95 %% %.2f à %.2f), p bilatéral = %.4g\n', bilan.lmm_score, bilan.p_lmm_score);
    catch err
        fprintf('    LMM impossible (%s).\n', err.message);
    end
    lam = 0.5;
    ok6 = ~isnan(Rt.p0_sup);
    [Ga6, anim6] = findgroups(Rt.animal(ok6));
    n_sup_a = splitapply(@sum, double(Rt.p_minp(ok6) > lam), Ga6);
    e_sup_a = splitapply(@sum, Rt.p0_sup(ok6), Ga6);
    bilan.pi1 = 1 - sum(n_sup_a) / sum(e_sup_a);
    pi1_boot = zeros(nB, 1);
    for b = 1:nB
        tir = randi(numel(anim6), numel(anim6), 1);
        pi1_boot(b) = 1 - sum(n_sup_a(tir)) / sum(e_sup_a(tir));
    end
    bilan.pi1_ic = prctile(pi1_boot, [2.5 97.5]);
    fprintf(['A6. Proportion estimée d''assemblées clonales : %.1f %% (IC 95 %% en rééchantillonnant ' ...
        'les animaux : %.1f à %.1f %%) ; un minimum\n'], 100 * bilan.pi1, 100 * bilan.pi1_ic);
else
    fprintf('A5-A6 : score_minp absent (refaire l''étape 3 avec ce pipeline).\n');
end

%% ---------- B. Test secondaire : test global par champ (sur demande) ----
if faire_test_champ
    [Gs, ~] = findgroups(R.session_u);
    z_sess   = splitapply(@(z) z(1), R.z_champ, Gs);
    anim_ses = splitapply(@(x) x(1), R.animal, Gs);
    [Ga, ~]  = findgroups(anim_ses);
    z_animal = splitapply(@(z) mean(z, 'omitnan'), z_sess, Ga);   % champs sans assemblée testable : NaN
    fprintf('B.  z moyen du test par champ : %.2f ; Wilcoxon > 0 sur les animaux : p = %.4f\n', ...
        mean(z_animal, 'omitnan'), signrank(z_animal, 0, 'tail', 'right'));
end

%% ---------- C. CLASSEMENT de chaque assemblée (FDR) ---------------------
R.q = nan(height(R), 1);
R.q(R.testable) = bh_qvalues(R.p_minp(R.testable));
R.categorie = repmat({'non_testable'}, height(R), 1);
R.categorie(R.testable & R.q <  fdr_seuil)                                = {'biaisee'};
R.categorie(R.testable & R.q >= fdr_seuil & R.p_minp >  seuil_non_biaise) = {'non_biaisee'};
R.categorie(R.testable & R.q >= fdr_seuil & R.p_minp <= seuil_non_biaise) = {'indeterminee'};
R.categorie(R.exclue_criteres) = {'exclue_criteres'};
fprintf('\nC.  Classement (FDR %.0f %%) :\n', 100 * fdr_seuil);
disp(groupcounts(R, 'categorie'))
disp(groupcounts(R, {'animal', 'categorie'}))
m_b = strcmp(R.categorie, 'biaisee');
bilan.fdr                 = fdr_seuil;
bilan.n_biaisees          = sum(m_b);
bilan.n_sessions_biaisees = numel(unique(R.session_u(m_b)));
bilan.n_animaux_biaisees  = numel(unique(R.animal(m_b)));
bilan.ids_biaisees        = strcat(R.session_u(m_b), '#', ...
    arrayfun(@(a) sprintf('%d', a), R.assemblee(m_b), 'UniformOutput', false));
if any(m_b)
    [coul_u, ~, ic_c] = unique(R.couleur_dominante(m_b));
    n_coul = accumarray(ic_c, 1);
    [n_coul, ordre_c] = sort(n_coul, 'descend');
    bilan.couleurs_biaisees = strjoin(arrayfun(@(k) sprintf('%s %d', coul_u{ordre_c(k)}, n_coul(k)), ...
        1:numel(n_coul), 'UniformOutput', false), ', ');
end

%% ---------- C ter. FDR sur les assemblées qui pouvaient gagner (10 nov.) --
% Sensibilité, demandée par l'utilisateur. Le FDR de C porte sur toutes les
% assemblées testables, dont beaucoup ne peuvent pas atteindre un petit p
% (peu de cellules colorées, couleur fréquente). Ici, Benjamini-Hochberg sur
% les seules assemblées dont le meilleur p possible (p_min_atteignable) est
% sous config.agreg.fdr_filtre_pmin. Ce tri ne regarde pas le résultat
% (taille de l'assemblée et couleurs du champ seulement) : il reste valide.
% Le classement principal (C) ne change pas.
filtre_pmin = 0.001;
if isfield(config.agreg, 'fdr_filtre_pmin'), filtre_pmin = config.agreg.fdr_filtre_pmin; end
R.fdr_filtre = R.testable & ~R.exclue_criteres & R.p_min_atteignable < filtre_pmin;
R.q_filtre = nan(height(R), 1);
if any(R.fdr_filtre), R.q_filtre(R.fdr_filtre) = bh_qvalues(R.p_minp(R.fdr_filtre)); end
m_f = R.fdr_filtre & R.q_filtre < fdr_seuil;
bilan.filtre_pmin        = filtre_pmin;
bilan.n_filtre           = sum(R.fdr_filtre);
bilan.n_biaisees_filtre  = sum(m_f);
bilan.n_animaux_filtre   = numel(unique(R.animal(m_f)));
bilan.n_nouvelles_filtre = sum(m_f & ~m_b);
fprintf(['\nC ter. FDR %.0f %% sur les %d assemblées à p_min_atteignable < %g (sur %d testables) : ' ...
    '%d biaisée(s), dont %d absente(s) du classement principal\n'], 100 * fdr_seuil, bilan.n_filtre, ...
    filtre_pmin, sum(R.testable & ~R.exclue_criteres), bilan.n_biaisees_filtre, bilan.n_nouvelles_filtre);
if any(m_f)
    disp(R(m_f, {'animal', 'session', 'assemblee', 'n_total', 'n_colorees', 'couleur_dominante', ...
        'p_minp', 'q', 'q_filtre'}))
end

%% ---------- C bis. FDR sur toutes les paires assemblée x couleur --------
P = R.p_brut;
P(~R.testable, :) = NaN;
ok = ~isnan(P);
q_paire = nan(size(P));
q_paire(ok) = bh_qvalues(P(ok));
[ia, ic] = find(q_paire < fdr_seuil);
lin = sub2ind(size(P), ia, ic);
noms_couleur = {'R', 'G', 'B', 'RG', 'RB', 'GB'};
paires = table(R.animal(ia), R.champ(ia), R.session(ia), R.assemblee(ia), ...
    reshape(noms_couleur(ic), [], 1), R.n_colorees(ia), P(lin), q_paire(lin), ...
    'VariableNames', {'animal', 'champ', 'session', 'assemblee', 'couleur', 'n_colorees', 'p_brut', 'q'});
paires = sortrows(paires, 'q');
fprintf('\nC bis. FDR sur %d paires assemblée x couleur : %d paire(s) avec q < %.2f\n', ...
    sum(ok(:)), height(paires), fdr_seuil);
disp(paires)
writetable(paires, fullfile(dossier_racine, 'paires_assemblee_couleur_FDR.csv'));

%% ---------- D et E (sur demande : config.agreg.analyses_exploratoires) --
% Attention : le résultat de ces modèles est la significativité ; une grande
% assemblée a plus de cellules colorées, donc plus de puissance : l'effet de
% la taille y est en partie mécanique.
if faire_exploratoire
    %% ---------- D. CARACTÉRISATION : biaisées contre non biaisées -----------
    % La taille est TOUJOURS dans le modèle. Jamais de comparaison sur la
    % composition en couleurs (circulaire).
    C = R(ismember(R.categorie, {'biaisee', 'non_biaisee'}), :);
    C.biaisee = double(strcmp(C.categorie, 'biaisee'));
    a_score = any(~isnan(C.score_sequence));
    a_cond  = numel(unique(C.condition)) > 1;
    termes  = {'log_taille'};
    if a_score && a_cond
        termes{end+1} = 'score_sequence*condition_c';
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
    if ~a_score
        fprintf('    (score_sequence absent : la séquentialité n''est pas dans le modèle)\n');
    end

    %% ---------- E. ROBUSTESSE ------------------------------------------------
    Rt.mlogp = -log10(Rt.p_minp);
    formule_E = ['mlogp ~ 1 + ' strjoin(termes, ' + ') ' + (1|animal_c) + (1|champ_c)'];
    try
        lme2 = fitlme(Rt, formule_E);
        fprintf('\nE1. %s\n', formule_E);
        disp(lme2.Coefficients)
    catch err
        fprintf('E1. LMM impossible (%s).\n', err.message);
    end
    for s = [0.05 0.10 0.20]
        fprintf('E2. FDR %.0f %% : %d assemblées biaisées\n', 100 * s, sum(R.q < s));
    end
    if any(~isnan(Rt.compacite_pct))
        compacte = Rt.compacite_pct <= 5;
        autre = ~compacte & ~isnan(Rt.compacite_pct);
        fprintf('E3. Assemblées très compactes (percentile <= 5) : %d / %d à p < %.2f (%.1f %%)\n', ...
            sum(Rt.sig(compacte)), sum(compacte), alpha, 100 * mean(Rt.sig(compacte)));
        fprintf('    Autres assemblées                         : %d / %d à p < %.2f (%.1f %%)\n', ...
            sum(Rt.sig(autre)), sum(autre), alpha, 100 * mean(Rt.sig(autre)));
        Rt.compacte = double(compacte);
        try
            glme3 = fitglme(Rt(~isnan(Rt.compacite_pct), :), ...
                'sig ~ 1 + compacte + log_taille + (1|animal_c) + (1|champ_c)', 'Distribution', 'Binomial');
            disp(glme3.Coefficients)
        catch err
            fprintf('    GLMM compacité impossible (%s).\n', err.message);
        end
    end
    if any(~isnan(Rt.ratio_etalement))
        try
            glme3b = fitglme(Rt(~isnan(Rt.ratio_etalement), :), ...
                'sig ~ 1 + ratio_etalement + log_taille + (1|animal_c) + (1|champ_c)', 'Distribution', 'Binomial');
            fprintf('E3b. sig ~ ratio_etalement + log_taille + (1|animal) + (1|champ)\n');
            disp(glme3b.Coefficients)
        catch err
            fprintf('E3b. GLMM impossible (%s).\n', err.message);
        end
    end
end

%% ---------- F. FIGURE : p_minp en fonction de la taille ------------------
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
fig = figure('Color', 'w', 'Visible', 'off');
h_leg = gobjects(0); noms_leg = {};
h_leg(end+1) = semilogy(Rf.n_colorees, Rf.p_min_atteignable, '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 10);
hold on
noms_leg{end+1} = 'p minimale atteignable';
cats = {'biaisee', 'indeterminee', 'non_biaisee', 'non_testable'};
noms_cats = {sprintf('biaisée (q FDR < %.2f)', fdr_seuil_figure), ...
    sprintf('indéterminée (q \\geq %.2f, p \\leq %.1f)', fdr_seuil_figure, seuil_non_biaise), ...
    sprintf('non biaisée (p > %.1f)', seuil_non_biaise), 'non testable'};
cols = [0.85 0.1 0.1; 0.95 0.6 0.1; 0.2 0.4 0.8; 0.5 0.5 0.5];
for k = 1:numel(cats)
    ix = strcmp(cf, cats{k});
    if ~any(ix), continue; end
    h_leg(end+1) = semilogy(Rf.n_colorees(ix), Rf.p_minp(ix), 'o', ...
        'MarkerFaceColor', cols(k, :), 'MarkerEdgeColor', 'none', 'MarkerSize', 6); %#ok<AGROW>
    noms_leg{end+1} = noms_cats{k}; %#ok<AGROW>
end
h_leg(end+1) = yline(alpha, '--k');
noms_leg{end+1} = sprintf('seuil p = %.2f', alpha);
hold off
xlabel('Cellules colorées dans l''assemblée'); ylabel('p_{min-p}');
legend(h_leg, noms_leg, 'Location', 'southwest');
title({sprintf('Biais clonal des assemblées : %d biaisées sur %d testables (FDR %.0f %%)', ...
    n_biais_fig, sum(R.testable), 100 * fdr_seuil_figure), ['Variante ' nom_variante]}, ...
    'Interpreter', 'none');
exportgraphics(fig, fullfile(dossier_racine, 'biais_clonal_p_vs_taille.png'), 'Resolution', 200);
close(fig)

%% ---------- Export -------------------------------------------------------
writetable(R(:, ~ismember(R.Properties.VariableNames, {'animal_c', 'champ_c', 'condition_c'})), ...
    fullfile(dossier_racine, 'resultats_biais_clonal.csv'));
writetable(T_animal, fullfile(dossier_racine, 'resultats_par_animal.csv'));
save(fullfile(dossier_racine, 'resultats_biais_clonal.mat'), 'R', 'T_animal', 'bilan');

%% ---------- G. Figure des analyses complémentaires ---------------------
% figures\figures_biais_clonal.m (script : il lit resultats_biais_clonal.mat
% de cette variante et reprend les seuils ci-dessus). Il écrit
% clonal_bias_complementary_analyses.png dans le dossier de la variante ;
% figure_assemblees_par_session écrit assemblies_per_session.png.
% Seule la figure du test principal reste ouverte.
figures_avant = findall(groot, 'Type', 'figure');
try
    figures_biais_clonal
catch err
    warning('figures_biais_clonal non exécuté : %s', err.message);
end
% Toutes les sessions : assemblées trouvées, biaisées et couleur (assemblies_per_session.png)
try
    figure_assemblees_par_session(dossier_racine, fdr_seuil_figure, alpha, nom_variante);
catch err
    warning('figure_assemblees_par_session non exécutée : %s', err.message);
end
% Assemblées biaisées de chaque animal, un rond par assemblée (biased_assemblies_per_animal.png)
try
    figure_biaisees_par_animal(dossier_racine, fdr_seuil_figure, alpha, nom_variante);
catch err
    warning('figure_biaisees_par_animal non exécutée : %s', err.message);
end
% Figures pour l'article : tailles par animal, quantification du biais (PNG + PDF)
try
    figures_article(dossier_racine, fdr_seuil_figure, alpha, nom_variante);
catch err
    warning('figures_article non exécutée : %s', err.message);
end
toutes = findall(groot, 'Type', 'figure');
nouvelles = toutes(arrayfun(@(h) ~any(h == figures_avant), toutes));
principal = contains(nom_variante, 'PRINCIPAL');
for k = 1:numel(nouvelles)
    nouvelles(k).Name = ['Clonal bias: ' nom_variante];
    if ~principal
        close(nouvelles(k));
    end
end
end


%% =========================================================================
function r = avec_composition(C, v)
% résultats d'une session pour la variante v, avec la composition en
% couleurs de l'urne de référence (Kcol, R G B RG RB GB) sur chaque ligne
r = C.res_variantes{v};
if isfield(C, 'Kcol_variantes') && numel(C.Kcol_variantes) >= v && ~isempty(C.Kcol_variantes{v})
    r.Kcol_urne = repmat(double(C.Kcol_variantes{v}(:)'), height(r), 1);
else
    r.Kcol_urne = nan(height(r), 6);
end
end


%% =========================================================================
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


%% =========================================================================
function lignes = verifier_etape1(sessions)
% L'étape 1 a-t-elle les mêmes réglages pour toutes les sessions ? (seuil
% des SCE, réglages des assemblées du K retenu). Sinon les sessions ne sont
% pas comparables : avertissement dans la console et dans le résumé.
nS = height(sessions);
regles = repmat({''}, nS, 1);
asm = cell(nS, 1);
for s = 1:nS
    ds = sessions.dossier{s};
    if ismember('dossier_source', sessions.Properties.VariableNames), ds = sessions.dossier_source{s}; end
    f = fullfile(ds, 'pretraitement.mat');
    if ~isfile(f), continue; end
    try
        P = load(f, 'opts');
        regles{s} = texte_regle_sce(P.opts);
        fch = fullfile(sessions.dossier{s}, 'choix_K.mat');
        if isfile(fch)
            C = load(fch, 'best_K', 'statut_choix');
            fa = fichier_assemblees(sessions.dossier{s}, ds, C.best_K);   % v1 ou v2 selon l'analyse
            if strcmp(C.statut_choix, 'ok') && isfile(fa) && ismember('config_asm', who('-file', fa))
                A = load(fa, 'config_asm');
                asm{s} = A.config_asm;
            end
        end
    catch
    end
end
lignes = {};
lues = ~cellfun(@isempty, regles);
[u, ~, ic] = unique(regles(lues));
noms = sessions.nom(lues);
if numel(u) == 1
    lignes{end+1} = sprintf('Seuil des SCE : le même pour les %d sessions (%s).', sum(lues), u{1});
elseif numel(u) > 1
    lignes{end+1} = sprintf(['ATTENTION : le seuil des SCE de l''étape 1 n''est PAS le même pour toutes ' ...
        'les sessions (%d règles) : refaire l''étape 1 avec une seule règle avant de conclure.'], numel(u));
    for k = 1:numel(u)
        ex = noms(ic == k);
        lignes{end+1} = sprintf('   %-40s %2d session(s), ex. %s', u{k}, numel(ex), ...
            strjoin(ex(1:min(2, end)), ', ')); %#ok<AGROW>
    end
end
% Réglages des assemblées : versions différentes ?
a_lus = find(~cellfun(@isempty, asm));
if ~isempty(a_lus)
    ref = asm{a_lus(1)};
    differents = {};
    n_diff = 0;
    for s = a_lus(2:end)'
        if ~isequaln(asm{s}, ref)
            n_diff = n_diff + 1;
            champs = union(fieldnames(ref), fieldnames(asm{s}));
            for c = 1:numel(champs)
                if ~isfield(ref, champs{c}) || ~isfield(asm{s}, champs{c}) || ...
                        ~isequaln(ref.(champs{c}), asm{s}.(champs{c}))
                    differents{end+1} = champs{c}; %#ok<AGROW>
                end
            end
        end
    end
    if n_diff == 0
        lignes{end+1} = sprintf('Réglages des assemblées : les mêmes pour les %d sessions.', numel(a_lus));
    else
        lignes{end+1} = sprintf(['ATTENTION : réglages des assemblées différents pour %d session(s) sur %d ' ...
            '(%s) : refaire l''étape 1 de ces sessions.'], n_diff, numel(a_lus), ...
            strjoin(unique(differents), ', '));
    end
end
fprintf('\n%s\n', strjoin(lignes, sprintf('\n')));
end


function t = texte_regle_sce(opts)
regle = 'percentile';
if isfield(opts, 'sce_regle'), regle = opts.sce_regle; end
switch regle
    case 'percentile'
        pc = NaN; if isfield(opts, 'sce_percentile'), pc = opts.sce_percentile; end
        t = sprintf('percentile %g, minimum %g', pc, opts.sce_n_cells_threshold);
    case 'fixe'
        t = sprintf('fixe %g', opts.sce_n_cells_threshold);
    case 'fdr'
        t = sprintf('fdr %g %%, décalage %g s', 100 * opts.sce_fdr, opts.sce_decalage_fdr_s);
    otherwise
        t = regle;
end
if isfield(opts, 'pretraitement_ancien') && opts.pretraitement_ancien
    t = [t ', ancien prétraitement'];
end
end


function L = ecrire_resume(bilans, noms_variantes, config, fichier, lignes_etape1, lignes_controles)
% Les deux phrases que les résultats permettent d'écrire, avec les tests.
if nargin < 6, lignes_controles = {}; end
L = {};
L{end+1} = '=================================================================';
L{end+1} = 'RÉSUMÉ : ce que les résultats permettent de dire';
L{end+1} = '=================================================================';
b = bilans{1};
L{end+1} = sprintf('Test principal : %s (référence des couleurs / permutation)', ...
    strrep(noms_variantes{1}, '_', ' / '));
L{end+1} = '';
L{end+1} = '1. Y a-t-il un biais clonal des assemblées en général ?';
L{end+1} = ['   « ' phrase_1(b, config) ' »'];
L{end+1} = sprintf(['   Robustesse : excès de %s points (IC 95 %% en rééchantillonnant les animaux : ' ...
    '%s à %s) ; en retirant chaque animal tour à tour, l''excès va de %s à %s points.'], ...
    nb(100 * (b.n_sig - b.n_attendu) / max(b.n_testables, 1), '%.1f'), nb(b.exces_ic(1), '%.1f'), ...
    nb(b.exces_ic(2), '%.1f'), nb(b.exces_sans_animal(1), '%.1f'), nb(b.exces_sans_animal(2), '%.1f'));
L{end+1} = sprintf(['   (Wilcoxon : avec %d animaux, la plus petite p possible est %s ; ' ...
    'le GLMM utilise toutes les assemblées.)'], b.n_animaux, nb(b.p_wilcoxon_min, '%.3f'));
L = [L, lignes_toutes_p(b)];
L{end+1} = '';
L{end+1} = '2. Lesquelles ?';
L{end+1} = ['   « ' phrase_2(b) ' »'];
if isfield(b, 'n_filtre')
    L{end+1} = sprintf(['   Sensibilité : FDR %g %% sur les seules assemblées qui pouvaient gagner (meilleur p possible ' ...
        '< %g : %d assemblées) : %d biaisée(s), dans %d animal(aux), dont %d absente(s) de la liste ci-dessus.'], ...
        100 * b.fdr, b.filtre_pmin, b.n_filtre, b.n_biaisees_filtre, b.n_animaux_filtre, b.n_nouvelles_filtre);
end
for v = 2:numel(bilans)
    c = bilans{v};
    communes = numel(intersect(c.ids_biaisees, b.ids_biaisees));
    L{end+1} = ''; %#ok<AGROW>
    L{end+1} = sprintf('Contrôle %s :', strrep(noms_variantes{v}, '_', ' / ')); %#ok<AGROW>
    L{end+1} = ['   1. ' phrase_1(c, config)]; %#ok<AGROW>
    L = [L, lignes_toutes_p(c)]; %#ok<AGROW>
    L{end+1} = sprintf('   2. %s Dont %d aussi biaisée(s) au test principal (sur %d).', ...
        phrase_2(c), communes, b.n_biaisees); %#ok<AGROW>
end
if ~isempty(lignes_controles)
    L{end+1} = '';
    L = [L, lignes_controles(:)'];
end
L{end+1} = '';
L{end+1} = 'Étape 1 :';
for k = 1:numel(lignes_etape1)
    L{end+1} = ['   ' lignes_etape1{k}]; %#ok<AGROW>
end
L{end+1} = '';
L{end+1} = ['Détail des tests : console de l''étape 4 (A1 Wilcoxon, A2 GLMM, A3 IC, A4 sans chaque ' ...
    'animal, C FDR) et dossiers agregation\<variante>\.'];
fid = fopen(fichier, 'w', 'n', 'UTF-8');
if fid < 0
    warning('Résumé non écrit : %s', fichier);
    return
end
fprintf(fid, '%s\n', L{:});
fclose(fid);
fprintf('\nRésumé écrit : %s (affiché à la fin de l''analyse)\n', fichier);
end


function t = phrase_1(b, config)
if b.n_testables == 0
    t = 'Aucune assemblée testable.';
    return
end
t = sprintf(['%s %% des assemblées testables (%d / %d) contiennent une couleur en excès ' ...
    '(p_minp < %s), contre %s %% attendu par hasard (n = %d animaux ; Wilcoxon signé sur les ' ...
    'animaux : %s ; GLMM binomial, animal et session en effets aléatoires : rapport des cotes ' ...
    '%s [IC 95 %% %s à %s], %s).'], ...
    nb(100 * b.n_sig / b.n_testables, '%.1f'), b.n_sig, b.n_testables, nb(config.clonal.seuilstat, '%g'), ...
    nb(100 * b.n_attendu / b.n_testables, '%.1f'), b.n_animaux, pv(b.p_wilcoxon), ...
    nb(b.OR, '%.2f'), nb(b.OR_ic(1), '%.2f'), nb(b.OR_ic(2), '%.2f'), pv(b.p_glmm));
end


function L = lignes_toutes_p(b)
% A5 (toutes les p) et A6 (proportion d'assemblées clonales), si calculés
L = {};
if b.n_score == 0, return; end
L{end+1} = sprintf(['   Avec toutes les p (score continu, 0 sans effet) : score moyen %s (Wilcoxon signé ' ...
    'sur les animaux : %s ; LMM, animal et session en effets aléatoires, ddl de Satterthwaite : %s [IC 95 %% %s à %s], %s).'], ...
    nb(b.score_moyen, '%.2f'), pv(b.p_wil_score), nb(b.lmm_score(1), '%.2f'), nb(b.lmm_score(2), '%.2f'), ...
    nb(b.lmm_score(3), '%.2f'), pv(b.p_lmm_score));
L{end+1} = sprintf(['   Proportion estimée d''assemblées clonales : %s %% (IC 95 %% en rééchantillonnant ' ...
    'les animaux : %s à %s %%) ; c''est un minimum.'], nb(100 * b.pi1, '%.1f'), ...
    nb(100 * b.pi1_ic(1), '%.1f'), nb(100 * b.pi1_ic(2), '%.1f'));
end


function t = phrase_2(b)
if b.n_testables == 0
    t = 'Aucune assemblée testable.';
elseif b.n_biaisees == 0
    t = sprintf(['Aucune assemblée biaisée après correction (Benjamini-Hochberg sur les %d ' ...
        'assemblées testables, FDR %g %%).'], b.n_testables, 100 * b.fdr);
else
    t = sprintf(['%d assemblée(s) biaisée(s) (Benjamini-Hochberg sur les %d assemblées testables, ' ...
        'FDR %g %%), dans %d session(s) de %d animal(aux) ; couleur en excès : %s.'], ...
        b.n_biaisees, b.n_testables, 100 * b.fdr, b.n_sessions_biaisees, b.n_animaux_biaisees, ...
        b.couleurs_biaisees);
end
end


function t = nb(x, fmt)
% nombre en texte, virgule décimale
if isempty(x) || isnan(x), t = 'NA'; return; end
t = strrep(sprintf(fmt, x), '.', ',');
end


function t = pv(p)
% « p = 0,034 » ou « p < 0,001 »
if isempty(p) || isnan(p), t = 'p = NA';
elseif p < 0.001, t = 'p < 0,001';
else, t = ['p = ' nb(p, '%.3f')];
end
end
