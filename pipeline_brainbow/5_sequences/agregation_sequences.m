function T = agregation_sequences(sessions, config)
%AGREGATION_SEQUENCES  Étape 5, toutes les sessions.
%   1. Séquences reproductibles pendant la course ? Par session : rho de
%      Spearman entre le décalage des cellules de la séquence sur une moitié
%      des courses et sur l'autre. Par animal : moyenne de ses sessions.
%      Wilcoxon unilatéral sur les animaux (rho > 0).
%   2. Même couleur plus proche dans la séquence ? (test de l'utilisateur)
%      Par session : distance moyenne dans la séquence des paires de même
%      couleur - couleurs différentes (delta, sans correction). Par animal :
%      moyenne de ses sessions. Wilcoxon unilatéral sur les animaux
%      (delta < 0). Contrôle : la distance dans la séquence suit-elle la
%      distance dans l'espace (rho_espace) ?
%   Écrit dans <PathSave>\sequences_agregation\ : sequences_sessions.csv,
%   sequences_agregation.mat, resume_sequences.txt, sequences.png,
%   sequence_colour.png.

opts_s = options_sequences(config);
dossier = fullfile(config.PathSave, strrep(dossier_sequences(opts_s), 'sequences', 'sequences_agregation'));   % 20 nov. : '_cascade' si signal CASCADE
if ~isfolder(dossier), mkdir(dossier); end
statuts = {};
validations = {};
lignes = cell(0, 22);
% 20 nov. : assemblées et séquence (sequence_assemblees)
champs_a = {'asm_n_assemblees', 'asm_n_enrichies', 'asm_frac_seq_assemblees', 'asm_frac_seq_toutes', ...
    'asm_n_cellules_paires', 'asm_n_meme', 'asm_delta', 'asm_pct', 'asm_p', 'n_robustes'};
lignes_a = zeros(0, numel(champs_a));
champs_c = {'coul_n_colorees', 'coul_n_meme', 'coul_n_diff', 'coul_d_meme', 'coul_d_diff', 'coul_delta', ...
    'coul_pct', 'coul_p', 'coul_rho_espace'};
for s = 1:height(sessions)
    f = fullfile(sessions.dossier{s}, dossier_sequences(opts_s), 'sequences.mat');
    if ~isfile(f), continue; end
    vars = intersect([{'NCell', 'n_courses', 'duree_course_s', 'n_sequence', 'n_sequence_tout', 'n_deux', ...
        'n_deux_hasard', 'rho', 'p_rho', 'statut', 'validation', 'n_validees'}, champs_c], who('-file', f));
    R = load(f, vars{:});
    for c = [{'n_sequence', 'n_sequence_tout', 'n_deux', 'n_deux_hasard', 'rho', 'p_rho', 'n_validees'}, champs_c]
        if ~isfield(R, c{1}), R.(c{1}) = NaN; end     % trop peu de courses, ou ancienne session sans couleur
    end
    if ~isfield(R, 'validation'), R.validation = 'moitie (1 coupure 50/50)'; end
    validations{end+1, 1} = R.validation; %#ok<AGROW>
    lignes(end+1, :) = {sessions.id_animal{s}, sessions.nom{s}, R.NCell, R.n_courses, R.duree_course_s, ...
        R.n_sequence, 100 * R.n_sequence / R.NCell, R.n_sequence_tout, R.n_deux, R.n_deux_hasard, ...
        R.rho, R.p_rho, R.coul_n_colorees, R.coul_n_meme, R.coul_n_diff, R.coul_d_meme, R.coul_d_diff, ...
        R.coul_delta, R.coul_pct, R.coul_p, R.coul_rho_espace, R.n_validees}; %#ok<AGROW>
    statuts{size(lignes, 1), 1} = R.statut; %#ok<AGROW>
    va = intersect(champs_a, who('-file', f));
    Ra = struct();
    if ~isempty(va), Ra = load(f, va{:}); end
    la = nan(1, numel(champs_a));
    for c = 1:numel(champs_a)
        if isfield(Ra, champs_a{c}), la(c) = double(Ra.(champs_a{c})); end
    end
    lignes_a(end+1, :) = la; %#ok<AGROW>
end
T = table();
if isempty(lignes)
    fprintf('Séquences : aucune session avec des résultats (lancer l''étape 5).\n');
    return
end
T = cell2table(lignes, 'VariableNames', {'animal', 'session', 'n_cellules', 'n_courses', 'duree_course_s', ...
    'n_sequence', 'pct_sequence', 'n_sequence_tout', 'n_retrouvees', 'n_retrouvees_hasard', 'rho', ...
    'p_rho', 'n_colorees_seq', 'n_paires_meme', 'n_paires_diff', 'dist_meme', 'dist_diff', 'delta_couleur', ...
    'pct_couleur', 'p_couleur', 'rho_espace', 'n_validees'});
T.validation = validations;
for c = 1:numel(champs_a)
    T.(champs_a{c}) = lignes_a(:, c);
end
T.statut = statuts;
writetable(T, fullfile(dossier, 'sequences_sessions.csv'));

% Par animal, puis entre animaux
Tok = T(isfinite(T.rho), :);
bilan = struct('n_sessions', height(Tok), 'n_sessions_total', height(T), 'n_sig', sum(Tok.p_rho < 0.05), ...
    'n_animaux', 0, 'rho_moyen', NaN, 'rho_sem', NaN, 'p_wilcoxon', NaN, 'n_pos', 0);
A = table();
if ~isempty(Tok)
    [G, anim] = findgroups(Tok.animal);
    A = table(anim, splitapply(@mean, Tok.rho, G), splitapply(@numel, Tok.rho, G), ...
        splitapply(@mean, Tok.pct_sequence, G), 'VariableNames', {'animal', 'rho', 'n_sessions', 'pct_sequence'});
    bilan.n_animaux = height(A);
    bilan.rho_moyen = mean(A.rho);
    bilan.rho_sem = std(A.rho) / sqrt(height(A));
    bilan.n_pos = sum(A.rho > 0);
    if height(A) >= 2
        bilan.p_wilcoxon = signrank(A.rho, 0, 'tail', 'right');
    end
end
% Couleurs : même couleur plus proche dans la séquence ?
Tc = T(isfinite(T.delta_couleur), :);
bc = struct('n_sessions', height(Tc), 'n_sig', sum(Tc.p_couleur < 0.05), 'n_animaux', 0, ...
    'delta_moyen', NaN, 'delta_sem', NaN, 'pct_moyen', NaN, 'p_wilcoxon', NaN, 'n_neg', 0, ...
    'rho_espace_moyen', NaN, 'p_espace', NaN);
Ac = table();
if ~isempty(Tc)
    [G, anim] = findgroups(Tc.animal);
    Ac = table(anim, splitapply(@mean, Tc.dist_meme, G), splitapply(@mean, Tc.dist_diff, G), ...
        splitapply(@mean, Tc.delta_couleur, G), splitapply(@mean, Tc.pct_couleur, G), ...
        splitapply(@(x) mean(x, 'omitnan'), Tc.rho_espace, G), splitapply(@numel, Tc.delta_couleur, G), ...
        'VariableNames', {'animal', 'dist_meme', 'dist_diff', 'delta', 'pct', 'rho_espace', 'n_sessions'});
    bc.n_animaux = height(Ac);
    bc.delta_moyen = mean(Ac.delta);
    bc.delta_sem = std(Ac.delta) / sqrt(height(Ac));
    bc.pct_moyen = mean(Ac.pct);
    bc.n_neg = sum(Ac.delta < 0);
    bc.rho_espace_moyen = mean(Ac.rho_espace, 'omitnan');
    if height(Ac) >= 2
        bc.p_wilcoxon = signrank(Ac.delta, 0, 'tail', 'left');
        re = Ac.rho_espace(isfinite(Ac.rho_espace));
        if numel(re) >= 2, bc.p_espace = signrank(re, 0); end
    end
end
% Assemblées (20 nov.) : enrichies en cellules de la séquence ? regroupées dans la séquence ?
Ta = T(isfinite(T.asm_n_assemblees) & T.asm_n_assemblees > 0, :);
ba = struct('n_sessions', height(Ta), 'n_assemblees', sum(Ta.asm_n_assemblees), ...
    'n_enrichies', sum(Ta.asm_n_enrichies), 'n_animaux', 0, 'frac_asm', NaN, 'frac_toutes', NaN, ...
    'p_enrichissement', NaN, 'n_plus', 0, 'n_sessions_reg', 0, 'n_sig_reg', 0, 'n_animaux_reg', 0, ...
    'pct_moyen', NaN, 'p_regroupement', NaN, 'n_neg', 0);
Aa = table();
if ~isempty(Ta)
    [G, anim] = findgroups(Ta.animal);
    Aa = table(anim, splitapply(@(x) mean(x, 'omitnan'), Ta.asm_frac_seq_assemblees, G), ...
        splitapply(@mean, Ta.asm_frac_seq_toutes, G), splitapply(@numel, Ta.asm_frac_seq_toutes, G), ...
        'VariableNames', {'animal', 'frac_asm', 'frac_toutes', 'n_sessions'});
    ba.n_animaux = height(Aa);
    ba.frac_asm = mean(Aa.frac_asm, 'omitnan');
    ba.frac_toutes = mean(Aa.frac_toutes);
    ba.n_plus = sum(Aa.frac_asm > Aa.frac_toutes);
    if height(Aa) >= 2
        ba.p_enrichissement = signrank(Aa.frac_asm, Aa.frac_toutes, 'tail', 'right');
    end
    Tr = Ta(isfinite(Ta.asm_delta), :);
    ba.n_sessions_reg = height(Tr);
    ba.n_sig_reg = sum(Tr.asm_p < 0.05);
    if ~isempty(Tr)
        [G, anim] = findgroups(Tr.animal);
        Ar = table(anim, splitapply(@mean, Tr.asm_delta, G), splitapply(@mean, Tr.asm_pct, G), ...
            'VariableNames', {'animal', 'delta', 'pct'});
        ba.n_animaux_reg = height(Ar);
        ba.pct_moyen = mean(Ar.pct);
        ba.n_neg = sum(Ar.delta < 0);
        if height(Ar) >= 2
            ba.p_regroupement = signrank(Ar.delta, 0, 'tail', 'left');
        end
        Aa = outerjoin(Aa, Ar, 'Keys', 'animal', 'MergeKeys', true);
    end
end
kfold = any(startsWith(T.validation, {'5-fold', '2-fold'}));
L = [resume(T, bilan, opts_s, kfold), resume_couleur(bc, kfold, opts_s), resume_assemblees(ba)];
if strcmp(opts_s.signal, 'cascade')
    L = [{sprintf('SIGNAL : probabilité de spike de CASCADE (%s), pas le calcium.', opts_s.cascade_modele), ''}, L];
end
fid = fopen(fullfile(dossier, 'resume_sequences.txt'), 'w', 'n', 'UTF-8');
if fid >= 0
    fprintf(fid, '%s\n', L{:});
    fclose(fid);
end
fprintf('\n%s\n', L{:});
date_sequences = datestr(now, 'yyyy-mm-dd HH:MM'); %#ok<NASGU>
lignes_resume = L; %#ok<NASGU>
save(fullfile(dossier, 'sequences_agregation.mat'), 'T', 'A', 'bilan', 'Ac', 'bc', 'Aa', 'ba', 'lignes_resume', ...
    'date_sequences');
try
    figure_sequences(dossier);
catch err
    warning('figure_sequences non exécutée : %s', err.message);
end
try
    figure_sequence_couleur(dossier);
catch err
    warning('figure_sequence_couleur non exécutée : %s', err.message);
end
end


function L = resume(T, b, opts, kfold)
if nargin < 4, kfold = false; end
L = {};
if strcmp(opts.selection, 'faux_debuts')
    L{end+1} = sprintf(['SÉQUENCES PENDANT LA COURSE (étape 5, sans couleur ; corrélation croisée avec les débuts ' ...
        'de course sur le signal continu, décalage jusqu''à %g s ; cellule retenue si son pic dépasse celui de ' ...
        '%d faux débuts tirés pendant la course, FDR %g %%, et décalage < 0) :'], opts.max_lag_s, ...
        opts.n_faux_debuts, 100 * opts.fdr_sequence);
else
    L{end+1} = sprintf(['SÉQUENCES PENDANT LA COURSE (étape 5, sans couleur ; méthode de sequence_2026_01_27 : ' ...
        'corrélation croisée avec les débuts de course, seuil d''Otsu, décalage < 0) :']);
end
L{end+1} = sprintf('   %d sessions sur %d analysées (au moins %d courses) ; %d animaux.', ...
    b.n_sessions, b.n_sessions_total, opts.min_courses, b.n_animaux);
if b.n_sessions == 0, return; end
Tok = T(isfinite(T.rho), :);
L{end+1} = sprintf(['   Cellules de la séquence : %s %% des cellules en médiane [%s-%s] ; courses par session : ' ...
    'médiane %s [%s-%s].'], nb(median(Tok.pct_sequence), '%.0f'), nb(min(Tok.pct_sequence), '%.0f'), ...
    nb(max(Tok.pct_sequence), '%.0f'), nb(median(Tok.n_courses), '%.0f'), nb(min(Tok.n_courses), '%.0f'), ...
    nb(max(Tok.n_courses), '%.0f'));
if kfold
    L{end+1} = sprintf(['   Validation croisée répétée : 5-fold x 20 (100 coupures 80 / 20) ; sessions de moins de ' ...
        '%d courses : 50 / 50 x 50. rho d''une session = moyenne des 100 coupures, p par permutation.'], ...
        opts.kfold_min_courses);
    L{end+1} = sprintf(['   « Rangées sur les courses d''apprentissage, les cellules gardent leur ordre sur les ' ...
        'courses de test : rho de Spearman moyen des animaux = %s ± %s (sem), %d/%d animaux > 0, Wilcoxon ' ...
        'unilatéral sur les animaux, %s ; %d/%d sessions à p < 0,05. »'], nb(b.rho_moyen, '%.2f'), ...
        nb(b.rho_sem, '%.2f'), b.n_pos, b.n_animaux, pv(b.p_wilcoxon), b.n_sig, b.n_sessions);
    L{end+1} = sprintf(['   Cellules retrouvées dans la séquence des courses de test (moyenne des coupures, somme ' ...
        'des sessions) : %.0f (%.0f attendues si les deux choix étaient indépendants) ; cellules validées ' ...
        '(retenues et retrouvées dans >= %.0f %% des coupures) : %d.'], sum(Tok.n_retrouvees), ...
        sum(Tok.n_retrouvees_hasard), 100 * opts.kfold_seuil_valide, sum(Tok.n_validees, 'omitnan'));
    if strcmp(opts.cellules_tests, 'robustes')
        L{end+1} = sprintf(['   Cellules robustes (validées ET position stable entre apprentissage et test, FDR ' ...
            '%g %%) : %d au total, médiane %s par session [%s-%s].'], 100 * opts.stabilite_fdr, ...
            sum(Tok.n_robustes, 'omitnan'), nb(median(Tok.n_robustes, 'omitnan'), '%.0f'), ...
            nb(min(Tok.n_robustes), '%.0f'), nb(max(Tok.n_robustes), '%.0f'));
    end
    return
end
L{end+1} = sprintf(['   « Rangées sur la moitié des courses, les cellules gardent leur ordre sur l''autre ' ...
    'moitié : rho de Spearman moyen des animaux = %s ± %s (sem), %d/%d animaux > 0, Wilcoxon unilatéral sur ' ...
    'les animaux, %s ; %d/%d sessions à p < 0,05. »'], nb(b.rho_moyen, '%.2f'), nb(b.rho_sem, '%.2f'), ...
    b.n_pos, b.n_animaux, pv(b.p_wilcoxon), b.n_sig, b.n_sessions);
L{end+1} = sprintf(['   Cellules retrouvées dans la séquence de l''autre moitié : %d (%d attendues si les deux ' ...
    'choix étaient indépendants).'], sum(Tok.n_retrouvees), round(sum(Tok.n_retrouvees_hasard)));
end


function L = resume_couleur(b, kfold, opts)
if nargin < 2, kfold = false; end
L = {''};
if kfold && nargin >= 3 && strcmp(opts.cellules_tests, 'robustes')
    quel_ordre = ['cellules robustes : validées par la validation croisée ET de position stable entre ' ...
        'apprentissage et test (FDR 5 %), rangées par leur décalage moyen sur les courses de test'];
elseif kfold
    quel_ordre = 'cellules validées par la validation croisée, rangées par leur décalage moyen sur les courses de test';
else
    quel_ordre = 'ordre de toutes les courses';
end
L{end+1} = ['MÊME COULEUR PLUS PROCHE DANS LA SÉQUENCE ? (' quel_ordre ' ; distance = écart de ' ...
    'position dans la séquence, en fraction de sa longueur ; paires de cellules colorées de la séquence, sans ' ...
    'correction) :'];
if b.n_sessions == 0
    L{end+1} = '   aucune session testable.';
    return
end
L{end+1} = sprintf(['   « Les paires de même couleur sont %s %% %s dans la séquence que les paires de couleurs ' ...
    'différentes (moyenne des animaux), %d/%d animaux dans ce sens, Wilcoxon unilatéral sur les animaux %s ; ' ...
    '%d/%d sessions à p < 0,05 (couleurs mélangées). »'], nb(abs(b.pct_moyen), '%.1f'), ...
    texte_sens(b.pct_moyen), b.n_neg, b.n_animaux, pv(b.p_wilcoxon), b.n_sig, b.n_sessions);
L{end+1} = sprintf(['   Contrôle : distance dans la séquence contre distance dans l''espace, rho de Spearman moyen ' ...
    'des animaux = %s (Wilcoxon bilatéral contre 0 : %s) ; proche de 0 = la séquence ne suit pas l''espace.'], ...
    nb(b.rho_espace_moyen, '%.3f'), pv(b.p_espace));
end


function L = resume_assemblees(b)
L = {''};
L{end+1} = ['ASSEMBLÉES (étape 2, SCE au repos) ET SÉQUENCE (mêmes cellules et même ordre que le test de ' ...
    'couleur ; « dans la séquence » = parmi ces cellules) :'];
if b.n_sessions == 0
    L{end+1} = '   aucune session avec des assemblées.';
    return
end
L{end+1} = sprintf(['   « Part des cellules d''assemblées dans la séquence : %s %% contre %s %% de toutes les ' ...
    'cellules (moyenne des animaux), %d/%d animaux au-dessus, Wilcoxon unilatéral apparié sur les animaux %s ; ' ...
    '%d/%d assemblées enrichies en cellules de la séquence (hypergéométrique, FDR 5 %% par session), %d ' ...
    'sessions. »'], nb(100 * b.frac_asm, '%.0f'), nb(100 * b.frac_toutes, '%.0f'), b.n_plus, b.n_animaux, ...
    pv(b.p_enrichissement), b.n_enrichies, b.n_assemblees, b.n_sessions);
if b.n_sessions_reg == 0
    L{end+1} = '   Regroupement : aucune session testable.';
    return
end
L{end+1} = sprintf(['   « Les paires de cellules de la même assemblée sont %s %% %s dans la séquence que les ' ...
    'paires d''assemblées différentes (moyenne des animaux), %d/%d animaux dans ce sens, Wilcoxon unilatéral ' ...
    'sur les animaux %s ; %d/%d sessions à p < 0,05 (assemblées mélangées). »'], nb(abs(b.pct_moyen), '%.1f'), ...
    texte_sens(b.pct_moyen), b.n_neg, b.n_animaux_reg, pv(b.p_regroupement), b.n_sig_reg, b.n_sessions_reg);
end


function t = texte_sens(x)
if x < 0, t = 'plus proches'; else, t = 'plus éloignées'; end
end


function t = nb(x, fmt)
if isnan(x), t = 'NA'; else, t = sprintf(fmt, x); end
end


function t = pv(p)
if isnan(p)
    t = 'p = NA';
elseif p < 0.001
    t = 'p < 0,001';
else
    t = ['p = ' strrep(sprintf('%.3f', p), '.', ',')];
end
end
