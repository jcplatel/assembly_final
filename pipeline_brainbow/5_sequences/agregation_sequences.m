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

dossier = fullfile(config.PathSave, 'sequences_agregation');
if ~isfolder(dossier), mkdir(dossier); end
statuts = {};
lignes = cell(0, 21);
champs_c = {'coul_n_colorees', 'coul_n_meme', 'coul_n_diff', 'coul_d_meme', 'coul_d_diff', 'coul_delta', ...
    'coul_pct', 'coul_p', 'coul_rho_espace'};
for s = 1:height(sessions)
    f = fullfile(sessions.dossier{s}, 'sequences', 'sequences.mat');
    if ~isfile(f), continue; end
    vars = intersect([{'NCell', 'n_courses', 'duree_course_s', 'n_sequence', 'n_sequence_tout', 'n_deux', ...
        'n_deux_hasard', 'rho', 'p_rho', 'statut'}, champs_c], who('-file', f));
    R = load(f, vars{:});
    for c = [{'n_sequence', 'n_sequence_tout', 'n_deux', 'n_deux_hasard', 'rho', 'p_rho'}, champs_c]
        if ~isfield(R, c{1}), R.(c{1}) = NaN; end     % trop peu de courses, ou ancienne session sans couleur
    end
    lignes(end+1, :) = {sessions.id_animal{s}, sessions.nom{s}, R.NCell, R.n_courses, R.duree_course_s, ...
        R.n_sequence, 100 * R.n_sequence / R.NCell, R.n_sequence_tout, R.n_deux, R.n_deux_hasard, ...
        R.rho, R.p_rho, R.coul_n_colorees, R.coul_n_meme, R.coul_n_diff, R.coul_d_meme, R.coul_d_diff, ...
        R.coul_delta, R.coul_pct, R.coul_p, R.coul_rho_espace}; %#ok<AGROW>
    statuts{size(lignes, 1), 1} = R.statut; %#ok<AGROW>
end
T = table();
if isempty(lignes)
    fprintf('Séquences : aucune session avec des résultats (lancer l''étape 5).\n');
    return
end
T = cell2table(lignes, 'VariableNames', {'animal', 'session', 'n_cellules', 'n_courses', 'duree_course_s', ...
    'n_sequence', 'pct_sequence', 'n_sequence_tout', 'n_retrouvees', 'n_retrouvees_hasard', 'rho', ...
    'p_rho', 'n_colorees_seq', 'n_paires_meme', 'n_paires_diff', 'dist_meme', 'dist_diff', 'delta_couleur', ...
    'pct_couleur', 'p_couleur', 'rho_espace'});
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
L = [resume(T, bilan, options_sequences(config)), resume_couleur(bc)];
fid = fopen(fullfile(dossier, 'resume_sequences.txt'), 'w', 'n', 'UTF-8');
if fid >= 0
    fprintf(fid, '%s\n', L{:});
    fclose(fid);
end
fprintf('\n%s\n', L{:});
date_sequences = datestr(now, 'yyyy-mm-dd HH:MM'); %#ok<NASGU>
lignes_resume = L; %#ok<NASGU>
save(fullfile(dossier, 'sequences_agregation.mat'), 'T', 'A', 'bilan', 'Ac', 'bc', 'lignes_resume', ...
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


function L = resume(T, b, opts)
L = {};
L{end+1} = sprintf(['SÉQUENCES PENDANT LA COURSE (étape 5, sans couleur ; méthode de sequence_2026_01_27 : ' ...
    'corrélation croisée avec les débuts de course, seuil d''Otsu, décalage < 0) :']);
L{end+1} = sprintf('   %d sessions sur %d analysées (au moins %d courses) ; %d animaux.', ...
    b.n_sessions, b.n_sessions_total, opts.min_courses, b.n_animaux);
if b.n_sessions == 0, return; end
Tok = T(isfinite(T.rho), :);
L{end+1} = sprintf(['   Cellules de la séquence : %s %% des cellules en médiane [%s-%s] ; courses par session : ' ...
    'médiane %s [%s-%s].'], nb(median(Tok.pct_sequence), '%.0f'), nb(min(Tok.pct_sequence), '%.0f'), ...
    nb(max(Tok.pct_sequence), '%.0f'), nb(median(Tok.n_courses), '%.0f'), nb(min(Tok.n_courses), '%.0f'), ...
    nb(max(Tok.n_courses), '%.0f'));
L{end+1} = sprintf(['   « Rangées sur la moitié des courses, les cellules gardent leur ordre sur l''autre ' ...
    'moitié : rho de Spearman moyen des animaux = %s ± %s (sem), %d/%d animaux > 0, Wilcoxon unilatéral sur ' ...
    'les animaux, %s ; %d/%d sessions à p < 0,05. »'], nb(b.rho_moyen, '%.2f'), nb(b.rho_sem, '%.2f'), ...
    b.n_pos, b.n_animaux, pv(b.p_wilcoxon), b.n_sig, b.n_sessions);
L{end+1} = sprintf(['   Cellules retrouvées dans la séquence de l''autre moitié : %d (%d attendues si les deux ' ...
    'choix étaient indépendants).'], sum(Tok.n_retrouvees), round(sum(Tok.n_retrouvees_hasard)));
end


function L = resume_couleur(b)
L = {''};
L{end+1} = ['MÊME COULEUR PLUS PROCHE DANS LA SÉQUENCE ? (ordre de toutes les courses ; distance = écart de ' ...
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
