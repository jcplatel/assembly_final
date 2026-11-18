function T = agregation_hierarchique(sessions, config)
%AGREGATION_HIERARCHIQUE  Étape 3c, toutes les sessions : jusqu'à quelle
%   taille les groupes de l'arbre hiérarchique ont-ils une couleur en excès ?
%   Par animal, pour chaque arbre, chaque hasard (champ, blocs) et chaque
%   coupe (groupes de 5, 10, 20, 40 cellules en moyenne ; 0 = clusters FCA) :
%   groupes avec une couleur en excès (somme des sessions) / groupes testés,
%   observé et attendu (couleurs mélangées), et groupes significatifs après
%   FDR. Test par coupe : Wilcoxon unilatéral sur les animaux (observé -
%   attendu > 0) : 4 tests par arbre, sans correction.
%   Écrit dans <PathSave>\hierarchique_agregation\ : hierarchique_sessions.csv,
%   hierarchique_animaux.csv, hierarchique_agregation.mat,
%   resume_hierarchique.txt, hierarchique.png.

opts = options_hierarchique(config);
dossier = fullfile(config.PathSave, 'hierarchique_agregation');
if ~isfolder(dossier), mkdir(dossier); end
morceaux = {};
for s = 1:height(sessions)
    f = fullfile(sessions.dossier{s}, 'hierarchique', 'hierarchique.mat');
    if ~isfile(f), continue; end
    H = load(f, 'T');
    if isfield(H, 'T') && ~isempty(H.T)
        morceaux{end+1, 1} = H.T; %#ok<AGROW>
    end
end
T = table();
if isempty(morceaux)
    fprintf('Hiérarchique : aucune session avec des résultats (lancer main_hierarchique).\n');
    return
end
T = vertcat(morceaux{:});
writetable(T, fullfile(dossier, 'hierarchique_sessions.csv'));

% par animal : somme des sessions
[G, arbre, perm, niveau, animal] = findgroups(T.arbre, T.permutation, T.niveau, T.animal);
A = table(arbre, perm, niveau, animal, splitapply(@sum, T.n_groupes, G), splitapply(@sum, T.observes, G), ...
    splitapply(@sum, T.attendus, G), splitapply(@sum, T.n_fdr, G), splitapply(@numel, T.n_groupes, G), ...
    'VariableNames', {'arbre', 'permutation', 'niveau', 'animal', 'n_groupes', 'observes', 'attendus', ...
    'n_fdr', 'n_sessions'});
A = A(A.n_groupes > 0, :);
A.pct_observe = 100 * A.observes ./ A.n_groupes;
A.pct_attendu = 100 * A.attendus ./ A.n_groupes;
writetable(A, fullfile(dossier, 'hierarchique_animaux.csv'));

% entre animaux, par coupe
[G2, arbre, perm, niveau] = findgroups(A.arbre, A.permutation, A.niveau);
B = table(arbre, perm, niveau, splitapply(@mean, A.pct_observe, G2), ...
    splitapply(@(x) std(x) / sqrt(numel(x)), A.pct_observe, G2), splitapply(@mean, A.pct_attendu, G2), ...
    splitapply(@numel, A.pct_observe, G2), splitapply(@sum, A.n_groupes, G2), splitapply(@sum, A.n_fdr, G2), ...
    splitapply(@(o, e) wilcoxon_droite(o - e), A.pct_observe, A.pct_attendu, G2), ...
    'VariableNames', {'arbre', 'permutation', 'niveau', 'pct_observe', 'pct_sem', ...
    'pct_attendu', 'n_animaux', 'n_groupes', 'n_fdr', 'p_wilcoxon'});
B = sortrows(B, {'arbre', 'permutation', 'niveau'});
L = resume(B, opts);
fid = fopen(fullfile(dossier, 'resume_hierarchique.txt'), 'w', 'n', 'UTF-8');
if fid >= 0
    fprintf(fid, '%s\n', L{:});
    fclose(fid);
end
fprintf('\n%s\n', L{:});
lignes_resume = L; %#ok<NASGU>
date_hier = datestr(now, 'yyyy-mm-dd HH:MM'); %#ok<NASGU>
save(fullfile(dossier, 'hierarchique_agregation.mat'), 'T', 'A', 'B', 'opts', 'lignes_resume', 'date_hier');
try
    figure_hierarchique(dossier);
catch err
    warning('figure_hierarchique non exécutée : %s', err.message);
end
end


function p = wilcoxon_droite(d)
d = d(isfinite(d));
p = NaN;
if numel(d) >= 2 && any(d ~= 0)
    p = signrank(d, 0, 'tail', 'right');
end
end


function L = resume(B, opts)
noms_arbre = struct('sce', 'co-participation aux SCE', 'sttc', 'STTC résiduel au repos (+/- 200 ms)', ...
    'sttc_100ms', 'STTC résiduel au repos (+/- 100 ms)', 'dff', 'Pearson du dF/F résiduel au repos', ...
    'fca', 'arbre FCA (Feldt et al. 2009)', 'fca_amas', 'clusters FCA seuls');
noms_perm = struct('globale', 'couleurs mélangées dans le champ', 'blocs', ...
    sprintf('couleurs mélangées par blocs de %g µm', opts.taille_bloc_um));
L = {};
L{end+1} = ['CLUSTERING HIÉRARCHIQUE (étape 3c, exploratoire) : arbre coupé en groupes de ~5, ~10, ~20, ~40 ' ...
    sprintf('cellules ; groupes avec une couleur en excès (min des p hypergéométriques < %g), observés ', opts.alpha) ...
    'contre attendus (couleurs mélangées), % des groupes testés, moyenne des animaux ; p = Wilcoxon unilatéral ' ...
    sprintf('sur les animaux, sans correction ; FDR = groupes significatifs (p_minP, FDR %g %%) :', 100 * opts.fdr_q)];
for a = unique(B.arbre)'
    for q = unique(B.permutation)'
        b = B(strcmp(B.arbre, a{1}) & strcmp(B.permutation, q{1}), :);
        t = arrayfun(@(i) sprintf('%s : %.0f %% contre %.0f %% (%s ; FDR %d/%d)', texte_niveau(b.niveau(i)), ...
            b.pct_observe(i), b.pct_attendu(i), pv(b.p_wilcoxon(i)), b.n_fdr(i), b.n_groupes(i)), 1:height(b), ...
            'UniformOutput', false);
        L{end+1} = sprintf('   %s, %s : %s', champ(noms_arbre, a{1}), champ(noms_perm, q{1}), strjoin(t, ' ; ')); %#ok<AGROW>
    end
end
end


function t = champ(s, nom)
if isfield(s, nom), t = s.(nom); else, t = nom; end
end


function t = texte_niveau(n)
if n == 0, t = 'clusters FCA'; else, t = sprintf('~%g cellules', n); end
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
