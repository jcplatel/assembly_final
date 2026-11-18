function C = agregation_coeurs(sessions, config)
%AGREGATION_COEURS  Étape 3d : verdict « cœur clonal » de chaque assemblée,
%   toutes les sessions (coeurs_assemblees).
%   FDR de Benjamini-Hochberg (opts.fdr_q) sur toutes les assemblées
%   testables, séparément pour p_champ et p_voisines. Verdict (revu après
%   simulation, 7 nov. : une assemblée compacte sans aucun choix clonal
%   passe p_champ une fois sur deux, p_voisines 8 % du temps ; un cœur de 4
%   sœurs dans 20 cellules passe p_voisines 28 % du temps) :
%     'pas testable'       : moins de 3 membres colorés, ou aucune paire de
%                            membres colorés assez proche ;
%     'cœur clonal'        : p_voisines significatif après FDR : plus de
%                            sœurs que des voisines n'en donneraient ;
%     'contient des sœurs' : seulement p_champ significatif (peut venir de la
%                            seule compacité de l'assemblée) ;
%     'pas de preuve'      : sinon (pas « pas clonale » : beaucoup de cœurs
%                            réels restent invisibles).
%   Le score continu z_voisines est à préférer pour relier le cœur clonal à
%   une autre propriété des assemblées (séquences…).
%   Écrit dans <PathSave>\soeurs_agregation\ : coeurs_assemblees.csv,
%   coeurs_agregation.mat, resume_coeurs.txt, coeurs_assemblees.png.

opts = options_soeurs(config);
dossier = fullfile(config.PathSave, 'soeurs_agregation');
if ~isfolder(dossier), mkdir(dossier); end
morceaux = {};
for s = 1:height(sessions)
    f = fullfile(sessions.dossier{s}, 'soeurs', 'coeurs.mat');
    if ~isfile(f), continue; end
    R = load(f, 'C');
    if ~isempty(R.C), morceaux{end+1, 1} = R.C; end %#ok<AGROW>
end
C = table();
if isempty(morceaux)
    fprintf('Cœurs clonaux : aucune assemblée (étapes 2-3 faites ? main_soeurs lancé ?).\n');
    return
end
C = vertcat(morceaux{:});
t = C.testable;
C.fdr_champ = false(height(C), 1);
C.fdr_voisines = false(height(C), 1);
if any(t)
    C.fdr_champ(t) = bh_fdr_mask(C.p_champ(t)', opts.fdr_q)';
    C.fdr_voisines(t) = bh_fdr_mask(C.p_voisines(t)', opts.fdr_q)';
end
C.verdict = repmat({'pas de preuve'}, height(C), 1);
C.verdict(~t) = {'pas testable'};
C.verdict(t & C.fdr_champ) = {'contient des sœurs'};
C.verdict(t & C.fdr_voisines) = {'cœur clonal'};
writetable(C, fullfile(dossier, 'coeurs_assemblees.csv'));

[G, anim] = findgroups(C.animal);
P = table(anim, splitapply(@numel, C.testable, G), splitapply(@sum, C.testable, G), ...
    splitapply(@(v) sum(strcmp(v, 'cœur clonal')), C.verdict, G), ...
    splitapply(@(v) sum(strcmp(v, 'contient des sœurs')), C.verdict, G), ...
    'VariableNames', {'animal', 'n_assemblees', 'n_testables', 'n_coeur', 'n_soeurs_seulement'});

L = {};
L{end+1} = sprintf(['CŒUR CLONAL DE CHAQUE ASSEMBLÉE (score = paires de sœurs attendues parmi ses membres ; ' ...
    'verdict : plus de sœurs que des fausses assemblées de voisines, FDR %g %% sur les assemblées ' ...
    'testables) :'], 100 * opts.fdr_q);
L{end+1} = sprintf(['   %d assemblées, %d testables : %d avec un cœur clonal (%.0f %% des testables) ; %d ' ...
    'contiennent des sœurs sans plus que leurs voisines ; %d sans preuve.'], height(C), sum(t), ...
    sum(P.n_coeur), 100 * sum(P.n_coeur) / max(1, sum(t)), sum(P.n_soeurs_seulement), ...
    sum(strcmp(C.verdict, 'pas de preuve')));
L{end+1} = ['   Par animal (cœur / testables) : ' strjoin(arrayfun(@(i) sprintf('%s %d/%d', P.animal{i}, ...
    P.n_coeur(i), P.n_testables(i)), 1:height(P), 'UniformOutput', false), ' ; ')];
L{end+1} = '   « Pas de preuve » ne veut pas dire « pas clonale » (voir coeurs_assemblees.csv, score z_voisines).';
fid = fopen(fullfile(dossier, 'resume_coeurs.txt'), 'w', 'n', 'UTF-8');
if fid >= 0
    fprintf(fid, '%s\n', L{:});
    fclose(fid);
end
fprintf('\n%s\n', L{:});
lignes_resume = L; %#ok<NASGU>
save(fullfile(dossier, 'coeurs_agregation.mat'), 'C', 'P', 'opts', 'lignes_resume');
try
    figure_coeurs(dossier);
catch err
    warning('figure_coeurs non exécutée : %s', err.message);
end
end
