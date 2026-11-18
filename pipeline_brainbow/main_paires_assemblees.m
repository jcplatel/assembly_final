function [T, A] = main_paires_assemblees(config)
%% =========================================================================
%  MAIN_PAIRES_ASSEMBLEES : « mesurer au lieu de compter » (10 nov., accepté
%  par l'utilisateur : « ok rajoute cette analyse, on va voir ce que ça donne »)
%  -------------------------------------------------------------------------
%  L'étape 3 fait un test par assemblée (p_minp) puis l'étape 4 compte les
%  assemblées significatives. Une assemblée seule a peu de cellules colorées :
%  même avec 3 sœurs elle est rarement significative, et son information est
%  perdue. Ici, UNE mesure par session, toutes les assemblées ensemble :
%    1. dans chaque assemblée, nombre de paires de cellules de même couleur ;
%    2. une paire de la couleur c compte 1 / (part de c parmi les cellules
%       colorées du champ) : deux cellules d'une couleur rare ensemble sont
%       plus surprenantes que deux d'une couleur fréquente ;
%    3. somme sur les assemblées de la session (Q) ;
%    4. même somme après mélange des couleurs entre les cellules colorées
%       du champ (attendu par hasard) ; excès = Q / attendu - 1.
%  Unité : l'animal (moyenne de ses sessions), Wilcoxon unilatéral.
%  Mêmes choix que l'étape 3 : urne = cellules colorées du champ (couleurs
%  1 à 6, doublons et ROI hors taille retirés), couleurs de moins de
%  config.clonal.nmin_couleur cellules ignorées (poids 0), assemblées au K
%  retenu par l'étape 2 (consensus compris).
%  Deux mélanges : 'globale' (test) et 'blocs' (contrôle spatial : couleurs
%  mélangées entre cellules d'un même carré, comme la variante champ_blocs).
%  Deux jeux d'assemblées : 'valides' = critères de l'étape 4 (au moins
%  config.agreg.min_cellules_asm cellules et min_SCE_asm activations) : le
%  test principal ; 'toutes' : sensibilité.
%  Aussi, sans poids : paires de même couleur en plus par assemblée.
%
%  Limites : ne dit pas QUELLE assemblée est clonale ; trompé comme le test
%  de l'étape 3 si la couleur est liée à la qualité du signal.
%  Simulation Python (8 animaux x 4 sessions x 6 assemblées de 20 cellules,
%  3 sœurs par assemblée) : détecté 96 fois sur 100, contre 71 en comptant
%  les assemblées à p_minp < 0,05 ; sans sœurs : 0 faux positif sur 24.
%
%  Lancer :   main_paires_assemblees
%             config = config_brainbow(); config.seulement = {'444175'}; main_paires_assemblees(config)
%  Il faut les étapes 1 et 2 (pas l'étape 3). Ne modifie rien d'autre.
%  Écrit dans <PathSave>\paires_assemblees\ : paires_assemblees_sessions.csv,
%  paires_assemblees_animaux.csv, paires_assemblees.mat,
%  resume_paires_assemblees.txt, paires_assemblees.png.
%  =========================================================================

if nargin < 1, config = config_brainbow(); end
addpath(genpath(fileparts(mfilename('fullpath'))));
etat_rng = rng;
remettre_rng = onCleanup(@() rng(etat_rng)); %#ok<NASGU>
config = dossiers_analyse(config);
sessions = lire_table_sessions(config);
a_traiter = 1:height(sessions);
if isfield(config, 'seulement') && ~isempty(config.seulement)
    seulement = cellstr(config.seulement);
    a_traiter = find(ismember(sessions.nom, seulement) | ismember(sessions.id_animal, seulement))';
end
nPerm = 2000;
if isfield(config, 'paires_asm') && isfield(config.paires_asm, 'nPerm'), nPerm = config.paires_asm.nPerm; end
nmin = config.clonal.nmin_couleur;
min_cel = config.agreg.min_cellules_asm;
min_sce = config.agreg.min_SCE_asm;
dossier = fullfile(config.PathSave, 'paires_assemblees');
if ~isfolder(dossier), mkdir(dossier); end
fprintf('=== PAIRES DE MÊME COULEUR DANS LES ASSEMBLÉES : %d session(s), %d mélanges ===\n', numel(a_traiter), nPerm);

jeux = {'valides', 'toutes'};
modes = {'globale', 'blocs'};
lignes = cell(0, 16);
for k = a_traiter
    session = table2struct(sessions(k, :));
    try
        [asm, PCl, K_s] = charger_assemblees(session);
        if isempty(asm), continue; end
        path = avec_separateur(session.dossier_suite2p);
        [outline_gcampx, ~, gcamp_centroid] = charger_geometrie(path);
        colorcell = charger_couleurs(path, session.fichier_colorcell);
        colorcell = couleurs_sans_doublons(colorcell, session);
        colorcell = colorcell(:);
        NCell = numel(colorcell);
        verifier_alignement(NCell, colorcell, outline_gcampx, session.nom);

        % urne : cellules colorées du champ
        idx = find(colorcell >= 1 & colorcell <= 6);
        N = numel(idx);
        cc = colorcell(idx);
        Kc = accumarray(cc, 1, [6 1])';
        poids = zeros(1, 6);
        ok_c = Kc >= nmin;
        poids(ok_c) = N ./ Kc(ok_c);
        pos = zeros(NCell, 1);
        pos(idx) = 1:N;

        % assemblées : membres colorés, critères de l'étape 4
        nA = numel(asm);
        M = zeros(nA, N);
        n_tot = zeros(nA, 1);
        for a = 1:nA
            m = unique(asm{a}(:));
            m = m(m >= 1 & m <= NCell);
            n_tot(a) = numel(m);
            p = pos(m);
            M(a, p(p > 0)) = 1;
        end
        n_act = nan(nA, 1);
        if ~isempty(PCl)
            if size(PCl, 1) ~= nA && size(PCl, 2) == nA, PCl = PCl'; end
            if size(PCl, 1) == nA, n_act = sum(PCl > 0, 2); end
        end
        valide = n_tot >= min_cel;
        if min_sce > 0, valide = valide & (n_act >= min_sce); end
        choix = {valide, true(nA, 1)};

        % blocs (comme la variante champ_blocs de l'étape 3)
        taille_bloc_px = config.clonal.taille_bloc_px;
        if ~isnan(config.clonal.taille_bloc_um) && ~isnan(session.um_par_px)
            taille_bloc_px = config.clonal.taille_bloc_um / session.um_par_px;
        end
        [~, ~, num_bloc] = unique(floor(gcamp_centroid(idx, :) / taille_bloc_px), 'rows');
        groupes = {{(1:N)'}, accumarray(num_bloc, (1:N)', [], @(v) {v})};

        X = zeros(N, 6);
        X(sub2ind([N 6], (1:N)', cc)) = 1;
        [Q_obs, S_obs] = sommes(M, X, poids, choix);
        for g = 1:2
            fixer_flux_aleatoire(session.graine, 30 + g, 1);
            Qp = zeros(nPerm, 2);
            Sp = zeros(nPerm, 2);
            for p = 1:nPerm
                c2 = cc;
                gr = groupes{g};
                for b = 1:numel(gr)
                    v = gr{b};
                    c2(v) = cc(v(randperm(numel(v))));
                end
                Xp = zeros(N, 6);
                Xp(sub2ind([N 6], (1:N)', c2)) = 1;
                [Qp(p, :), Sp(p, :)] = sommes(M, Xp, poids, choix);
            end
            for j = 1:2
                nAj = sum(choix{j});
                if nAj == 0, continue; end
                att = mean(Qp(:, j));
                lignes(end+1, :) = {session.id_animal, session.nom, K_s, jeux{j}, modes{g}, nAj, ...
                    sum(sum(M(choix{j}, :), 2)), N, Q_obs(j), att, std(Qp(:, j)), ...
                    (1 + sum(Qp(:, j) >= Q_obs(j) * (1 - 1e-12))) / (nPerm + 1), ...
                    S_obs(j), mean(Sp(:, j)), (S_obs(j) - mean(Sp(:, j))) / nAj, Q_obs(j) / att - 1}; %#ok<AGROW>
            end
        end
        fprintf('[%s] K = %d : %d assemblées (%d valides)\n', session.nom, K_s, nA, sum(valide));
    catch err
        fprintf(2, '[%s] : ÉCHEC : %s\n', session.nom, err.message);
    end
end
T = cell2table(lignes, 'VariableNames', {'animal', 'session', 'K', 'assemblees', 'melange', 'n_assemblees', ...
    'n_membres_colores', 'n_colorees_champ', 'Q_obs', 'Q_attendu', 'Q_sd', 'p_session', 'paires_obs', ...
    'paires_attendues', 'paires_en_plus_par_assemblee', 'exces'});
A = table();
if isempty(T)
    fprintf('Aucune assemblée trouvée (étapes 1 et 2 faites ?).\n');
    return
end
T.z = (T.Q_obs - T.Q_attendu) ./ T.Q_sd;
writetable(T, fullfile(dossier, 'paires_assemblees_sessions.csv'));

% Animaux : moyenne des sessions ; Wilcoxon unilatéral
B = table();
for j = 1:2
    for g = 1:2
        t = T(strcmp(T.assemblees, jeux{j}) & strcmp(T.melange, modes{g}) & isfinite(T.exces), :);
        if isempty(t), continue; end
        [G, anim] = findgroups(t.animal);
        ex = splitapply(@mean, t.exces, G);
        pl = splitapply(@mean, t.paires_en_plus_par_assemblee, G);
        ns = splitapply(@numel, t.exces, G);
        A = [A; table(repmat(jeux(j), numel(anim), 1), repmat(modes(g), numel(anim), 1), anim, ns, ex, pl, ...
            'VariableNames', {'assemblees', 'melange', 'animal', 'n_sessions', 'exces', ...
            'paires_en_plus_par_assemblee'})]; %#ok<AGROW>
        p_w = NaN;
        if numel(ex) >= 2 && any(ex ~= 0), p_w = signrank(ex, 0, 'tail', 'right'); end
        z_c = sum(t.Q_obs - t.Q_attendu) / sqrt(sum(t.Q_sd .^ 2));
        B = [B; table(jeux(j), modes(g), height(t), sum(t.n_assemblees), numel(anim), sum(ex > 0), ...
            mean(ex), median(ex), p_w, sum(t.p_session < 0.05), z_c, 0.5 * erfc(z_c / sqrt(2)), mean(pl), ...
            'VariableNames', {'assemblees', 'melange', 'n_sessions', 'n_assemblees', 'n_animaux', ...
            'n_animaux_positifs', 'exces_moyen', 'exces_median', 'p_wilcoxon', 'n_sessions_p05', ...
            'z_sessions', 'p_sessions', 'paires_en_plus_par_assemblee'})]; %#ok<AGROW>
    end
end
writetable(A, fullfile(dossier, 'paires_assemblees_animaux.csv'));

% Résumé
L = {};
L{end+1} = ['PAIRES DE MÊME COULEUR DANS LES ASSEMBLÉES (main_paires_assemblees ; une mesure par session, ' ...
    'toutes les assemblées ensemble ; paires pondérées par la rareté de la couleur) :'];
for r = 1:height(B)
    if strcmp(B.assemblees{r}, 'valides') && strcmp(B.melange{r}, 'globale')
        role = 'TEST PRINCIPAL';
    elseif strcmp(B.melange{r}, 'blocs')
        role = 'contrôle spatial';
    else
        role = 'sensibilité';
    end
    L{end+1} = sprintf(['   Assemblées %s, mélange %s (%s) : %d sessions, %d assemblées, %d animaux ; excès ' ...
        '%+.1f %% (moyenne des animaux ; médiane %+.1f %%), %d/%d animaux > 0, Wilcoxon p = %s ; sessions ' ...
        'combinées z = %.2f (p = %s), %d/%d sessions à p < 0,05 ; sans poids : %+.2f paire de même couleur ' ...
        'en plus par assemblée.'], B.assemblees{r}, B.melange{r}, role, B.n_sessions(r), B.n_assemblees(r), ...
        B.n_animaux(r), 100 * B.exces_moyen(r), 100 * B.exces_median(r), B.n_animaux_positifs(r), ...
        B.n_animaux(r), texte_p(B.p_wilcoxon(r)), B.z_sessions(r), texte_p(B.p_sessions(r)), ...
        B.n_sessions_p05(r), B.n_sessions(r), B.paires_en_plus_par_assemblee(r)); %#ok<AGROW>
end
L{end+1} = sprintf(['   Assemblées valides = critères de l''étape 4 (au moins %d cellules et %d activations). ' ...
    'Ne dit pas quelle assemblée est clonale ; « sessions combinées » suppose les sessions indépendantes ' ...
    '(optimiste) : la généralisation repose sur le Wilcoxon des animaux.'], min_cel, min_sce);
fid = fopen(fullfile(dossier, 'resume_paires_assemblees.txt'), 'w', 'n', 'UTF-8');
if fid >= 0
    fprintf(fid, '%s\n', L{:});
    fclose(fid);
end
fprintf('\n%s\n', L{:});
save(fullfile(dossier, 'paires_assemblees.mat'), 'T', 'A', 'B', 'nPerm', 'min_cel', 'min_sce');
try
    figure_paires_assemblees(dossier);
catch err
    warning('figure_paires_assemblees non exécutée : %s', err.message);
end
end


function [Q, S] = sommes(M, X, poids, choix)
% paires de même couleur dans les assemblées : Q pondéré, S sans poids,
% pour chaque jeu d'assemblées (choix{j} = assemblées gardées)
cnt = M * X;
pm = cnt .* (cnt - 1) / 2;
q = pm * poids';
s = sum(pm(:, poids > 0), 2);
Q = zeros(1, numel(choix));
S = zeros(1, numel(choix));
for j = 1:numel(choix)
    Q(j) = sum(q(choix{j}));
    S(j) = sum(s(choix{j}));
end
end


function [asm, PCl, K_s] = charger_assemblees(session)
% assemblées au K retenu par l'étape 2 (0 = consensus), avec leurs activations
asm = {};
PCl = [];
K_s = NaN;
f = fullfile(session.dossier, 'choix_K.mat');
if ~isfile(f), return; end
c = load(f);
if ~isfield(c, 'statut_choix') || ~strcmp(c.statut_choix, 'ok'), return; end
K_s = c.best_K;
fa = fichier_assemblees(session.dossier, session.dossier_source, K_s);   % v1 ou v2 selon l'analyse
if ~isfile(fa), return; end
S = load(fa, 'assemblyortho', 'PCl_final');
asm = S.assemblyortho;
if isfield(S, 'PCl_final'), PCl = S.PCl_final; end
end


function t = texte_p(p)
if isnan(p), t = 'NA'; elseif p < 0.001, t = '< 0.001'; else, t = sprintf('%.3f', p); end
end
