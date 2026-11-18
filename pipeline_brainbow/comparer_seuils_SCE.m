function [T_sessions, T_asm] = comparer_seuils_SCE(lancer, sessions_test, percentiles, liste_K, K_commun)
%COMPARER_SEUILS_SCE  Étape 1 (+ choix du K) avec le seuil des SCE au 95e et
%   au 99e percentile, sur quelques sessions, dans deux dossiers séparés ;
%   puis comparaison des assemblées, SANS les couleurs.
%
%   comparer_seuils_SCE              % lance les calculs puis compare (long)
%   comparer_seuils_SCE(false)       % compare seulement (calculs déjà faits)
%   comparer_seuils_SCE(true, {}, [95 99], 4:16)     % moins de K : plus rapide
%   comparer_seuils_SCE(true, {'444119_220919_plane0_sans_indice'})
%
%   Dossiers de l'étape 1 : <racine>\test_SCE_p95 et <racine>\test_SCE_p99
%   (racine = dossier qui contient resultats). Réglage 'valide' sauf le
%   seuil des SCE : règle 'percentile', minimum absolu 5 cellules. Étapes 1
%   et 2 seulement (aucune couleur) ; figures PNG coupées (plus rapide).
%   Les résultats se regardent aussi dans l'interface : « Step 1 folder » =
%   un de ces deux dossiers. Si le calcul s'arrête, relancer : les sessions
%   déjà faites sont gardées en cache.
%
%   Sessions par défaut : choisies sur la taille du champ seulement (89 à
%   887 cellules, réparties sur les 41 sessions de moins de 2000 cellules),
%   dont 444119_220919_plane0 (session de l'ancienne figure).
%
%   Comparaison, pour chaque session et chaque percentile :
%     seuil, n_SCE, pct_hasard  seuil des SCE, nombre de SCE, part attendue par hasard
%     K_choisi, n_asm, tailles  K retenu à l'étape 2 et ses assemblées
%     Kc_n_asm, Kc_tailles      assemblées au même K pour tous (K_commun, défaut 9) :
%                               effet du seuil seul, sans l'effet du choix du K
%   Pour chaque assemblée (au K retenu et au K commun) :
%     taille
%     activite_rel   transients de ses cellules / médiane du champ ; >> 1 :
%                    assemblée de cellules très actives (à surveiller)
%     meilleur_J     meilleur indice de Jaccard avec une assemblée de l'autre
%                    percentile (même niveau : K retenu ou K commun)
%     retrouvee      meilleur_J >= 0,3
%   Écrit <racine>\comparaison_seuils_SCE.xlsx (feuilles 'sessions' et 'assemblees').

if nargin < 1 || isempty(lancer), lancer = true; end
if nargin < 2 || isempty(sessions_test)
    sessions_test = {'443695_220728_plane0_sans_indice', '444113_220930_plane1_sans_indice', ...
        '444175_221121_plane0_sans_indice', '444175_221125_plane1_sans_indice', ...
        '444119_220919_plane0_sans_indice', '444152_221129_plane1_sans_indice'};
end
sessions_test = cellstr(sessions_test);
if nargin < 3 || isempty(percentiles), percentiles = [95 99]; end
config0 = config_brainbow();
if nargin < 4 || isempty(liste_K), liste_K = config0.liste_K; end
if nargin < 5 || isempty(K_commun), K_commun = 9; end
seuil_J = 0.3;
dossier_pipeline = fileparts(mfilename('fullpath'));
addpath(genpath(dossier_pipeline));

racine = fileparts(regexprep(config0.PathSave, '[\\/]+$', ''));
nP = numel(percentiles);
dossiers = arrayfun(@(p) fullfile(racine, sprintf('test_SCE_p%g', p)), percentiles, 'UniformOutput', false);
noms_p = arrayfun(@(p) sprintf('p%g', p), percentiles, 'UniformOutput', false);

%% 1. Étapes 1 et 2 pour chaque percentile
if lancer
    t_total = tic;
    for ip = 1:nP
        cfg = config_brainbow();
        cfg.PathSave = dossiers{ip};
        cfg.analyse = '';
        cfg.pre.sce_regle = 'percentile';
        cfg.pre.sce_percentile = percentiles(ip);
        cfg.pre.sce_minimum = 5;                    % minimum des règles calculées (sce_n_cells_threshold : règle 'fixe' seulement)
        cfg.liste_K = liste_K;
        cfg.choixK.k_start = min(liste_K);
        cfg.seulement = sessions_test;
        cfg.faire_etape1 = true;
        cfg.faire_etape2 = true;
        cfg.faire_etape3 = false;
        cfg.faire_etape3b = false;
        cfg.faire_etape4 = false;
        cfg.faire_figures = false;
        fprintf('\n===== Seuil des SCE au percentile %g : %s =====\n', percentiles(ip), dossiers{ip});
        main_brainbow(cfg);
    end
    fprintf('\nCalculs terminés (%.1f h).\n', toc(t_total) / 3600);
end

%% 2. Lecture des résultats (aucune couleur)
nS = numel(sessions_test);
R = cell(nS, nP);
for s = 1:nS
    for ip = 1:nP
        R{s, ip} = lire_session(fullfile(dossiers{ip}, sessions_test{s}), K_commun);
    end
end

%% 3. Comparaison
niveaux = {'K retenu', 'K commun'};
champs_niveau = {'asm_choisi', 'asm_commun'};
lignes_s = cell(0, 13);
lignes_a = cell(0, 11);
fprintf('\n===== Comparaison des seuils des SCE (sans couleur) =====\n');
fprintf('Assemblée « retrouvée » : Jaccard >= %.1f avec une assemblée de l''autre percentile.\n', seuil_J);
for s = 1:nS
    nom = sessions_test{s};
    ncell = NaN;
    for ip = 1:nP
        if R{s, ip}.ok, ncell = R{s, ip}.n_cellules; break; end
    end
    fprintf('\n%s (%d cellules)\n', nom, ncell);
    for ip = 1:nP
        r = R{s, ip};
        if ~r.ok
            fprintf('  %-4s : %s\n', noms_p{ip}, r.message);
            continue
        end
        autres = setdiff(1:nP, ip);
        non_retr = zeros(1, 2);
        for iv = 1:2
            A = r.(champs_niveau{iv});
            if iv == 1, K = r.K_choisi; else, K = K_commun; end
            for c = 1:numel(A)
                cel = A{c}(:);
                act = median(r.activite(cel)) / r.activite_mediane;
                J = 0;
                for io = autres
                    if R{s, io}.ok
                        B = R{s, io}.(champs_niveau{iv});
                        for b = 1:numel(B)
                            J = max(J, jaccard(cel, B{b}(:)));
                        end
                    end
                end
                retrouvee = J >= seuil_J;
                non_retr(iv) = non_retr(iv) + ~retrouvee;
                lignes_a(end+1, :) = {nom, ncell, percentiles(ip), niveaux{iv}, K, c, numel(cel), ...
                    round(act, 2), round(J, 2), retrouvee, sprintf('%d ', sort(cel(:))')}; %#ok<AGROW>
            end
        end
        lignes_s(end+1, :) = {nom, ncell, percentiles(ip), r.seuil, r.n_SCE, round(r.pct_hasard, 1), ...
            r.K_choisi, numel(r.asm_choisi), tailles(r.asm_choisi), non_retr(1), ...
            numel(r.asm_commun), tailles(r.asm_commun), non_retr(2)}; %#ok<AGROW>
        fprintf(['  %-4s : seuil %d cellules, %d SCE (%.0f %% hasard) | K retenu %d : %d assemblée(s) [%s] ' ...
            '| K = %d : %d [%s]\n'], noms_p{ip}, r.seuil, r.n_SCE, r.pct_hasard, r.K_choisi, ...
            numel(r.asm_choisi), tailles(r.asm_choisi), K_commun, numel(r.asm_commun), tailles(r.asm_commun));
    end
    % assemblées d'un percentile absentes de l'autre
    for ip = 1:nP
        for iv = 1:2
            m = strcmp(lignes_a(:, 1), nom) & cell2mat(lignes_a(:, 3)) == percentiles(ip) & ...
                strcmp(lignes_a(:, 4), niveaux{iv}) & ~cell2mat(lignes_a(:, 10));
            if any(m)
                txt = cellfun(@(t, a, j) sprintf('%d cellules (activité x%.1f, J %.2f)', t, a, j), ...
                    lignes_a(m, 7), lignes_a(m, 8), lignes_a(m, 9), 'UniformOutput', false);
                fprintf('     seulement au %s (%s) : %s\n', noms_p{ip}, niveaux{iv}, strjoin(txt, ' ; '));
            end
        end
    end
end

T_sessions = cell2table(lignes_s, 'VariableNames', {'session', 'n_cellules', 'percentile', 'seuil', ...
    'n_SCE', 'pct_hasard', 'K_choisi', 'n_asm', 'tailles', 'non_retrouvees', ...
    'Kc_n_asm', 'Kc_tailles', 'Kc_non_retrouvees'});
T_asm = cell2table(lignes_a, 'VariableNames', {'session', 'n_cellules', 'percentile', 'niveau', 'K', ...
    'assemblee', 'taille', 'activite_rel', 'meilleur_J', 'retrouvee', 'cellules'});

%% 4. Bilan sur toutes les sessions
fprintf('\n===== Bilan (%d sessions) =====\n', nS);
for iv = 1:2
    fprintf('%s :\n', niveaux{iv});
    for ip = 1:nP
        if isempty(T_asm), break; end
        m = T_asm.percentile == percentiles(ip) & strcmp(T_asm.niveau, niveaux{iv});
        mr = m & ~T_asm.retrouvee;
        fprintf(['  %-4s : %3d assemblées, taille médiane %4.1f, activité médiane x%.2f | ' ...
            'absentes de l''autre percentile : %d (taille médiane %s)\n'], noms_p{ip}, sum(m), ...
            median(T_asm.taille(m)), median(T_asm.activite_rel(m)), sum(mr), texte_median(T_asm.taille(mr)));
    end
end
fprintf(['Lecture : des assemblées présentes au %s seulement, petites et d''activité normale, ' ...
    'plaident pour le %s ;\nsi elles sont faites de cellules très actives (activité >> 1) ou ' ...
    'instables, le %s suffit.\n'], noms_p{1}, noms_p{1}, noms_p{end});

fichier = fullfile(racine, 'comparaison_seuils_SCE.xlsx');
try
    if isfile(fichier), delete(fichier); end
    writetable(T_sessions, fichier, 'Sheet', 'sessions');
    writetable(T_asm, fichier, 'Sheet', 'assemblees');
    fprintf('\nTableaux : %s\n', fichier);
catch err
    warning('comparaison_seuils_SCE.xlsx non écrit (ouvert dans Excel ?) : %s', err.message);
end
end


function r = lire_session(d, K_commun)
% Résultats des étapes 1 et 2 d'une session (aucune couleur)
r = struct('ok', false, 'message', '', 'n_cellules', NaN, 'seuil', NaN, 'n_SCE', NaN, ...
    'pct_hasard', NaN, 'K_choisi', NaN, 'activite', [], 'activite_mediane', NaN);
r.asm_choisi = {};
r.asm_commun = {};
fpre = fullfile(d, 'pretraitement.mat');
if ~isfile(fpre)
    r.message = 'pas encore calculée (étape 1)';
    return
end
P = load(fpre, 'Raster', 'TRace', 'opts', 'diagnostic');
r.ok = true;
r.n_cellules = size(P.Raster, 1);
r.n_SCE = numel(P.TRace);
if isfield(P.opts, 'sce_n_cells_threshold_calcule'), r.seuil = P.opts.sce_n_cells_threshold_calcule; end
if isfield(P, 'diagnostic') && isfield(P.diagnostic, 'n_SCE_hasard')
    r.pct_hasard = 100 * P.diagnostic.n_SCE_hasard / max(r.n_SCE, 1);
end
r.activite = full(sum(P.Raster > 0, 2));            % transients au repos de chaque cellule
actives = r.activite(r.activite > 0);
if isempty(actives), r.activite_mediane = 1; else, r.activite_mediane = median(actives); end
fch = fullfile(d, 'choix_K.mat');
if isfile(fch)
    C = load(fch, 'best_K', 'statut_choix');
    if strcmp(C.statut_choix, 'ok'), r.K_choisi = C.best_K; end
end
r.asm_choisi = lire_assemblees(d, r.K_choisi);
r.asm_commun = lire_assemblees(d, K_commun);
end


function A = lire_assemblees(d, K)
A = {};
if isnan(K), return; end
f = fullfile(d, sprintf('k%d', K), 'assemblees.mat');
if ~isfile(f), return; end
S = load(f, 'assemblyortho');
if isfield(S, 'assemblyortho'), A = S.assemblyortho; end
end


function J = jaccard(a, b)
J = numel(intersect(a, b)) / max(numel(union(a, b)), 1);
end


function t = tailles(A)
if isempty(A), t = '-'; return; end
t = strjoin(arrayfun(@(c) sprintf('%d', numel(A{c})), 1:numel(A), 'UniformOutput', false), ' ');
end


function t = texte_median(x)
if isempty(x), t = '-'; else, t = sprintf('%.1f', median(x)); end
end
