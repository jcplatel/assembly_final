function T = comparer_analyses(racines)
%% =========================================================================
%  COMPARER_ANALYSES : tableau de toutes les analyses déjà faites (tous les
%  dossiers d'étape 1 et toutes les analyses dedans), avec leurs réglages
%  (seuil des SCE, règle du K, K utilisés) et leur résultat clonal
%  (11 nov., demande de l'utilisateur : « trouver les paramètres qui me
%  donnent le plus d'assemblées biaisées »)
%  -------------------------------------------------------------------------
%  NE RECALCULE ET NE MODIFIE RIEN. Cherche, sous chaque dossier donné, tous
%  les fichiers agregation\<variante>\resultats_biais_clonal.mat écrits par
%  l'étape 4, et lit pour chacun :
%   - le résultat (bilan de l'étape 4) : assemblées testables, à p_minp <
%     0,05, attendues par hasard, excès en points, biaisées après FDR,
%     Wilcoxon et GLMM, score sur toutes les p, proportion estimée ;
%   - les réglages : règle et seuil des SCE (opts_utilisees de
%     pretraitement.mat, toutes les sessions), règle du K
%     (best_k_par_session.mat), K retenus (colonne K des résultats ;
%     0 = consensus), K calculés à l'étape 1 (dossiers k<K>).
%
%  ATTENTION, à lire avant de s'en servir : garder le réglage qui donne le
%  plus d'assemblées biaisées, c'est choisir d'après le résultat. Sans aucun
%  effet réel, le meilleur de N réglages paraît quand même meilleur que le
%  hasard, et ses p ne valent plus rien. Ce tableau sert à voir si l'effet
%  est STABLE d'un réglage à l'autre (c'est un argument) et quels réglages
%  le font varier ; le réglage final doit se justifier sans les couleurs
%  (qualité des SCE, des assemblées), et le nombre de réglages essayés doit
%  être dit.
%
%  Lancer :   comparer_analyses('E:\Data\Aurelie\analysis\Final\assembly\nocues')
%             comparer_analyses({'E:\...\nocues', 'D:\autre_dossier'})
%  Sans argument : le dossier qui contient config.PathSave.
%  Écrit comparaison_analyses.csv et resume_comparaison_analyses.txt dans le
%  premier dossier donné ; le tableau est aussi affiché, rangé par excès
%  (test principal champ_globale d'abord).
%  =========================================================================

addpath(genpath(fileparts(mfilename('fullpath'))));
if nargin < 1 || isempty(racines)
    config = config_brainbow();
    racines = {fileparts(regexprep(char(config.PathSave), '[\\/]+$', ''))};
end
racines = cellstr(racines);
fichiers = {};
for r = 1:numel(racines)
    L = dir(fullfile(racines{r}, '**', 'resultats_biais_clonal.mat'));
    for k = 1:numel(L)
        fichiers{end+1} = fullfile(L(k).folder, L(k).name); %#ok<AGROW>
    end
end
fichiers = unique(fichiers, 'stable');
fprintf('=== COMPARAISON DES ANALYSES : %d résultat(s) trouvé(s) ===\n', numel(fichiers));
T = table();
if isempty(fichiers), return; end

cache_etape1 = containers.Map('KeyType', 'char', 'ValueType', 'any');
lignes = cell(numel(fichiers), 1);
for k = 1:numel(fichiers)
    try
        f = fichiers{k};
        d_var = fileparts(f);
        [d_agreg, variante] = fileparts(d_var);
        d_analyse = fileparts(d_agreg);
        [d_parent, nom_analyse] = fileparts(d_analyse);
        [d_etape1, nom_parent] = fileparts(d_parent);
        if strcmpi(nom_parent, 'analyses')
            analyse = nom_analyse;
        else
            d_etape1 = d_analyse;
            analyse = '';
        end
        [~, nom_etape1] = fileparts(d_etape1);
        S = load(f);
        b = struct();
        if isfield(S, 'bilan'), b = S.bilan; end
        info = dir(f);

        % réglages de l'étape 1 (une fois par dossier)
        if ~isKey(cache_etape1, d_etape1), cache_etape1(d_etape1) = reglages_etape1(d_etape1); end
        E = cache_etape1(d_etape1);

        % règle du K et K retenus
        regle_K = '';
        fk = fullfile(d_analyse, 'best_k_par_session.mat');
        if isfile(fk)
            Bk = load(fk);
            if isfield(Bk, 'regle_choix_K'), regle_K = texte_regle(Bk.regle_choix_K); end
        end
        K_txt = ''; n_sessions = NaN; n_asm = NaN; taille_med = NaN;
        if isfield(S, 'R') && istable(S.R) && ~isempty(S.R)
            R = S.R;
            n_asm = height(R);
            if ismember('n_total', R.Properties.VariableNames), taille_med = median(R.n_total); end
            if ismember('session_u', R.Properties.VariableNames)
                [~, premier] = unique(R.session_u);
                n_sessions = numel(premier);
                if ismember('K', R.Properties.VariableNames)
                    Ks = R.K(premier);
                    if all(Ks == 0)
                        K_txt = 'consensus';
                    else
                        K_txt = sprintf('%g [%g-%g]', median(Ks), min(Ks), max(Ks));
                    end
                end
            end
        end
        nt = champ(b, 'n_testables'); ns = champ(b, 'n_sig'); na = champ(b, 'n_attendu');
        nfdr = champ(b, 'n_biaisees');
        % anciens résultats sans bilan : recalcul depuis le tableau des assemblées
        % (p_minp, testable, taux_nul ; Benjamini-Hochberg à 10 % sur les testables)
        if isnan(nt) && isfield(S, 'R') && istable(S.R) && all(ismember({'p_minp', 'testable'}, S.R.Properties.VariableNames))
            R = S.R;
            t = logical(R.testable) & isfinite(R.p_minp);
            if ismember('exclue_criteres', R.Properties.VariableNames), t = t & ~logical(R.exclue_criteres); end
            nt = sum(t);
            ns = sum(R.p_minp(t) < 0.05);
            if ismember('taux_nul', R.Properties.VariableNames), na = sum(R.taux_nul(t), 'omitnan'); end
            pv = sort(R.p_minp(t));
            kmax = find(pv(:)' <= 0.10 * (1:nt) / max(nt, 1), 1, 'last');
            nfdr = 0;
            if ~isempty(kmax), nfdr = kmax; end
            if ismember('animal', R.Properties.VariableNames), b.n_animaux = numel(unique(R.animal(t))); end
            b.couleurs_biaisees = '(recalculé)';
        end
        pi1_ic = champ2(b, 'pi1_ic'); OR_ic = champ2(b, 'OR_ic');
        lignes{k} = table({nom_etape1}, {analyse}, {variante}, {E.sce}, {E.K_calcules}, {regle_K}, {K_txt}, ...
            n_sessions, champ(b, 'n_animaux'), n_asm, taille_med, nt, ns, 100 * ns / max(nt, 1), ...
            100 * na / max(nt, 1), 100 * (ns - na) / max(nt, 1), nfdr, ...
            {champ_txt(b, 'couleurs_biaisees')}, champ(b, 'p_wilcoxon'), champ(b, 'OR'), OR_ic(1), ...
            champ(b, 'p_glmm'), champ(b, 'score_moyen'), champ(b, 'p_wil_score'), 100 * champ(b, 'pi1'), ...
            100 * pi1_ic(1), 100 * pi1_ic(2), {datestr(info.datenum, 'yyyy-mm-dd HH:MM')}, {d_analyse}, ...
            'VariableNames', {'etape1', 'analyse', 'variante', 'seuil_SCE', 'K_calcules', 'regle_K', ...
            'K_retenus', 'n_sessions', 'n_animaux', 'n_assemblees', 'taille_mediane', 'n_testables', ...
            'n_p05', 'pct_p05', 'pct_attendu', 'exces_points', 'n_biaisees_FDR', 'couleurs_FDR', ...
            'p_wilcoxon', 'OR', 'OR_ic_bas', 'p_glmm', 'score', 'p_score', 'pct_estime', 'pct_estime_bas', ...
            'pct_estime_haut', 'date_etape4', 'dossier'});
    catch err
        fprintf(2, '%s : ÉCHEC : %s\n', fichiers{k}, err.message);
    end
end
T = vertcat(lignes{~cellfun(@isempty, lignes)});
if isempty(T), return; end
principal = strcmp(T.variante, 'champ_globale');
[~, ordre] = sortrows(table(~principal, -T.exces_points));
T = T(ordre, :);
writetable(T, fullfile(racines{1}, 'comparaison_analyses.csv'));

% Résumé
P = T(strcmp(T.variante, 'champ_globale'), :);
L = {};
L{end+1} = sprintf('COMPARAISON DES ANALYSES (comparer_analyses) : %d résultat(s), dont %d au test principal (champ_globale).', ...
    height(T), height(P));
L{end+1} = ['   ATTENTION : choisir le réglage qui donne le plus d''assemblées biaisées, c''est choisir d''après ' ...
    'le résultat ; ses p ne valent plus. À lire comme un tableau de stabilité.'];
if ~isempty(P)
    L{end+1} = sprintf(['   Test principal, d''un réglage à l''autre : excès médian %+.1f points [%+.1f à %+.1f] ; ' ...
        '%d réglage(s) sur %d avec un excès > 0 ; assemblées biaisées après FDR : médiane %g [%g-%g].'], ...
        median(P.exces_points, 'omitnan'), min(P.exces_points), max(P.exces_points), sum(P.exces_points > 0), ...
        height(P), median(P.n_biaisees_FDR, 'omitnan'), min(P.n_biaisees_FDR), max(P.n_biaisees_FDR));
    L{end+1} = '   Test principal, rangé par excès (points au-dessus du hasard) :';
    for r = 1:height(P)
        L{end+1} = sprintf(['      %-22s | %-32s | SCE %-26s | K %-40s %-12s | %3d testables, %4.1f %% à p < 0,05 ' ...
            '(attendu %3.1f) = %+5.1f pts | FDR %d | Wilcoxon %s, GLMM %s, score %s | estimé %s %%'], ...
            P.etape1{r}, texte_vide(P.analyse{r}), P.seuil_SCE{r}, P.regle_K{r}, P.K_retenus{r}, P.n_testables(r), ...
            P.pct_p05(r), P.pct_attendu(r), P.exces_points(r), P.n_biaisees_FDR(r), texte_p(P.p_wilcoxon(r)), ...
            texte_p(P.p_glmm(r)), texte_p(P.p_score(r)), nb(P.pct_estime(r))); %#ok<AGROW>
    end
end
fid = fopen(fullfile(racines{1}, 'resume_comparaison_analyses.txt'), 'w', 'n', 'UTF-8');
if fid >= 0
    fprintf(fid, '%s\n', L{:});
    fclose(fid);
end
fprintf('\n%s\n', L{:});
fprintf('\nTableau complet (toutes les variantes) : %s\n', fullfile(racines{1}, 'comparaison_analyses.csv'));
end


function E = reglages_etape1(d_etape1)
% seuil des SCE (toutes les sessions du dossier) et K calculés
E = struct('sce', '?', 'K_calcules', '?');
L = dir(fullfile(d_etape1, '*', 'pretraitement.mat'));
textes = {};
Kmin = Inf; Kmax = -Inf;
for k = 1:numel(L)
    try
        P = load(fullfile(L(k).folder, L(k).name), 'opts_utilisees');
        o = P.opts_utilisees;
        regle = 'fixe';
        if isfield(o, 'sce_regle'), regle = char(o.sce_regle); end
        switch regle
            case 'percentile'
                t = sprintf('percentile %g', val(o, 'sce_percentile'));
                if isfield(o, 'sce_minimum'), t = sprintf('%s (min %g)', t, o.sce_minimum); end
            case 'fdr'
                t = sprintf('fdr %g %%, décalage %g s', 100 * val(o, 'sce_fdr'), val(o, 'sce_decalage_fdr_s'));
                if isfield(o, 'sce_minimum'), t = sprintf('%s (min %g)', t, o.sce_minimum); end
            otherwise
                t = sprintf('fixe %g', val(o, 'sce_n_cells_threshold'));
        end
        textes{end+1} = t; %#ok<AGROW>
    catch
    end
    D = dir(fullfile(L(k).folder, 'k*'));
    for j = 1:numel(D)
        v = sscanf(D(j).name, 'k%d');
        if D(j).isdir && ~isempty(v), Kmin = min(Kmin, v); Kmax = max(Kmax, v); end
    end
end
if ~isempty(textes)
    [u, ~, ic] = unique(textes);
    if numel(u) == 1
        E.sce = u{1};
    else
        nb_u = accumarray(ic, 1);
        [~, i] = max(nb_u);
        E.sce = sprintf('MIXTE (%d règles ; %s : %d/%d)', numel(u), u{i}, nb_u(i), numel(textes));
    end
end
if isfinite(Kmin), E.K_calcules = sprintf('%d-%d', Kmin, Kmax); end
end


function t = texte_regle(p)
% règle du K (config.choixK tel que sauvé par l'étape 2) : règle, K fixé,
% réglages du consensus, K proposés
t = '';
if ischar(p) || isstring(p), t = char(p); return; end
if ~isstruct(p), return; end
if isfield(p, 'regle'), t = char(p.regle); end
if isfield(p, 'K_fixe') && ~isempty(p.K_fixe), t = sprintf('K fixe %g', p.K_fixe); end
if strcmp(t, 'consensus') && isfield(p, 'consensus') && isstruct(p.consensus)
    c = p.consensus;
    if isfield(c, 'seuil') && isfield(c, 'min_K')
        t = sprintf('consensus (seuil %g, min_K %g)', c.seuil, c.min_K);
    end
end
if isfield(p, 'liste_K') && ~isempty(p.liste_K)
    t = sprintf('%s, K %g-%g', t, min(p.liste_K), max(p.liste_K));
end
end


function v = val(o, nom)
v = NaN;
if isfield(o, nom) && isnumeric(o.(nom)) && isscalar(o.(nom)), v = double(o.(nom)); end
end


function v = champ(b, nom)
v = NaN;
if isfield(b, nom) && isnumeric(b.(nom)) && ~isempty(b.(nom)), v = double(b.(nom)(1)); end
end


function v = champ2(b, nom)
v = [NaN NaN];
if isfield(b, nom) && isnumeric(b.(nom)) && numel(b.(nom)) >= 2, v = double(b.(nom)(1:2)); end
end


function t = champ_txt(b, nom)
t = '';
if isfield(b, nom) && (ischar(b.(nom)) || isstring(b.(nom))), t = char(b.(nom)); end
end


function t = texte_vide(t)
if isempty(t), t = '(principale)'; end
end


function t = texte_p(p)
if isnan(p), t = 'NA'; elseif p < 0.001, t = '<0.001'; else, t = sprintf('%.3f', p); end
end


function t = nb(x)
if isnan(x), t = 'NA'; else, t = sprintf('%.0f', x); end
end
