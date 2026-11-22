function T = registre_analyses(racine, mode)
%REGISTRE_ANALYSES  Tableau de toutes les analyses (22 nov. 2026, demandé par
%   l'utilisateur : « les noms des analyses deviennent incompréhensibles »).
%
%   T = registre_analyses()                     racine par défaut (dossier assembly)
%   T = registre_analyses(racine)
%   T = registre_analyses(racine, 'rescanner')  relit tout (ignore le cache)
%   C = registre_analyses(racine, 'codes')      seulement les codes (rapide, sans parcours)
%
%   Une analyse = un dossier d'étape 1 (<racine>\nocues\<dossier> ou
%   <racine>\cues\<dossier>, avec des <session>\pretraitement.mat) et ses
%   étapes 2 à 4 : principale (dans le dossier d'étape 1, « main ») ou
%   <dossier d'étape 1>\analyses\<nom>.
%
%   Chaque analyse reçoit un CODE stable (A001, A002...) : gardé dans
%   <racine>\registre_analyses.mat ; une nouvelle analyse prend le numéro
%   suivant, rien n'est renuméroté. Le tableau est aussi écrit dans
%   <racine>\registre_analyses.csv.
%
%   Colonnes : code, condition, date, dossier d'étape 1, analyse, version
%   des assemblées, réglages (résumé), étapes faites, sessions, fin du
%   dernier lancement, SCE (total), assemblées, testables, % p < 0,05
%   (p hypergéométrique non corrigée), % p_minp < 0,05, nombre q < 0,05 et
%   q < 0,10 (FDR de Benjamini-Hochberg de l'étape 4, test champ / globale).
%
%   Sources : parametres_analyse.txt (réglages), log_analyse_*.txt (dernier
%   lancement), best_k_par_session.mat, <session>\choix_K.mat,
%   <session>\clonal\stats_clonales.mat, agregation\champ_globale\
%   resultats_biais_clonal.csv, taille de Race dans pretraitement.mat
%   (nombre de SCE, lu sans charger le fichier).
%   Cache : une analyse n'est relue que si ses fichiers ont changé.

if nargin < 1 || isempty(racine)
    c = config_brainbow();
    racine = fileparts(fileparts(c.PathSave));     % ...\assembly
end
if nargin < 2, mode = ''; end
f_reg = fullfile(racine, 'registre_analyses.mat');
R = struct('codes', table(cell(0, 1), cell(0, 1), 'VariableNames', {'chemin', 'code'}), ...
    'cache', struct('cle', {}, 'signature', {}, 'ligne', {}), ...
    'cache_sce', struct('cle', {}, 'signature', {}, 'n', {}));
if isfile(f_reg)
    try
        S = load(f_reg, 'R');
        R = S.R;
    catch
    end
end
if strcmpi(mode, 'codes')
    T = R.codes;
    return
end
version_cache = 2;      % à augmenter quand le contenu d'une ligne change (cache des lignes relu)
if ~isfield(R, 'version_cache') || R.version_cache ~= version_cache
    R.cache = R.cache([]);
    R.version_cache = version_cache;
end
if strcmpi(mode, 'rescanner')
    R.cache = R.cache([]);
    R.cache_sce = R.cache_sce([]);
end

conditions = {'nocues', 'no cues'; 'cues', 'cues'};
lignes = {};
for ic = 1:size(conditions, 1)
    dc = fullfile(racine, conditions{ic, 1});
    L = dir(dc);
    L = L([L.isdir] & ~ismember({L.name}, {'.', '..'}));
    for i = 1:numel(L)
        d1 = fullfile(dc, L(i).name);
        pre = dir(fullfile(d1, '*', 'pretraitement.mat'));
        if isempty(pre), continue; end          % pas un dossier d'étape 1
        [n_sce, R] = sce_etape1(d1, pre, R);
        % analyse principale, puis analyses\<nom>
        noms = {''};
        A = dir(fullfile(d1, 'analyses'));
        A = A([A.isdir] & ~ismember({A.name}, {'.', '..'}));
        noms = [noms, sort({A.name})]; %#ok<AGROW>
        for j = 1:numel(noms)
            if isempty(noms{j}), da = d1; else, da = fullfile(d1, 'analyses', noms{j}); end
            [ligne, R] = ligne_analyse(da, d1, noms{j}, conditions{ic, 2}, numel(pre), n_sce, R);
            if ~isempty(ligne), lignes(end+1, :) = ligne; end %#ok<AGROW>
        end
    end
end

noms_col = {'Code', 'Condition', 'Date', 'Step1_folder', 'Analysis', 'Version', 'Settings', ...
    'Steps', 'Sessions', 'Last_run', 'SCE', 'Assemblies', 'Testable', 'Pct_p05', ...
    'Pct_pminp05', 'N_q05', 'N_q10', 'Step1_path', 'Analysis_path'};
if isempty(lignes)
    T = cell2table(cell(0, numel(noms_col)), 'VariableNames', noms_col);
else
    T = cell2table(lignes, 'VariableNames', noms_col);
    T = sortrows(T, 'Code', 'descend');                 % les plus récentes en haut
end
try
    save(f_reg, 'R');
catch err
    warning('registre_analyses : %s non écrit (%s)', f_reg, err.message);
end
try
    writetable(T, fullfile(racine, 'registre_analyses.csv'));
catch
end
end


%% ------------------------------------------------------------------------
function [ligne, R] = ligne_analyse(da, d1, nom, condition, n_etape1, n_sce, R)
% une ligne du tableau ; relue seulement si les fichiers ont changé
ligne = {};
sig = signature_analyse(da);
if isempty(nom) && sig(2) == 0 && n_etape1 == 0, return; end
cle = lower(da);
ic = find(strcmp({R.cache.cle}, cle), 1);
[code, R] = code_de(R, da);
if ~isempty(ic) && isequal(R.cache(ic).signature, sig)
    ligne = R.cache(ic).ligne;
else
    ligne = lire_analyse(da, d1, nom, condition, n_etape1);
    if isempty(ic), ic = numel(R.cache) + 1; end
    R.cache(ic).cle = cle;
    R.cache(ic).signature = sig;
    R.cache(ic).ligne = ligne;
end
ligne{1} = code;
ligne{11} = n_sce;          % SCE : cache à part (étape 1)
end


function sig = signature_analyse(da)
% dates et nombre des fichiers qui décident de la ligne
fichiers = {'parametres_analyse.txt', 'best_k_par_session.mat', ...
    fullfile('agregation', 'champ_globale', 'resultats_biais_clonal.csv'), ...
    fullfile('agregation', 'statut_sessions.csv')};
sig = zeros(1, numel(fichiers) + 3);
for k = 1:numel(fichiers)
    f = dir(fullfile(da, fichiers{k}));
    if ~isempty(f), sig(k) = f(1).datenum; end
end
Lg = dir(fullfile(da, 'log_analyse_*.txt'));
if ~isempty(Lg), sig(end - 2) = max([Lg.datenum]) + numel(Lg) + max([Lg.bytes]) / 1e9; end
sig(end - 1) = numel(dir(fullfile(da, '*', 'choix_K.mat')));
sig(end) = numel(dir(fullfile(da, '*', 'clonal', 'stats_clonales.mat')));
end


function [code, R] = code_de(R, da)
i = find(strcmpi(R.codes.chemin, da), 1);
if isempty(i)
    num = 0;
    if ~isempty(R.codes.code)
        num = max(cellfun(@(c) str2double(c(2:end)), R.codes.code));
    end
    code = sprintf('A%03d', num + 1);
    R.codes = [R.codes; {da, code}];
else
    code = R.codes.code{i};
end
end


function [n, R] = sce_etape1(d1, pre, R)
% nombre total de SCE de l'étape 1 (taille de Race, sans charger)
cle = lower(d1);
sig = [numel(pre), sum([pre.datenum])];
ic = find(strcmp({R.cache_sce.cle}, cle), 1);
if ~isempty(ic) && isequal(R.cache_sce(ic).signature, sig)
    n = R.cache_sce(ic).n;
    return
end
n = 0;
for k = 1:numel(pre)
    try
        w = whos('-file', fullfile(pre(k).folder, pre(k).name), 'Race');
        if ~isempty(w), n = n + w.size(2); end
    catch
    end
end
if isempty(ic), ic = numel(R.cache_sce) + 1; end
R.cache_sce(ic).cle = cle;
R.cache_sce(ic).signature = sig;
R.cache_sce(ic).n = n;
end


function ligne = lire_analyse(da, d1, nom, condition, n_etape1)
[~, n1] = fileparts(d1);
if isempty(nom), texte_nom = 'main'; else, texte_nom = nom; end
version = 'v1';
if ~isempty(regexpi(nom, '(^v2$|^v2_|_v2$)', 'once')), version = 'v2'; end

% réglages : parametres_analyse.txt (sinon : fichiers de l'étape 1 et 2)
P = lire_parametres(fullfile(da, 'parametres_analyse.txt'));
if isempty(fieldnames(P)), P = parametres_secours(da, d1); end
reglages = texte_reglages(P);

% étapes
n_k = numel(dir(fullfile(da, '*', 'choix_K.mat')));
n_cl = numel(dir(fullfile(da, '*', 'clonal', 'stats_clonales.mat')));
f_csv = fullfile(da, 'agregation', 'champ_globale', 'resultats_biais_clonal.csv');
e4 = isfile(f_csv);
etapes = sprintf('1:%d', n_etape1);
if n_k > 0, etapes = sprintf('%s 2:%d', etapes, n_k); end
if n_cl > 0, etapes = sprintf('%s 3:%d', etapes, n_cl); end
if e4, etapes = [etapes ' 4']; end

% sessions analysées : statut de l'étape 4, sinon choix du K
n_sess = sprintf('%d', max(n_k, 0));
f_st = fullfile(da, 'agregation', 'statut_sessions.csv');
if isfile(f_st)
    try
        S = readtable(f_st, 'TextType', 'string', 'Delimiter', ',');
        n_ok = sum(startsWith(S.statut, "ok"));
        n_sess = sprintf('%d/%d', n_ok, height(S));
    catch
    end
end

% dernier lancement
[fin, date_log] = dernier_lancement(da);
date_txt = date_log;
if isfield(P, 'ecrit_le') && ~isempty(P.ecrit_le), date_txt = P.ecrit_le; end
if isempty(date_txt)
    f = dir(da);
    if ~isempty(f), date_txt = datestr(max([f.datenum]), 'yyyy-mm-dd HH:MM'); end %#ok<DATST>
end

% biais clonal (test principal champ / globale)
n_asm = NaN; n_test = NaN; p05 = NaN; pm05 = NaN; q05 = NaN; q10 = NaN;
if e4
    try
        B = readtable(f_csv, 'TextType', 'string');
        v = B.Properties.VariableNames;
        n_asm = height(B);
        t = true(height(B), 1);
        if ismember('testable', v), t = logique(B.testable); end
        if ismember('exclue_criteres', v), t = t & ~logique(B.exclue_criteres); end
        n_test = sum(t);
        if n_test > 0
            if ismember('p_initial', v), p05 = 100 * mean(B.p_initial(t) < 0.05); end
            if ismember('p_minp', v), pm05 = 100 * mean(B.p_minp(t) < 0.05); end
            if ismember('q', v), q05 = sum(B.q(t) < 0.05); q10 = sum(B.q(t) < 0.10); end
        end
    catch
    end
end
ligne = {'', condition, date_txt, n1, texte_nom, version, reglages, etapes, n_sess, fin, ...
    NaN, n_asm, n_test, arrondi(p05), arrondi(pm05), q05, q10, d1, da};
end


function x = arrondi(x)
x = round(x * 10) / 10;
end


function t = logique(x)
if islogical(x), t = x; return; end
if isnumeric(x), t = x ~= 0 & ~isnan(x); return; end
t = ismember(lower(string(x)), ["1", "true", "oui"]);
end


function [fin, date_txt] = dernier_lancement(da)
% fin du dernier journal : fini, fini avec erreurs, arrêté, en cours
fin = '-'; date_txt = '';
Lg = dir(fullfile(da, 'log_analyse_*.txt'));
if isempty(Lg), return; end
[~, k] = max([Lg.datenum]);
f = fullfile(da, Lg(k).name);
date_txt = datestr(Lg(k).datenum, 'yyyy-mm-dd HH:MM'); %#ok<DATST>
try
    txt = fileread(f);
catch
    return
end
n_err = numel(regexp(txt, 'FAILED|ÉCHEC|échec', 'match'));
etapes = regexp(txt, '(?m)^(?:Step|Étape) (\d+b?) -', 'tokens');   % pas « Step 1 folder: »
etapes = unique(cellfun(@(c) c{1}, etapes, 'UniformOutput', false));
texte_etapes = '';
if ~isempty(etapes), texte_etapes = [' (steps ' strjoin(etapes, ',') ')']; end
if ~isempty(regexp(txt, '=== (END|FIN)', 'once'))
    if n_err > 0
        fin = sprintf('done, %d error(s)%s', n_err, texte_etapes);
    else
        fin = ['done' texte_etapes];
    end
elseif (now - Lg(k).datenum) < 15 / 1440 %#ok<TNOW1>
    fin = ['running?' texte_etapes];
else
    fin = ['stopped' texte_etapes];
end
end


function P = lire_parametres(f)
% lignes « clé = valeur » de parametres_analyse.txt (première occurrence)
P = struct();
if ~isfile(f), return; end
try
    L = splitlines(string(fileread(f)));
catch
    return
end
for i = 1:numel(L)
    tok = regexp(L(i), '^\s+([A-Za-z_][\w\.]*)\s+=\s+(.*?)\s*$', 'tokens', 'once');
    if isempty(tok), continue; end
    c = matlab.lang.makeValidName(strrep(char(tok(1)), '.', '_'));
    if ~isfield(P, c), P.(c) = char(tok(2)); end
end
ec = regexp(L, '^Écrit le\s*:\s*(.*)$', 'tokens', 'once');
ec = ec(~cellfun(@isempty, ec));
if ~isempty(ec), P.ecrit_le = char(ec{1}); end
end


function P = parametres_secours(da, d1)
% sans parametres_analyse.txt : opts de l'étape 1, config.asm, règle du K
P = struct();
try
    pre = dir(fullfile(d1, '*', 'pretraitement.mat'));
    S = load(fullfile(pre(1).folder, pre(1).name), 'opts');
    P = ajouter(P, S.opts);
catch
end
try
    fa = dir(fullfile(d1, '*', 'k*', 'assemblees.mat'));
    S = load(fullfile(fa(1).folder, fa(1).name), 'config_asm');
    P = ajouter(P, S.config_asm);
catch
end
try
    S = load(fullfile(da, 'best_k_par_session.mat'), 'regle_choix_K');
    r = S.regle_choix_K;
    if isfield(r, 'regle'), P.regle = r.regle; end
    if isfield(r, 'k_start'), P.k_start = num2str(r.k_start); end
    if isfield(r, 'K_fixe'), P.K_fixe = mat2str(r.K_fixe); end
    if isfield(r, 'liste_K'), P.liste_K = mat2str(r.liste_K); end
    if isfield(r, 'consensus') && isfield(r.consensus, 'seuil'), P.consensus_seuil = num2str(r.consensus.seuil); end
catch
end
end


function P = ajouter(P, s)
ch = fieldnames(s);
for i = 1:numel(ch)
    v = s.(ch{i});
    if ischar(v) || isstring(v), v = char(v);
    elseif (isnumeric(v) || islogical(v)) && numel(v) <= 30
        v = mat2str(v);
    else
        continue
    end
    if ~isfield(P, ch{i}), P.(ch{i}) = v; end
end
end


function t = texte_reglages(P)
% résumé court des réglages qui changent d'une analyse à l'autre
g = @(c, def) champ(P, c, def);
morceaux = {};
det = g('detection', '');
if contains(det, 'cascade')
    m = g('cascade_modele', '');
    tok = regexp(m, '(GC8m?)_EXC_(\d+(?:\.\d+)?Hz)_smoothing(\d+ms)', 'tokens', 'once');
    if ~isempty(tok), m = sprintf('%s %s %s', tok{:}); end
    morceaux{end+1} = ['CASCADE ' m];
else
    morceaux{end+1} = 'transients';
end
regle = g('sce_regle', '');
switch regle
    case 'fdr', morceaux{end+1} = ['SCE fdr ' g('sce_fdr', '?')];
    case 'percentile', morceaux{end+1} = ['SCE pct ' g('sce_percentile', '?')];
    case 'fixe', morceaux{end+1} = ['SCE ' g('sce_n_cells_threshold', '?') ' cells'];
end
mpd = g('MinPeakDistancesce', '');
if ~isempty(mpd), morceaux{end+1} = ['MPD ' mpd]; end
if strcmp(g('retirer_hyperactives', 'false'), 'true')
    h = 'hyperactive out';
    if ~isempty(g('hyperactives_seuil_hz', '')), h = [h ' ' g('hyperactives_seuil_hz', '') ' Hz']; end
    morceaux{end+1} = h;
end
cl = g('seuil_clusters', '');
if ~isempty(cl), morceaux{end+1} = ['clusters ' cl ' ' g('percentile_seuil', '')]; end
q = g('FDR_Q', '');
if ~isempty(q), morceaux{end+1} = ['FDR_Q ' q]; end
rk = g('regle', '');
kf = g('K_fixe', '');
lk = g('liste_K', '');
if ~isempty(kf) && ~strcmp(kf, '[]')
    morceaux{end+1} = ['K = ' kf];
elseif strcmp(rk, 'consensus')
    ks = g('k_start', '');
    plage = lk;
    if ~isempty(ks) && ~isempty(lk)
        try
            v = eval(lk); v = v(v >= str2double(ks)); %#ok<EVLDIR>
            plage = sprintf('%d-%d', min(v), max(v));
        catch
        end
    end
    morceaux{end+1} = sprintf('consensus %s K %s', g('consensus_seuil', ''), plage);
elseif ~isempty(rk)
    morceaux{end+1} = sprintf('K: %s %s', rk, lk);
end
mc = g('min_cellules_asm', '');
if ~isempty(mc), morceaux{end+1} = sprintf('>= %s cells, >= %s SCE', mc, g('min_SCE_asm', '?')); end
if strcmp(g('midp', ''), 'true'), morceaux{end+1} = 'mid-p'; end
t = strjoin(morceaux, ' | ');
end


function v = champ(P, c, def)
if isfield(P, c) && ~isempty(P.(c)), v = P.(c); else, v = def; end
end
