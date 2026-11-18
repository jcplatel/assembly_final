function lignes = parametres_analyse(sessions, config, fichier)
%PARAMETRES_ANALYSE  Tous les paramètres d'une analyse, de l'étape 1 à l'étape 4,
%   LUS DANS LES FICHIERS DE RÉSULTATS : ce qui a vraiment été utilisé pour
%   chaque session, pas seulement ce que dit config_brainbow.m aujourd'hui.
%
%   lignes = parametres_analyse(sessions, config)
%   parametres_analyse(sessions, config, 'E:\...\parametres_analyse.txt')
%
%   sessions : table des sessions (lire_table_sessions : dossier, dossier_source)
%   config   : avec les dossiers de l'analyse (dossiers_analyse)
%   Un paramètre qui n'a pas la même valeur pour toutes les sessions est
%   marqué « ≠ » avec chaque valeur et son nombre de sessions.
%   Écrit aussi par main_brainbow à la fin de chaque lancement
%   (<dossier de l'analyse>\parametres_analyse.txt) ; bouton « Analysis
%   parameters » de l'interface.

config = dossiers_analyse(config);
L = {};
L{end+1} = '=================================================================';
L{end+1} = 'PARAMÈTRES DE L''ANALYSE (lus dans les fichiers de résultats)';
L{end+1} = '=================================================================';
L{end+1} = sprintf('Version du code          : %s (lancée le %s)', version_pipeline(), datestr(now, 'yyyy-mm-dd HH:MM'));
v_asm = 'v1 (une assemblée par cellule)';
if isfield(config, 'version_assemblees') && strcmpi(config.version_assemblees, 'v2'), v_asm = 'v2 (cellules partagées)'; end
L{end+1} = sprintf('Version des assemblées   : %s', v_asm);
L{end+1} = sprintf('Dossier de l''étape 1     : %s', config.PathEtape1);
if isempty(config.analyse)
    L{end+1} = 'Analyse (étapes 2 à 4)   : main (dans le dossier de l''étape 1)';
else
    L{end+1} = sprintf('Analyse (étapes 2 à 4)   : %s (%s)', config.analyse, config.PathSave);
end
L{end+1} = sprintf('Table des sessions       : %s', config.fichier_sessions);
cond = 'toutes';
if isfield(config, 'conditions') && ~isempty(config.conditions), cond = strjoin(cellstr(config.conditions), ', '); end
L{end+1} = sprintf('Sessions retenues        : %d (inclure = 1 ; conditions : %s)', height(sessions), cond);
L{end+1} = sprintf('Écrit le                 : %s', datestr(now, 'yyyy-mm-dd HH:MM'));

%% Étape 1 : prétraitement, détection des SCE, assemblées
nS = height(sessions);
pre_opts = cell(nS, 1); asm_opts = cell(nS, 1);
seuil = nan(nS, 1); nsce = nan(nS, 1); Kmin = nan(nS, 1); Kmax = nan(nS, 1); dates1 = nan(nS, 1);
for s = 1:nS
    ds = sessions.dossier_source{s};
    f = fullfile(ds, 'pretraitement.mat');
    if ~isfile(f), continue; end
    try
        noms = intersect({'opts_utilisees', 'opts', 'TRace'}, who('-file', f));
        P = load(f, noms{:});
        if isfield(P, 'opts_utilisees'), pre_opts{s} = P.opts_utilisees; end
        if isfield(P, 'opts') && isfield(P.opts, 'sce_n_cells_threshold_calcule')
            seuil(s) = P.opts.sce_n_cells_threshold_calcule;
        end
        if isfield(P, 'TRace'), nsce(s) = numel(P.TRace); end
        d = dir(f); dates1(s) = d.datenum;
        Lk = dir(fullfile(ds, 'k*'));
        Ks = str2double(regexprep({Lk([Lk.isdir]).name}, '^k', ''));
        Ks = sort(Ks(~isnan(Ks)));
        if ~isempty(Ks)
            Kmin(s) = Ks(1); Kmax(s) = Ks(end);
            fk = fullfile(ds, sprintf('k%d', Ks(1)), 'assemblees.mat');
            if isfile(fk) && ismember('config_asm', who('-file', fk))
                A = load(fk, 'config_asm');
                asm_opts{s} = A.config_asm;
            end
        end
    catch err
        L{end+1} = sprintf('   (session %s illisible : %s)', sessions.nom{s}, err.message); %#ok<AGROW>
    end
end
faits = ~cellfun(@isempty, pre_opts);
L{end+1} = '';
L{end+1} = sprintf('--- Étape 1 : prétraitement et seuil des SCE (%d sessions calculées sur %d%s) ---', ...
    sum(faits), nS, texte_dates(dates1));
L = [L, comparer(pre_opts(faits), sessions.nom(faits))];
L{end+1} = sprintf('   seuil des SCE réellement utilisé : %s cellules ; nombre de SCE : %s', ...
    med_min_max(seuil), med_min_max(nsce));
L{end+1} = '';
L{end+1} = sprintf('--- Étape 1 : assemblées (k-means, surrogates, RACEORTHO) ; K calculés : %s ---', ...
    texte_K(Kmin, Kmax));
a_faits = ~cellfun(@isempty, asm_opts);
L = [L, comparer(asm_opts(a_faits), sessions.nom(a_faits))];

%% Étape 2 : choix du K
L{end+1} = '';
f2 = fullfile(config.PathSave, 'best_k_par_session.mat');
if isfile(f2)
    B = load(f2);
    d = dir(f2);
    L{end+1} = sprintf('--- Étape 2 : choix du K (%s) ---', datestr(d.datenum, 'yyyy-mm-dd HH:MM'));
    if isfield(B, 'regle_choix_K'), L = [L, comparer({B.regle_choix_K}, {'-'})]; end
    if isfield(B, 'best_k') && istable(B.best_k)
        Kr = B.best_k.K(strcmp(B.best_k.statut, 'ok'));
        if any(Kr == 0)
            L{end+1} = sprintf('   résultat : consensus entre les K pour %d sessions ; %d exclues', ...
                sum(Kr == 0), sum(~strcmp(B.best_k.statut, 'ok')));
        else
            L{end+1} = sprintf('   résultat : K retenus %s pour %d sessions ; %d exclues', med_min_max(Kr), ...
                numel(Kr), sum(~strcmp(B.best_k.statut, 'ok')));
        end
    end
else
    L{end+1} = '--- Étape 2 : pas encore faite dans cette analyse ---';
end

%% Étape 3 : statistiques clonales
L{end+1} = '';
oc = {}; noms_oc = {}; ids = {}; variantes = {}; dates3 = nan(nS, 1);
for s = 1:nS
    f = fullfile(sessions.dossier{s}, 'clonal', 'stats_clonales.mat');
    if ~isfile(f), continue; end
    try
        noms = intersect({'opts_clonal', 'id_execution', 'noms_variantes'}, who('-file', f));
        C = load(f, noms{:});
        d = dir(f); dates3(s) = d.datenum;
        if isfield(C, 'id_execution'), ids{end+1} = char(string(C.id_execution)); end %#ok<AGROW>
        if isfield(C, 'noms_variantes') && isempty(variantes), variantes = C.noms_variantes; end
        if isfield(C, 'opts_clonal'), oc{end+1} = C.opts_clonal; noms_oc{end+1} = sessions.nom{s}; end %#ok<AGROW>
    catch
    end
end
n3 = sum(~isnan(dates3));
if n3 == 0
    L{end+1} = '--- Étape 3 : pas encore faite dans cette analyse ---';
else
    L{end+1} = sprintf('--- Étape 3 : statistiques clonales (%d sessions%s ; %d exécution(s)) ---', ...
        n3, texte_dates(dates3), numel(unique(ids)));
    if ~isempty(variantes)
        L{end+1} = sprintf('   variantes : %s (la première = test principal)', strjoin(cellstr(variantes), ', '));
    end
    if isempty(oc)
        L{end+1} = '   réglages non enregistrés (étape 3 faite avant le 4 nov. : la refaire pour les voir)';
    else
        L = [L, comparer(oc, noms_oc)];
    end
end

%% Étape 4 : agrégation
L{end+1} = '';
f4 = fullfile(config.PathSave, 'agregation', 'parametres_etape4.mat');
if isfile(f4)
    G = load(f4);
    d = dir(f4);
    L{end+1} = sprintf('--- Étape 4 : agrégation (%s) ---', datestr(d.datenum, 'yyyy-mm-dd HH:MM'));
    if isfield(G, 'agreg'), L = [L, comparer({G.agreg}, {'-'})]; end
elseif isfolder(fullfile(config.PathSave, 'agregation'))
    L{end+1} = '--- Étape 4 : faite, mais réglages non enregistrés (avant le 4 nov. : la refaire pour les voir) ---';
else
    L{end+1} = '--- Étape 4 : pas encore faite dans cette analyse ---';
end

lignes = L;
if nargin >= 3 && ~isempty(fichier)
    fid = fopen(fichier, 'w', 'n', 'UTF-8');
    if fid >= 0
        fprintf(fid, '%s\n', L{:});
        fclose(fid);
    end
end
end


function L = comparer(structs, noms_sessions)
% Une ligne par paramètre ; « ≠ » si les sessions n'ont pas toutes la même valeur
L = {};
if isempty(structs)
    L{end+1} = '   (aucun fichier)';
    return
end
vals = cellfun(@aplatir, structs, 'UniformOutput', false);
champs = {};
for i = 1:numel(vals), champs = union(champs, vals{i}(:, 1), 'stable'); end
for c = 1:numel(champs)
    v = cell(numel(vals), 1);
    for i = 1:numel(vals)
        k = find(strcmp(vals{i}(:, 1), champs{c}), 1);
        if isempty(k), v{i} = '(absent)'; else, v{i} = vals{i}{k, 2}; end
    end
    [u, ~, ic] = unique(v, 'stable');
    if numel(u) == 1
        L{end+1} = sprintf('   %-34s = %s', champs{c}, u{1}); %#ok<AGROW>
    else
        morceaux = arrayfun(@(k) sprintf('%s (%d, ex. %s)', u{k}, sum(ic == k), ...
            noms_sessions{find(ic == k, 1)}), 1:numel(u), 'UniformOutput', false);
        L{end+1} = sprintf(' ≠ %-34s : %s', champs{c}, strjoin(morceaux, ' | ')); %#ok<AGROW>
    end
end
end


function P = aplatir(s, prefixe)
% struct -> {nom, valeur en texte} (sous-structures en a.b ; fonctions ignorées)
if nargin < 2, prefixe = ''; end
P = cell(0, 2);
noms = fieldnames(s);
for i = 1:numel(noms)
    v = s.(noms{i});
    nom = [prefixe noms{i}];
    if isa(v, 'function_handle'), continue; end
    if isstruct(v) && isscalar(v)
        P = [P; aplatir(v, [nom '.'])]; %#ok<AGROW>
    else
        P(end+1, :) = {nom, en_texte(v)}; %#ok<AGROW>
    end
end
end


function t = en_texte(v)
if ischar(v)
    t = v;
elseif isstring(v)
    t = char(strjoin(v, ', '));
elseif islogical(v) && isscalar(v)
    if v, t = 'true'; else, t = 'false'; end
elseif isnumeric(v) || islogical(v)
    v = double(v);
    if isempty(v)
        t = '[]';
    elseif isscalar(v)
        t = sprintf('%g', v);
    elseif isvector(v) && numel(v) > 3 && all(v == round(v)) && numel(unique(diff(v))) == 1
        pas = v(2) - v(1);
        if pas == 1, t = sprintf('%g:%g', v(1), v(end)); else, t = sprintf('%g:%g:%g', v(1), pas, v(end)); end
    else
        t = mat2str(v, 6);
    end
elseif iscell(v)
    if iscellstr(v) && size(v, 2) == 2 && size(v, 1) > 1
        t = strjoin(strcat(v(:, 1), '/', v(:, 2)), ', ');      % ex. variantes {référence, permutation}
    elseif iscellstr(v)
        t = ['{' strjoin(v(:)', ', ') '}'];
    else
        t = sprintf('{cellule %dx%d}', size(v, 1), size(v, 2));
    end
elseif isstruct(v)
    t = sprintf('[struct %dx%d]', size(v, 1), size(v, 2));
else
    t = sprintf('[%s]', class(v));
end
if numel(t) > 140, t = [t(1:137) '...']; end
end


function t = med_min_max(x)
x = x(~isnan(x));
if isempty(x), t = '-'; return; end
t = sprintf('médiane %g [%g - %g]', median(x), min(x), max(x));
end


function t = texte_K(Kmin, Kmax)
ok = ~isnan(Kmin);
if ~any(ok), t = '-'; return; end
paires = unique([Kmin(ok), Kmax(ok)], 'rows');
if size(paires, 1) == 1
    t = sprintf('%d à %d (toutes les sessions)', paires(1, 1), paires(1, 2));
else
    t = sprintf('≠ selon les sessions (%s)', strjoin(arrayfun(@(i) sprintf('%d-%d', paires(i, 1), ...
        paires(i, 2)), 1:size(paires, 1), 'UniformOutput', false), ', '));
end
end


function t = texte_dates(d)
d = d(~isnan(d));
if isempty(d), t = ''; return; end
if max(d) - min(d) < 1
    t = sprintf(', le %s', datestr(min(d), 'yyyy-mm-dd'));
else
    t = sprintf(', du %s au %s', datestr(min(d), 'yyyy-mm-dd'), datestr(max(d), 'yyyy-mm-dd'));
end
end
