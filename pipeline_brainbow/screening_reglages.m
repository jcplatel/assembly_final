function S = screening_reglages(config, opts)
%% =========================================================================
%  SCREENING_REGLAGES : essayer beaucoup de réglages sur quelques animaux et
%  chercher ceux qui donnent le plus d'assemblées biaisées en couleur
%  (13 nov., demande de l'utilisateur ; réglages à valider ENSUITE sur les
%  autres animaux et sur les sessions avec indice)
%  -------------------------------------------------------------------------
%  Trois étapes, toutes reprises là où elles se sont arrêtées si on relance :
%
%  1. 'etape1'    pour chaque réglage du PRÉTRAITEMENT (config.pre : seuil
%                 des transients, seuil et largeur des SCE, lissage...), un
%                 seul changé à la fois : étape 1 complète (détection, SCE,
%                 k-means et surrogates pour chaque K) des sessions du
%                 screening. C'EST LA PARTIE LONGUE.
%                 -> <dossier>\etape1\<réglage>\<session>\  (dossiers d'étape 1
%                 ordinaires : utilisables ensuite comme config.PathSave)
%  2. 'variantes' pour chaque étape 1, toutes les variantes d'assemblées,
%                 SANS refaire le k-means : clusters retenus (rang / max /
%                 tous, percentile), règle des membres (FDR_Q...), version 1
%                 ou 2, chaque K, la règle 'exces', les consensus.
%                 -> <dossier>\variantes\<réglage>\<session>.mat
%                 Aucune couleur n'est lue jusqu'ici.
%  3. 'agregation' lit les couleurs, calcule le biais clonal de chaque
%                 combinaison (même test que l'étape 3, champ_globale), les
%                 classe, et mesure ce que le hasard donnerait après le même
%                 tri (couleurs mélangées, screening refait).
%                 -> <dossier>\resultats\ : resume_screening.txt (à lire en
%                 premier), figure_screening.png, screening_top.csv,
%                 screening_facteurs.csv, screening_combinaisons.csv
%
%  Lancer :
%     screening_reglages                       % config_brainbow.m, 3 animaux par défaut
%     config = config_brainbow();  ...         % (PathSave = dossier d'étape 1 de la base)
%     screening_reglages(config)
%     screening_reglages(config, opts)
%
%  opts (tous facultatifs) :
%     animaux        {'446178', '444119', '444175'}   animaux du screening
%     dossier        <dossier parent de PathSave>\screening
%     etapes         {'etape1', 'variantes', 'agregation'} ; ex. {'agregation'}
%                    pour refaire seulement le classement (rapide) ;
%                    {'plan'} : ne calcule rien, affiche les sessions, les
%                    réglages et le nombre d'étapes 1 qui restent à faire
%     facteurs       {} = tous ; ex. {'sce_seuil', 'threshold_peak'} : seulement
%                    ces réglages du prétraitement (+ la base)
%     grille         grille modifiée (G = grille_screening(config), puis changer G)
%     pre_supplementaires   réglages de prétraitement en plus : struct array
%                    (nom, facteur, niveau, modif), ex. une combinaison de
%                    plusieurs changements (voir resume_screening.txt, point 5)
%     reutiliser_base true : si l'étape 1 de config.PathSave a exactement les
%                    réglages de base pour ces sessions et ces K, elle est lue
%                    telle quelle (rien n'y est écrit) au lieu d'être refaite
%     NShuf          config.asm.NShuf : permutations des membres et des activations
%     nNull          500 : mélanges des couleurs pour mesurer l'effet du tri
%     nTirages       20000 : tirages de la loi de p_minp (par champ)
%     min_testables  12 : une combinaison est classée si elle a au moins
%     min_par_animal 3    12 assemblées testables, dont 3 par animal (avec 1
%                    ou 2 assemblées, l'excès d'un animal vaut 0 ou +95 points :
%                    trop de bruit)
%
%  Les réglages essayés sont dans screening\grille_screening.m.
%  Durée : (nombre de réglages de prétraitement) x (durée de l'étape 1 des
%  sessions du screening). Le temps restant est affiché après chaque réglage.
%  On peut arrêter (Ctrl+C) et relancer ; et lancer {'agregation'} à tout
%  moment sur les réglages déjà finis.
%  =========================================================================

if nargin < 1 || isempty(config), config = config_brainbow(); end
if nargin < 2 || isempty(opts), opts = struct(); end
addpath(genpath(fileparts(mfilename('fullpath'))));

% L'étape 1 de la base est dans config.PathSave (jamais dans une analyse)
if isfield(config, 'PathEtape1') && ~isempty(config.PathEtape1)
    config.PathSave = config.PathEtape1;
    config = rmfield(config, 'PathEtape1');
    if isfield(config, 'PathSource'), config = rmfield(config, 'PathSource'); end
end
config.analyse = '';
config.version_assemblees = 'v1';
config.seulement = {};

racine_defaut = fullfile(fileparts(regexprep(char(config.PathSave), '[\\/]+$', '')), 'screening');
opts = defaut(opts, 'animaux', {'446178', '444119', '444175'});
opts = defaut(opts, 'dossier', racine_defaut);
opts = defaut(opts, 'etapes', {'etape1', 'variantes', 'agregation'});
opts = defaut(opts, 'facteurs', {});
opts = defaut(opts, 'grille', []);
opts = defaut(opts, 'pre_supplementaires', []);
opts = defaut(opts, 'reutiliser_base', true);
opts = defaut(opts, 'NShuf', config.asm.NShuf);
opts = defaut(opts, 'nNull', 500);
opts = defaut(opts, 'nTirages', 20000);
opts = defaut(opts, 'min_testables', 12);
opts = defaut(opts, 'min_par_animal', 3);
opts = defaut(opts, 'sessions', []);
opts.animaux = cellstr(opts.animaux);
opts.etapes = cellstr(opts.etapes);
opts.facteurs = cellstr(opts.facteurs);
racine = char(opts.dossier);
if ~isfolder(racine), mkdir(racine); end
fichier_journal = fullfile(racine, 'journal_screening.txt');

%% ---------- Sessions du screening ------------------------------------------
if ~isempty(opts.sessions)
    sessions = opts.sessions;                    % (tests) struct array déjà prêt
else
    T = lire_table_sessions(dossiers_analyse(config));
    T = T(ismember(T.id_animal, opts.animaux), :);
    if isempty(T)
        error('Screening : aucune session de la table pour les animaux %s.', strjoin(opts.animaux, ', '));
    end
    sessions = table2struct(T);
end
nS = numel(sessions);
manque = setdiff(opts.animaux, {sessions.id_animal});
if ~isempty(manque) && isempty(opts.sessions)
    warning('Screening : aucune session pour %s.', strjoin(manque, ', '));
end

%% ---------- Grille ----------------------------------------------------------
if isempty(opts.grille), G = grille_screening(config); else, G = opts.grille; end
if isfield(G.pre, 'dossier'), G.pre = rmfield(G.pre, 'dossier'); end
if ~isempty(opts.facteurs)
    G.pre = G.pre(strcmp({G.pre.facteur}, 'base') | ismember({G.pre.facteur}, opts.facteurs));
end
% réglages ajoutés avant : gardés d'un lancement à l'autre (fichier pre_supplementaires.mat)
fichier_sup = fullfile(racine, 'pre_supplementaires.mat');
sup = struct('nom', {}, 'facteur', {}, 'niveau', {}, 'modif', {});
if isfile(fichier_sup)
    Ssup = load(fichier_sup, 'sup');
    sup = Ssup.sup;
end
for i = 1:numel(opts.pre_supplementaires)
    x = opts.pre_supplementaires(i);
    x.nom = regexprep(char(x.nom), '[^A-Za-z0-9_\-]', 'p');
    sup(strcmp({sup.nom}, x.nom)) = [];
    sup(end+1) = struct('nom', x.nom, 'facteur', char(x.facteur), 'niveau', char(x.niveau), 'modif', x.modif); %#ok<AGROW>
end
if ~isempty(sup), save(fichier_sup, 'sup'); end
for i = 1:numel(sup)
    if any(strcmp({G.pre.nom}, sup(i).nom)), continue; end
    G.pre(end+1) = sup(i);
end
nP = numel(G.pre);
% Une seule population de cellules par screening (config.pre.population) :
% les couleurs sont mélangées parmi les cellules gardées, qui doivent être
% les mêmes pour tous les réglages comparés.
pops = cell(1, nP);
for p = 1:nP
    pre_p = fusion_struct(config.pre, G.pre(p).modif);
    pops{p} = 'toutes';
    if isfield(pre_p, 'population') && ~isempty(pre_p.population), pops{p} = char(pre_p.population); end
end
if numel(unique(pops)) > 1
    error(['Screening : les réglages n''ont pas tous la même population de cellules (%s). Faire un ' ...
        'screening par population (un dossier opts.dossier chacun), avec config.pre.population dans la base.'], ...
        strjoin(unique(pops), ', '));
end
asm = config.asm;
asm.NShuf = opts.NShuf;

% Dossier d'étape 1 de chaque réglage
[G.pre.dossier] = deal('');
for p = 1:nP
    G.pre(p).dossier = fullfile(racine, 'etape1', G.pre(p).nom);
end
i_base = find(strcmp({G.pre.nom}, 'base'), 1);
if opts.reutiliser_base && ~isempty(i_base)
    [ok, raison] = base_deja_faite(config, sessions, G.liste_K);
    if ok
        G.pre(i_base).dossier = char(config.PathSave);
        journal(fichier_journal, sprintf('Base : étape 1 déjà faite, lue dans %s (rien n''y est écrit).', ...
            G.pre(i_base).dossier));
    elseif ~isfolder(fullfile(G.pre(i_base).dossier, sessions(1).nom))
        journal(fichier_journal, sprintf(['Base : l''étape 1 de %s n''est pas réutilisable (%s) : ' ...
            'refaite dans le dossier du screening.'], char(config.PathSave), raison));
    end
end
save(fullfile(racine, 'grille.mat'), 'G', 'config');
journal(fichier_journal, sprintf(['=== SCREENING (pipeline %s) : %d sessions (%s), %d réglages de prétraitement, ' ...
    'K = %s ==='], version_pipeline(), nS, strjoin(unique({sessions.id_animal}), ', '), nP, mat2str(G.liste_K)));
journal(fichier_journal, sprintf(['Par prétraitement : %d règles de clusters x %d règles de membres x %d versions ' ...
    'x %d choix du K x %d critères = %d combinaisons'], size(G.clusters, 1), numel(G.membres), numel(G.versions), ...
    numel(G.choix_K), numel(G.min_cellules) * numel(G.min_SCE), size(G.clusters, 1) * numel(G.membres) * ...
    numel(G.versions) * numel(G.choix_K) * numel(G.min_cellules) * numel(G.min_SCE)));

% Ce qui reste à calculer à l'étape 1 (un réglage x une session = une « étape 1 de session »)
n_reste = 0;
for p = 1:nP
    if strcmp(G.pre(p).dossier, char(config.PathSave)), continue; end
    cfg = config;
    cfg.pre = fusion_struct(config.pre, G.pre(p).modif);
    cfg.liste_K = G.liste_K;
    for s = 1:nS
        n_reste = n_reste + ~etape1_a_jour(fullfile(G.pre(p).dossier, sessions(s).nom), cfg);
    end
end
journal(fichier_journal, sprintf(['Étape 1 : %d étape(s) 1 de session restent à calculer sur %d (%d réglages x %d ' ...
    'sessions). Durée = ce nombre x la durée de l''étape 1 d''une session (voir « Step 1 - ... done » dans ' ...
    'un log_analyse).'], n_reste, nP * nS, nP, nS));
if ismember('plan', opts.etapes)
    fprintf('Sessions : %s\n', strjoin({sessions.nom}, ', '));
    fprintf('Réglages du prétraitement : %s\n', strjoin({G.pre.nom}, ', '));
    fprintf('Dossier du screening : %s\n', racine);
    S = [];
    return
end

if any(ismember({'etape1', 'variantes'}, opts.etapes))
    try
        if isempty(gcp('nocreate')), parpool; end
    catch
    end
end

%% ---------- 1. Étape 1 de chaque réglage du prétraitement ------------------
if ismember('etape1', opts.etapes)
    durees = [];
    for p = 1:nP
        if strcmp(G.pre(p).dossier, char(config.PathSave)), continue; end   % base déjà faite
        cfg = config;
        cfg.pre = fusion_struct(config.pre, G.pre(p).modif);
        cfg.liste_K = G.liste_K;
        cfg.faire_figures = false;
        cfg.executer_scripts_anciens = false;
        cfg.forcer_pretraitement = false;
        cfg.forcer_assemblees = false;
        t0 = tic;
        calcule = false;
        n_echecs = 0;
        for s = 1:nS
            sess = sessions(s);
            sess.dossier = fullfile(G.pre(p).dossier, sess.nom);
            sess.dossier_source = sess.dossier;
            if etape1_a_jour(sess.dossier, cfg), continue; end
            calcule = true;
            ts = tic;
            try
                etape1_fonctionnel(sess, cfg);
                journal(fichier_journal, sprintf('Étape 1 - %s (%d/%d) - %s : faite (%s)', G.pre(p).nom, p, nP, ...
                    sess.nom, duree(toc(ts))));
            catch err
                n_echecs = n_echecs + 1;
                journal(fichier_journal, sprintf('Étape 1 - %s (%d/%d) - %s : ÉCHEC : %s', G.pre(p).nom, p, nP, ...
                    sess.nom, err.message));
            end
        end
        if calcule && n_echecs > 0
            journal(fichier_journal, sprintf(['Étape 1 - %s : %d session(s) en échec ; ce réglage ne sera pas ' ...
                'classé tant qu''elles manquent (relancer pour réessayer).'], G.pre(p).nom, n_echecs));
        elseif calcule
            durees(end+1) = toc(t0); %#ok<AGROW>
            reste = 0;
            for q = p + 1:nP
                for s = 1:nS
                    cfgq = cfg;
                    cfgq.pre = fusion_struct(config.pre, G.pre(q).modif);
                    if ~etape1_a_jour(fullfile(G.pre(q).dossier, sessions(s).nom), cfgq)
                        reste = reste + 1 / nS;
                    end
                end
            end
            journal(fichier_journal, sprintf(['Étape 1 - %s finie en %s. Reste environ %.1f réglage(s) de ' ...
                'prétraitement : %s.'], G.pre(p).nom, duree(durees(end)), reste, duree(reste * median(durees))));
        end
    end
end

%% ---------- 2. Variantes d'assemblées (sans couleur) -----------------------
if ismember('variantes', opts.etapes)
    for p = 1:nP
        cfg = config;
        cfg.pre = fusion_struct(config.pre, G.pre(p).modif);
        cfg.liste_K = G.liste_K;
        d_var = fullfile(racine, 'variantes', G.pre(p).nom);
        if ~isfolder(d_var), mkdir(d_var); end
        for s = 1:nS
            sess = sessions(s);
            d_s = fullfile(G.pre(p).dossier, sess.nom);
            if ~etape1_a_jour(d_s, cfg)
                journal(fichier_journal, sprintf('Variantes - %s - %s : étape 1 absente ou incomplète, sautée.', ...
                    G.pre(p).nom, sess.nom));
                continue
            end
            f = fullfile(d_var, [sess.nom '.mat']);
            info = dir(fullfile(d_s, 'pretraitement.mat'));
            signature = struct();
            signature.liste_K = G.liste_K;
            signature.clusters = G.clusters;
            signature.membres = G.membres;
            signature.versions = G.versions;
            signature.choix_K = G.choix_K;
            signature.asm = asm;
            signature.date_etape1 = info.datenum;
            if isfile(f) && memes_parametres(f, 'signature', signature), continue; end
            ts = tic;
            try
                R = screening_session(d_s, sess, G, asm);
                save(f, 'R', 'signature');
                journal(fichier_journal, sprintf('Variantes - %s (%d/%d) - %s : %d jeux d''assemblées (%s)', ...
                    G.pre(p).nom, p, nP, sess.nom, numel(R.ens), duree(toc(ts))));
            catch err
                if isfile(f), delete(f); end
                journal(fichier_journal, sprintf('Variantes - %s (%d/%d) - %s : ÉCHEC : %s', G.pre(p).nom, p, nP, ...
                    sess.nom, err.message));
            end
        end
    end
end

%% ---------- 3. Agrégation (couleurs) ----------------------------------------
S = [];
if ismember('agregation', opts.etapes)
    etat_rng = rng;
    remettre = onCleanup(@() rng(etat_rng)); %#ok<NASGU>
    S = screening_agregation(racine, G, sessions, config, opts);
end
end


%% =========================================================================
function opts = defaut(opts, champ, valeur)
if ~isfield(opts, champ) || (isempty(opts.(champ)) && ~isempty(valeur))
    opts.(champ) = valeur;
end
end


function s = fusion_struct(s, modif)
noms = fieldnames(modif);
for k = 1:numel(noms)
    s.(noms{k}) = modif.(noms{k});
end
end


function ok = etape1_a_jour(dossier_s, cfg)
% Étape 1 de cette session finie avec exactement ces réglages, pour tous les K
ok = false;
f = fullfile(dossier_s, 'pretraitement.mat');
if ~isfile(f) || ~memes_parametres(f, 'opts_utilisees', cfg.pre), return; end
for K = cfg.liste_K
    fK = fullfile(dossier_s, sprintf('k%d', K), 'assemblees.mat');
    if ~isfile(fK) || ~memes_parametres(fK, 'config_asm', cfg.asm), return; end
end
ok = true;
end


function [ok, raison] = base_deja_faite(config, sessions, liste_K)
% L'étape 1 de config.PathSave a-t-elle exactement les réglages de base
% pour toutes les sessions du screening et tous les K ? (sinon on ne la
% touche pas : elle serait refaite, donc écrasée)
ok = true; raison = '';
cfg = config;
cfg.liste_K = liste_K;
for s = 1:numel(sessions)
    d_s = fullfile(char(config.PathSave), sessions(s).nom);
    if ~etape1_a_jour(d_s, cfg)
        ok = false;
        if ~isfile(fullfile(d_s, 'pretraitement.mat'))
            raison = sprintf('%s : pas de prétraitement', sessions(s).nom);
        else
            raison = sprintf('%s : réglages différents de config, ou K manquants', sessions(s).nom);
        end
        return
    end
end
end


function journal(fichier, msg)
ligne = sprintf('[%s] %s', datestr(now, 'dd/mm HH:MM'), msg);
fprintf('%s\n', ligne);
try
    fid = fopen(fichier, 'a', 'n', 'UTF-8');
    fprintf(fid, '%s\n', ligne);
    fclose(fid);
catch
end
end


function t = duree(secondes)
h = floor(secondes / 3600);
m = floor((secondes - 3600 * h) / 60);
if h > 0
    t = sprintf('%d h %02d min', h, m);
elseif m > 0
    t = sprintf('%d min %02d s', m, round(secondes - 60 * m));
else
    t = sprintf('%.0f s', secondes);
end
end
