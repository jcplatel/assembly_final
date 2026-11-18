function statut = etape3c_hierarchique(session, config)
%% =========================================================================
%  ÉTAPE 3c : clustering HIÉRARCHIQUE des cellules, puis couleur de chaque
%  groupe de l'arbre (demandé le 7 nov.)
%  -------------------------------------------------------------------------
%  Idée : partir de petits groupes de cellules qui s'activent ensemble et
%  monter vers des groupes de plus en plus gros ; jusqu'à quelle taille
%  trouve-t-on des groupes avec une couleur en excès ?
%  1. Arbres SANS couleur, sur les cellules actives dans au moins
%     opts.min_sce SCE (les mêmes cellules pour tous les arbres) :
%       'sce'  : co-participation aux SCE ; similarité = corrélation des
%                lignes de Race (phi, corrige le niveau d'activité de
%                chaque cellule) ; distance = 1 - phi ; lien moyen ;
%       'sttc' : STTC au repos (fenêtre de config.paires.dt_sttc_s), distance
%                et niveau de chaque cellule retirés (résidus de l'étape 3b) ;
%                distance = max - résidu ; paire non mesurée (ROI qui se
%                chevauchent) = résidu 0 ; lien moyen ;
%       'sttc_100ms' (ou toute autre fenêtre 'sttc_<n>ms') : pareil, fenêtre
%                de coïncidence +/- n ms (ajouté le 7 nov.) ;
%       'dff'  : corrélation de Pearson du dF/F au repos, résidus comme à
%                l'étape 3b (avec le couplage de chaque cellule au signal
%                commun) ; distance = max - résidu ; lien moyen (7 nov.) ;
%       'fca'  : Functional Clustering Algorithm (Feldt et al. 2009 ; Feldt
%                Muldoon et al. 2013) sur les événements au repos (fca_arbre) :
%                fusions par significativité contre des surrogates décalés,
%                coupe quand plus rien n'est significatif ; arbre continué
%                jusqu'à la racine pour la comparaison par taille.
%  2. (7 nov., demandé par l'utilisateur : peu de niveaux pour ne pas
%     multiplier les tests) Chaque arbre est COUPÉ à opts.niveaux : groupes
%     de 5, 10, 20, 40 cellules en moyenne (nombre de groupes = cellules /
%     niveau ; cluster 'maxclust'). À un niveau, les groupes sont séparés.
%  3. À chaque niveau (test_couleur_groupes) : couleur la plus en excès de
%     chaque groupe d'au moins opts.min_colorees cellules colorées (min des
%     p hypergéométriques) ; UN test du niveau : nombre de groupes avec
%     p < alpha contre les couleurs mélangées sur la MÊME coupe, dans tout
%     le champ (principal) et dans des carrés de taille_bloc_um (contrôle
%     spatial) ; puis quels groupes : p_minP (6 couleurs corrigées par les
%     mélanges) et FDR sur les groupes du niveau.
%     Pour la FCA, aussi ses clusters seuls ('fca_amas', niveau 0 : groupes
%     d'au moins 3 cellules formés avant la coupe de la FCA).
%  Écrit <session>\hierarchique\hierarchique.mat et arbres.png.
%  =========================================================================

opts = options_hierarchique(config);
dossier = fullfile(session.dossier, 'hierarchique');
if ~isfolder(dossier), mkdir(dossier); end
fichier = fullfile(dossier, 'hierarchique.mat');
if isfile(fichier), delete(fichier); end

%% 1. Activité (aucune couleur)
pre = load(fullfile(session.dossier_source, 'pretraitement.mat'), 'Race', 'Tr1b', 'Raster', 'WinRest', ...
    'sampling_rate');
Race = double(pre.Race) > 0;
NCell = size(Race, 1);
nSCE = size(Race, 2);
n_sce = sum(Race, 2);
repos = pre.WinRest(:)';
Rrepos = pre.Raster(:, repos) > 0;                % événements au repos (pour la FCA)
fs = pre.sampling_rate;

path = avec_separateur(session.dossier_suite2p);
[outline_gcampx, ~, gcamp_centroid, ~, ~, Coomasqx, Coomasqy] = charger_geometrie(path);
um = session.um_par_px;
if isnan(um), um = config.um_par_px_defaut; end
xy_um = gcamp_centroid * um;

% mesures de paires demandées comme arbres ('sttc', 'sttc_<n>ms', 'dff') :
% résidus comme à l'étape 3b (distance et niveau de chaque cellule retirés)
mesures_arbres = opts.arbres(~ismember(opts.arbres, {'sce', 'fca'}));
RES = struct();
VAL = struct();
if ~isempty(mesures_arbres)
    popts = config.paires;
    popts.mesures = unique([popts.mesures(:)', intersect(mesures_arbres, {'sttc', 'dff'})], 'stable');
    fen = regexp(mesures_arbres, '^sttc_(\d+)ms$', 'tokens', 'once');
    fen = fen(~cellfun(@isempty, fen));
    dts = cellfun(@(t) str2double(t{1}) / 1000, fen);
    if ~isfield(popts, 'dt_sensibilite_s'), popts.dt_sensibilite_s = []; end
    popts.dt_sensibilite_s = unique([popts.dt_sensibilite_s(:)', dts(:)']);
    S = similarite_paires(pre.Tr1b, pre.Raster, pre.WinRest, pre.sampling_rate, popts);
    dist_um = sqrt((xy_um(:, 1) - xy_um(:, 1)') .^ 2 + (xy_um(:, 2) - xy_um(:, 2)') .^ 2);
    V0 = ~masques_chevauchent(Coomasqx, Coomasqy) & isfinite(dist_um);
    V0(logical(eye(NCell))) = false;
    for k = 1:numel(mesures_arbres)
        nom = mesures_arbres{k};
        if ~isfield(S, nom)
            error('Arbre ''%s'' : mesure inconnue (''sttc'', ''sttc_<n>ms'' ou ''dff'').', nom);
        end
        VAL.(nom) = V0 & isfinite(S.(nom));
        RES.(nom) = residus_paires(S.(nom), VAL.(nom), dist_um, popts.pas_distance_um, strcmp(nom, 'dff'));
    end
    clear S
end
clear pre

%% 2. Couleurs (seule partie qui les lit)
colorcell = charger_couleurs(path, session.fichier_colorcell);
colorcell = couleurs_sans_doublons(colorcell, session);    % doublons de l'étape 1 : non colorés
verifier_alignement(NCell, colorcell, outline_gcampx, session.nom);
col = colorcell(:);

%% 3. Cellules des arbres : les mêmes pour tous les arbres
garder = n_sce >= opts.min_sce & n_sce < nSCE & sum(Rrepos, 2) >= 1;
for nom = fieldnames(VAL)'
    garder = garder & any(VAL.(nom{1}), 2);
end
idx = find(garder);
N = numel(idx);
H = struct('nom', session.nom, 'animal', session.id_animal, 'opts', opts, 'NCell', NCell, ...
    'nSCE', nSCE, 'cellules', idx, 'couleur', col(idx), 'xy_um', xy_um(idx, :));
H.n_colorees = sum(H.couleur >= 1 & H.couleur <= 6);
if N < 10 || H.n_colorees < opts.min_colorees
    statut = sprintf('%d cellules dans au moins %d SCE (%d colorées) : trop peu, session sautée', N, ...
        opts.min_sce, H.n_colorees);
    H.statut = statut;
    save(fichier, '-struct', 'H');
    return
end

% groupes de mélange des couleurs (les mêmes pour tous les arbres)
groupes = cell(1, numel(opts.permutations));
for ip = 1:numel(opts.permutations)
    switch opts.permutations{ip}
        case 'globale'
            groupes{ip} = ones(N, 1);
        case 'blocs'
            [~, ~, groupes{ip}] = unique(floor(H.xy_um / opts.taille_bloc_um), 'rows');
        otherwise
            error('config.hier.permutations : ''globale'' ou ''blocs'' (reçu ''%s'').', opts.permutations{ip});
    end
end

%% 4. Arbres, coupes et tests
lignes = cell(0, 10);
H.arbres = struct([]);
for a = 1:numel(opts.arbres)
    nom_arbre = opts.arbres{a};
    F = [];
    switch nom_arbre
        case 'fca'
            ev = arrayfun(@(u) find(Rrepos(u, :)), idx(:)', 'UniformOutput', false);
            F = fca_arbre(ev, size(Rrepos, 2), max(1, round(opts.fca_sigma_s * fs)), opts.fca_nsurr, true, ...
                session.graine + 1000);
            Z = F.Z;
        otherwise
            if strcmp(nom_arbre, 'sce')
                D = 1 - corr(double(Race(idx, :))');
            else                                    % 'sttc', 'sttc_<n>ms', 'dff' : résidus
                Es = RES.(nom_arbre)(idx, idx);
                Vs = VAL.(nom_arbre)(idx, idx);
                Es(~Vs) = 0;
                D = max(Es(Vs)) - Es;
            end
            D = (D + D') / 2;
            D(~isfinite(D)) = max(D(isfinite(D)));
            D(1:N + 1:end) = 0;
            Z = linkage(squareform(D, 'tovector'), opts.liaison);
    end
    [M, ordre, x_noeud, hauteur] = noeuds_arbre(Z);
    taille_noeud = sum(M, 2);
    A = struct('nom', nom_arbre, 'Z', Z, 'ordre', ordre, 'x', x_noeud, 'hauteur', hauteur, ...
        'coupes', struct([]), 'fca', [], 'amas', struct([]));
    for iv = 1:numel(opts.niveaux)
        niv = opts.niveaux(iv);
        k = max(1, min(N, round(N / niv)));
        lab = cluster(Z, 'maxclust', k);
        k = max(lab);
        G = false(k, N);
        G(sub2ind([k N], lab(:), (1:N)')) = true;
        C = struct('niveau', niv, 'k', k, 'lab', lab, 'hauteur', hauteur_coupe(Z, k), ...
            'noeud', noeud_de_chaque_groupe(G, M, taille_noeud), 'tests', struct([]));
        for ip = 1:numel(opts.permutations)
            R = test_couleur_groupes(G, H.couleur, groupes{ip}, opts, session.graine + 100 * a + 10 * iv + ip);
            R.permutation = opts.permutations{ip};
            if isempty(C.tests), C.tests = R; else, C.tests(end+1) = R; end
            lignes(end+1, :) = ligne_test(nom_arbre, niv, R); %#ok<AGROW>
        end
        if isempty(A.coupes), A.coupes = C; else, A.coupes(end+1) = C; end
    end
    if ~isempty(F)
        A.fca = rmfield(F, 'Z');
        if ~isempty(F.noeuds_amas)
            G = M(F.noeuds_amas, :);
            C = struct('niveau', 0, 'k', size(G, 1), 'lab', [], 'hauteur', F.coupe + 0.5, ...
                'noeud', F.noeuds_amas(:), 'tests', struct([]));
            for ip = 1:numel(opts.permutations)
                R = test_couleur_groupes(G, H.couleur, groupes{ip}, opts, session.graine + 900 + ip);
                R.permutation = opts.permutations{ip};
                if isempty(C.tests), C.tests = R; else, C.tests(end+1) = R; end
                lignes(end+1, :) = ligne_test('fca_amas', 0, R); %#ok<AGROW>
            end
            A.amas = C;
        end
    end
    if isempty(H.arbres), H.arbres = A; else, H.arbres(end+1) = A; end
end
H.T = cell2table(lignes, 'VariableNames', {'arbre', 'permutation', 'niveau', 'n_groupes', 'taille_mediane', ...
    'observes', 'attendus', 'p', 'n_fdr', 'n_coupe'});
H.T.animal = repmat({session.id_animal}, height(H.T), 1);
H.T.session = repmat({session.nom}, height(H.T), 1);

t = H.T(strcmp(H.T.permutation, opts.permutations{1}), :);
morceaux = arrayfun(@(i) sprintf('%s ~%g : %d/%.1f (FDR %d)', t.arbre{i}, t.niveau(i), t.observes(i), ...
    t.attendus(i), t.n_fdr(i)), 1:height(t), 'UniformOutput', false);
statut = sprintf('%d cellules (%d colorées) ; groupes avec une couleur en excès, observés/attendus : %s', ...
    N, H.n_colorees, strjoin(morceaux, ' ; '));
i_fca = find(strcmp({H.arbres.nom}, 'fca'), 1);
if ~isempty(i_fca)
    statut = sprintf('%s ; FCA : %d fusions significatives, %d clusters', statut, H.arbres(i_fca).fca.coupe, ...
        numel(H.arbres(i_fca).fca.amas));
end
H.statut = statut;
save(fichier, '-struct', 'H');

if ~isfield(config, 'faire_figures') || config.faire_figures
    try
        figure_arbre_session(dossier);
    catch err
        warning('figure_arbre_session : %s', err.message);
    end
end
end


function L = ligne_test(nom, niv, R)
L = {nom, R.permutation, niv, R.n_groupes, median(R.m(R.teste)), R.obs, R.attendu, R.p, R.n_fdr, ...
    numel(R.teste)};
end


function h = hauteur_coupe(Z, k)
% hauteur où l'arbre est coupé en k groupes (entre deux fusions)
N = size(Z, 1) + 1;
i = N - k;                     % nombre de fusions faites sous la coupe
if i < 1
    h = Z(1, 3) / 2;
elseif i >= N - 1
    h = Z(end, 3);
else
    h = (Z(i, 3) + Z(i + 1, 3)) / 2;
end
end


function noeud = noeud_de_chaque_groupe(G, M, taille_noeud)
% noeud de l'arbre qui correspond à chaque groupe (0 = cellule seule)
noeud = zeros(size(G, 1), 1);
for c = 1:size(G, 1)
    nc = sum(G(c, :));
    if nc < 2, continue; end
    cand = find(taille_noeud == nc);
    for i = cand(:)'
        if isequal(M(i, :), G(c, :))
            noeud(c) = i;
            break
        end
    end
end
end
