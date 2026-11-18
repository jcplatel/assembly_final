function statut = etape3b_paires(session, config)
%% =========================================================================
%  ÉTAPE 3b : les cellules voisines de même couleur sont-elles plus
%  synchrones au repos que les voisines de couleurs différentes ?
%  -------------------------------------------------------------------------
%  N'utilise que le prétraitement (étape 1) : ni assemblées ni choix du K.
%  1. Activité seule (aucune couleur) : similarité de toutes les paires au
%     repos (STTC des événements, corrélation du ΔF/F), puis on retire ce qui
%     ne dépend pas de la lignée : distance, niveau propre de chaque cellule
%     (residus_paires).
%  2. Couleurs : même couleur contre couleurs différentes, par tranche de
%     distance et par couleur, contre des permutations des couleurs.
%  Écrit <session>\paires\paires.mat (res : une ligne par mesure x
%  permutation x couleur x tranche) et paires.png.
%  =========================================================================

opts = config.paires;
dossier = fullfile(session.dossier, 'paires');
if ~isfolder(dossier), mkdir(dossier); end
fichier = fullfile(dossier, 'paires.mat');
if isfile(fichier), delete(fichier); end

%% 1. Activité au repos (aucune couleur)
pre = load(fullfile(session.dossier_source, 'pretraitement.mat'), 'Tr1b', 'Raster', 'WinRest', 'sampling_rate', 'Race');
NCell = size(pre.Raster, 1);
S = similarite_paires(pre.Tr1b, pre.Raster, pre.WinRest, pre.sampling_rate, opts, pre.Race);
clear pre

%% 2. Géométrie : distances en µm, ROI qui se chevauchent exclues
path = avec_separateur(session.dossier_suite2p);
[outline_gcampx, ~, gcamp_centroid, ~, ~, Coomasqx, Coomasqy] = charger_geometrie(path);
um_par_px = session.um_par_px;
if isnan(um_par_px), um_par_px = config.um_par_px_defaut; end
xy_um = gcamp_centroid * um_par_px;
dist_um = sqrt((xy_um(:, 1) - xy_um(:, 1)') .^ 2 + (xy_um(:, 2) - xy_um(:, 2)') .^ 2);
chevauche = masques_chevauchent(Coomasqx, Coomasqy);
V0 = ~chevauche & isfinite(dist_um);
V0(logical(eye(NCell))) = false;

%% 3. Couleurs (seule partie qui les lit)
[colorcell, ~, source_couleur] = charger_couleurs(path, session.fichier_colorcell);
colorcell = couleurs_sans_doublons(colorcell, session);    % doublons de l'étape 1 : non colorés
verifier_alignement(NCell, colorcell, outline_gcampx, session.nom);
colorcell = colorcell(:);
colorees = colorcell >= 1 & colorcell <= 6;

% Diagnostic : les couleurs ont-elles la même activité au repos ?
taux_min = S.n_evenements / (S.duree_repos_s / 60);       % événements par minute de repos
p_taux_couleur = NaN;
if numel(unique(colorcell(colorees))) >= 2
    p_taux_couleur = kruskalwallis(taux_min(colorees), colorcell(colorees), 'off');
end

% Diagnostic : certaines couleurs sont-elles plus couplées au réseau ?
% couplage = STTC moyen de chaque cellule avec toutes les autres (sans couleur),
% puis comparaison entre couleurs (z dans la session, moyenne par couleur)
Sp = S.(opts.mesures{1});
Sp(~V0) = NaN;
couplage = mean(Sp, 2, 'omitnan');
couplage_z = (couplage - mean(couplage, 'omitnan')) / std(couplage, 'omitnan');
couplage_par_couleur = nan(1, 6);
for c = 1:6
    couplage_par_couleur(c) = mean(couplage_z(colorcell == c), 'omitnan');
end
p_couplage_couleur = NaN;
ok_c = colorees & ~isnan(couplage);
if numel(unique(colorcell(ok_c))) >= 2
    p_couplage_couleur = kruskalwallis(couplage(ok_c), colorcell(ok_c), 'off');
end
clear Sp

%% 4. Tests
res = {};
mesures = liste_mesures_paires(opts);
[bornes, cumulatif] = distances_paires(opts);     % disques (défaut) ou tranches
for im = 1:numel(mesures)
    mesure = mesures{im};
    C = S.(mesure);
    V = V0 & isfinite(C);
    E = residus_paires(C, V, dist_um, opts.pas_distance_um, strcmp(mesure, 'dff'));
    for ip = 1:numel(opts.permutations)
        mode_perm = opts.permutations{ip};
        groupes = groupes_permutation(mode_perm, xy_um, taux_min, couplage, opts);
        T = test_couleur_paires(E, V, dist_um, colorcell, bornes, opts.nPerm, groupes, ...
            opts.min_paires, session.graine + 100 * im + ip, cumulatif);
        n = height(T);
        % similarité brute moyenne de toutes les paires à cette distance (pour
        % exprimer delta en % : delta / similarité moyenne)
        T.similarite_moyenne = nan(n, 1);
        for i = 1:n
            switch T.forme{i}
                case 'disque'
                    m = V & dist_um <= T.dist_max_um(i);
                case 'anneau'
                    m = V & dist_um > T.dist_min_um(i) & dist_um <= T.dist_max_um(i);
                otherwise
                    m = V & dist_um >= T.dist_min_um(i) & dist_um < T.dist_max_um(i);
            end
            T.similarite_moyenne(i) = mean(C(m), 'omitnan');
        end
        T = [table(repmat({session.id_animal}, n, 1), repmat({session.id_champ}, n, 1), ...
            repmat({session.id_session}, n, 1), repmat({session.id_condition}, n, 1), ...
            repmat({mesure}, n, 1), repmat({mode_perm}, n, 1), ...
            'VariableNames', {'animal', 'champ', 'session', 'condition', 'mesure', 'permutation'}), T]; %#ok<AGROW>
        res{end+1, 1} = T; %#ok<AGROW>
    end
end
res = vertcat(res{:});
statut = 'ok';
config_paires = opts;
dt_trames = S.dt_trames;
save(fichier, 'res', 'statut', 'taux_min', 'p_taux_couleur', 'source_couleur', 'config_paires', ...
    'um_par_px', 'dt_trames', 'couplage', 'couplage_par_couleur', 'p_couplage_couleur');

%% 5. Résumé et figure
principal = strcmp(res.mesure, opts.mesures{1}) & strcmp(res.permutation, opts.permutations{1}) ...
    & strcmp(res.couleur, 'toutes') & ~strcmp(res.forme, 'anneau');
R1 = res(principal, :);
fprintf('   [%s] paires (%s, %s), toutes couleurs : ', session.nom, opts.mesures{1}, opts.permutations{1});
for i = 1:height(R1)
    fprintf('%g-%g µm : delta %+.4f p = %.3f | ', R1.dist_min_um(i), R1.dist_max_um(i), R1.delta(i), R1.p(i));
end
fprintf('\n');
if config.faire_figures
    figure_paires(res(~strcmp(res.forme, 'anneau'), :), opts, fullfile(dossier, 'paires.png'), session.nom);
end
end


function groupes = groupes_permutation(mode_perm, xy_um, taux_min, couplage, opts)
% Groupes à l'intérieur desquels les couleurs sont permutées
N = size(xy_um, 1);
switch mode_perm
    case 'globale'
        groupes = ones(N, 1);
    case 'blocs'        % carrés de taille_bloc_um : garde la densité locale de chaque couleur
        cases = floor(xy_um / opts.taille_bloc_um);
        [~, ~, groupes] = unique(cases, 'rows');
    case 'activite'     % 10 niveaux d'activité : garde le lien éventuel couleur / activité
        q = quantile(taux_min, 0.1:0.1:0.9);
        groupes = discretize(taux_min, unique([-Inf, q(:)', Inf]));
    case 'couplage'     % 20 niveaux de couplage au réseau (STTC moyen avec toutes les autres) :
                        % garde le lien éventuel couleur / couplage individuel
        q = quantile(couplage, (1:19) / 20);
        groupes = discretize(couplage, unique([-Inf, q(:)', Inf]));
    otherwise
        error('Permutation inconnue : %s (''globale'', ''blocs'', ''couplage'' ou ''activite'').', mode_perm);
end
groupes = groupes(:);
groupes(isnan(groupes)) = 0;          % cellules sans valeur : un groupe à part
end
