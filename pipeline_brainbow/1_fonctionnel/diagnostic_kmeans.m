function T = diagnostic_kmeans(dossier_session, K, liste_replicats, nSurrogates, liste_online)
%DIAGNOSTIC_KMEANS  Combien de réplicats k-means faut-il, et la phase en
%   ligne sert-elle à quelque chose ? Sur UNE session déjà prétraitée et UN K.
%
%   diagnostic_kmeans(dossier_session, 9)
%   diagnostic_kmeans(dossier_session, 9, [10 50 100], 20, [true false])
%
%   liste_online : phase en ligne testée (défaut [true false] : les deux).
%   Données : k-means avec chaque nombre de réplicats et 3 graines.
%     - inertie (somme des distances au centre : plus bas = meilleur optimum)
%     - silhouette moyenne et des 3 meilleurs clusters
%     - accord entre les 3 graines (indice de Rand ajusté, 1 = même partition) :
%       si l'accord est déjà proche de 1, plus de réplicats ne change rien
%     - accord_avec_online (lignes sans phase en ligne) : même partition
%       qu'avec la phase en ligne, même graine et mêmes réplicats ?
%   Surrogates : les MÊMES nSurrogates matrices mélangées, classées avec
%     chaque réglage. Si la silhouette de leurs meilleurs clusters ne monte
%     plus, des réplicats en plus ne changent pas le seuil.
%   Temps : durée mesurée, et estimation pour une session entière
%     (tous les K de config.liste_K, config.asm.kmeans_surrogate surrogates).
%   Aucune couleur n'est lue.

if nargin < 3 || isempty(liste_replicats), liste_replicats = [10 50 100 200]; end
if nargin < 4 || isempty(nSurrogates), nSurrogates = 20; end
if nargin < 5 || isempty(liste_online), liste_online = [true false]; end
config = config_brainbow();
MinSCE = config.asm.MinSCE;
null_kmeans = config.asm.null_kmeans;
if isempty(gcp('nocreate')), parpool; end
etat_rng = rng;
remettre_rng = onCleanup(@() rng(etat_rng)); %#ok<NASGU>

P = load(fullfile(dossier_session, 'pretraitement.mat'), 'Race');
Race = double(P.Race);
M = CovarMjc(Race);
fprintf('%s : %d cellules, %d SCE, K = %d, %d surrogates\n', dossier_session, ...
    size(Race, 1), size(Race, 2), K, nSurrogates);

nR = numel(liste_replicats);
nO = numel(liste_online);
nL = nR * nO;
col_online = false(nL, 1); col_rep = zeros(nL, 1);
inertie = nan(nL, 1); silh = nan(nL, 1); top3 = nan(nL, 1); accord = nan(nL, 1);
accord_on = nan(nL, 1);
t_donnees = nan(nL, 1); t_surr = nan(nL, 1);
sur_top = nan(nL, 3); sur_seuil = nan(nL, 1);
IDX_on = cell(nR, 3);        % partitions avec phase en ligne (pour accord_avec_online)
l = 0;
for io = 1:nO
    online = logical(liste_online(io));
    for ir = 1:nR
        l = l + 1;
        r = liste_replicats(ir);
        col_online(l) = online; col_rep(l) = r;
        % --- données, 3 graines
        IDXs = cell(1, 3); in_g = nan(1, 3); s_g = nan(1, 3); t_g = nan(1, 3);
        t0 = tic;
        for g = 1:3
            [IDXs{g}, sCl, ~, S, ~, sumd] = kmeansoptnew2(Race, r, 'var', K, M, MinSCE, g, online);
            in_g(g) = sum(sumd);
            s_g(g) = S;
            t_g(g) = mean(sCl(1:min(3, end)));
        end
        t_donnees(l) = toc(t0) / 3;
        inertie(l) = mean(in_g); silh(l) = mean(s_g); top3(l) = mean(t_g);
        accord(l) = mean([ari(IDXs{1}, IDXs{2}), ari(IDXs{1}, IDXs{3}), ari(IDXs{2}, IDXs{3})]);
        if online
            IDX_on(ir, :) = IDXs;
        elseif ~isempty(IDX_on{ir, 1})
            accord_on(l) = mean([ari(IDXs{1}, IDX_on{ir, 1}), ari(IDXs{2}, IDX_on{ir, 2}), ...
                ari(IDXs{3}, IDX_on{ir, 3})]);
        end

        % --- surrogates : mêmes matrices mélangées pour tous les réglages
        sCl_rnd = zeros(nSurrogates, K);
        t0 = tic;
        parfor i = 1:nSurrogates
            fixer_flux_aleatoire(5000, i, nSurrogates);
            sc = kmeansoptrndnew2(Race, r, K, MinSCE, null_kmeans, online);
            sCl_rnd(i, :) = sort(sc(:)', 'descend');
        end
        t_surr(l) = toc(t0);
        sur_top(l, :) = mean(sCl_rnd(:, 1:3), 1);
        sur_seuil(l) = prctile(sCl_rnd(:, 1), 95);
        fprintf('  phase en ligne %s, %d réplicats : fait (%.0f s)\n', onoff(online), r, ...
            t_donnees(l) * 3 + t_surr(l));
    end
end

texte_online = arrayfun(@onoff, col_online, 'UniformOutput', false);
T = table(texte_online, col_rep, inertie, silh, top3, accord, accord_on, ...
    sur_top(:, 1), sur_top(:, 2), sur_top(:, 3), sur_seuil, t_donnees, t_surr / nSurrogates, ...
    'VariableNames', {'online', 'replicats', 'inertie_donnees', 'silhouette_donnees', 'silh_3_meilleurs', ...
    'accord_graines_ARI', 'accord_avec_online', 'surr_silh_rang1', 'surr_silh_rang2', 'surr_silh_rang3', ...
    'surr_seuil95_rang1', 'temps_donnees_s', 'temps_par_surrogate_s'});
fprintf('\n');
disp(T)

% Estimation du temps d'une session (tous les K, kmeans_surrogate surrogates),
% même nombre de réplicats pour les données et les surrogates.
% Les surrogates tournent en parallèle : temps par surrogate mesuré ici.
nK = numel(config.liste_K);
nS = config.asm.kmeans_surrogate;
fprintf('Temps estimé pour une session (%d K, %d surrogates, mêmes réplicats partout), en minutes :\n', nK, nS);
for l = 1:nL
    fprintf('  phase en ligne %-3s, %4d réplicats : %6.0f min\n', T.online{l}, col_rep(l), ...
        nK * (t_donnees(l) + nS * t_surr(l) / nSurrogates) / 60);
end
fprintf(['\nLecture : sans phase en ligne, si inertie et silhouette sont presque les mêmes et ' ...
    'accord_avec_online proche de 1,\nla phase en ligne ne sert à rien. Réplicats : le plus petit ' ...
    'nombre où accord_graines_ARI et surr_silh_rang1 ne bougent plus.\n' ...
    'Garder le même réglage pour les données et les surrogates.\n']);
end


function t = onoff(x)
if x, t = 'on'; else, t = 'off'; end
end


function a = ari(x, y)
% Indice de Rand ajusté entre deux partitions
x = x(:); y = y(:);
[~, ~, ix] = unique(x);
[~, ~, iy] = unique(y);
C = accumarray([ix iy], 1);
n = numel(x);
somme = @(v) sum(v .* (v - 1) / 2);
sij = somme(C(:));
si = somme(sum(C, 2));
sj = somme(sum(C, 1));
attendu = si * sj / (n * (n - 1) / 2);
maxi = (si + sj) / 2;
a = (sij - attendu) / (maxi - attendu);
end
