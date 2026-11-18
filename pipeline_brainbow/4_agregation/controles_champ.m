function lignes = controles_champ(contenus, sessions, config, dossier_agreg)
%CONTROLES_CHAMP  Contrôles du champ, sans les assemblées (étape 4).
%   1. Les couleurs sont-elles regroupées dans l'espace ?
%      Par champ (regroupement_couleurs) : distance moyenne des paires de même
%      couleur / paires de couleurs différentes, comparée au mélange des
%      couleurs ; indice 1 = hasard, < 1 = même couleur plus proche.
%      Toutes couleurs ensemble, puis chaque couleur.
%   2. La couleur change-t-elle la détection ou l'activité ?
%      Pour chaque couleur : bruit (seuil de détection), F0, nombre de
%      transients au repos, participation aux SCE, part de membres
%      d'assemblées ; chaque fois en rapport à la moyenne des cellules
%      colorées du même champ (1 = pareil). Et les cellules hyperactives.
%   Unité : l'animal. Chaque session est son propre champ de vue (un même
%   plan revu un autre jour n'est en général pas à la même profondeur) :
%   sessions moyennées par animal ; Wilcoxon sur les animaux.
%   Pas de correction pour les comparaisons multiples dans la partie 2 :
%   c'est un contrôle (on cherche un biais), ne pas corriger est le choix
%   prudent.
%
%   contenus : stats_clonales.mat des sessions analysées (variable
%              cellules, écrite par l'étape 3 depuis le 5 nov.)
%   sessions : les lignes correspondantes de la table des sessions
%   Écrit <dossier_agreg>\controles_champ\ : regroupement_couleurs.csv,
%   detection_couleurs.csv, controles_champ.mat, et les figures
%   (figures_controles_champ). lignes : texte de resume_resultats.txt.

noms_coul = {'R', 'G', 'B', 'RG', 'RB', 'GB'};
nC = 6;
nmin = config.clonal.nmin_couleur;
nshuffle = config.clonal.nshuffle_distance;
lignes = {};

avec = cellfun(@(C) isfield(C, 'cellules'), contenus(:));
if ~any(avec)
    lignes = {['Contrôles du champ : relancer l''étape 3 (les données par cellule ' ...
        'qu''ils utilisent sont écrites depuis le 5 nov.).']};
    fprintf('\n%s\n', lignes{1});
    return
end
if ~all(avec)
    warning(['%d session(s) sans données par cellule (étape 3 d''avant le 5 nov.) : ' ...
        'ignorée(s) dans les contrôles du champ.'], sum(~avec));
end
contenus = contenus(avec);
sessions = sessions(avec, :);
dossier = fullfile(dossier_agreg, 'controles_champ');
if ~isfolder(dossier), mkdir(dossier); end
nS = numel(contenus);
nom    = cellstr(string(sessions.nom));
animal = cellstr(string(sessions.id_animal));
% Chaque session = un champ de vue (5 nov. : pas de regroupement par id_champ,
% un même plan revu un autre jour est en général à une autre profondeur)
champ  = nom;
if ismember('graine', sessions.Properties.VariableNames)
    graine = sessions.graine;
else
    graine = (1:nS)';
end

%% 1. Couleurs regroupées dans l'espace
% 12 nov. (demande de l'utilisateur) : mesure principale = paires de même
% couleur à moins de config.paires.rayon_principal_um (100 µm), observées /
% attendues après mélange des couleurs (> 1 = regroupées), la même mesure que
% les panneaux D et E ; l'indice de distance (< 1 = regroupées) est gardé
% dans le CSV et une ligne du résumé.
rayon = 100;
if isfield(config, 'paires') && isfield(config.paires, 'rayon_principal_um') && ...
        isfinite(config.paires.rayon_principal_um)
    rayon = config.paires.rayon_principal_um;
end
fprintf('\nContrôles du champ : couleurs regroupées dans l''espace (%d mélanges par session)...\n', nshuffle);
G = cell(nS, 1);
for s = 1:nS
    c = contenus{s}.cellules;
    G{s} = regroupement_couleurs(c.xy_um, c.couleur, nshuffle, graine(s) + 300, rayon);
end
Tsp = table(nom, animal, champ, cellfun(@(g) g.n_colorees, G), cellfun(@(g) g.indice, G), ...
    cellfun(@(g) g.z, G), cellfun(@(g) g.p, G), ...
    'VariableNames', {'session', 'animal', 'champ', 'n_colorees', 'indice', 'z', 'p'});
Tsp.obs_attendu = cellfun(@(g) g.obs_attendu, G);
Tsp.z_rayon     = cellfun(@(g) g.z_rayon, G);
Tsp.p_rayon     = cellfun(@(g) g.p_rayon, G);
Tsp.n_paires_rayon = cellfun(@(g) g.n_paires_rayon, G);
for k = 1:nC
    nk = cellfun(@(g) g.n_couleur(k), G);
    ik = cellfun(@(g) g.indice_c(k), G);
    pk = cellfun(@(g) g.p_c(k), G);
    ok_k = cellfun(@(g) g.obs_attendu_c(k), G);
    pok_k = cellfun(@(g) g.p_rayon_c(k), G);
    ik(nk < nmin) = NaN;              % couleur trop rare dans ce champ
    pk(nk < nmin) = NaN;
    ok_k(nk < nmin) = NaN;
    pok_k(nk < nmin) = NaN;
    Tsp.(['n_' noms_coul{k}]) = nk;
    Tsp.(['indice_' noms_coul{k}]) = ik;
    Tsp.(['p_' noms_coul{k}]) = pk;
    Tsp.(['obs_attendu_' noms_coul{k}]) = ok_k;
    Tsp.(['p_rayon_' noms_coul{k}]) = pok_k;
end
writetable(Tsp, fullfile(dossier, 'regroupement_couleurs.csv'));

% Toutes couleurs : champ (moyenne des sessions), puis animal
% (esp.* = mesure principale, le comptage à moins de rayon ; esp.dist_* =
% l'indice de distance d'avant le 12 nov.)
esp = struct();
esp.mesure = 'comptage';
esp.rayon_um = rayon;
[esp.animal, esp.animaux, esp.champ, esp.champs, esp.animal_du_champ] = ...
    par_animal(log(Tsp.obs_attendu), animal, champ);
[~, ~, z_c, champs_z] = par_animal(Tsp.z_rayon, animal, champ);
[trouve, loc] = ismember(esp.champs, champs_z);
z_champ = nan(numel(esp.champs), 1);
z_champ(trouve) = z_c(loc(trouve));
esp.z_champ = z_champ;
esp.n_champs_sig = sum(z_champ > 1.6449);           % z de la session : p < 0,05 unilatéral
esp.p_wilcoxon = test_animaux(esp.animal);
esp.indice_median = exp(median(esp.animal));
[esp.dist_animal] = par_animal(log(Tsp.indice), animal, champ);
esp.dist_p_wilcoxon = test_animaux(esp.dist_animal);
esp.dist_indice_median = exp(median(esp.dist_animal));
% Par couleur
esp.c_animal = cell(1, nC); esp.c_animaux = cell(1, nC);
esp.c_indice_median = nan(1, nC); esp.c_p = nan(1, nC); esp.c_n_champs = zeros(1, nC);
for k = 1:nC
    [esp.c_animal{k}, esp.c_animaux{k}, vc] = par_animal(log(Tsp.(['obs_attendu_' noms_coul{k}])), animal, champ);
    esp.c_n_champs(k) = numel(vc);
    if ~isempty(esp.c_animal{k})
        esp.c_indice_median(k) = exp(median(esp.c_animal{k}));
        esp.c_p(k) = test_animaux(esp.c_animal{k});
    end
end

%% 2. Détection et activité selon la couleur
mesures = {'bruit', 'F0', 'n_transients', 'n_SCE', 'membre'};
noms_mesures = {'noise (detection threshold)', 'F0', 'transients at rest', ...
    'SCE participation', 'assembly members'};
noms_mesures_fr = {'bruit', 'F0', 'transients au repos', 'participation aux SCE', 'membres d''assemblées'};
nM = numel(mesures);
lignes_det = cell(0, 5);                 % session, animal, champ, mesure, rapports (1 x 6)
rapports = nan(nS, nC, nM);
hyper_obs = zeros(1, nC); n_coul_tot = zeros(1, nC); n_hyper_tot = 0; n_cell_tot = 0; n_hyper_col = 0;
for s = 1:nS
    c = contenus{s}.cellules;
    col = c.couleur(:);
    est_col = col >= 1 & col <= nC;
    n_cell_tot = n_cell_tot + numel(col);
    if isfield(c, 'hyperactive')
        n_hyper_tot = n_hyper_tot + sum(c.hyperactive);
        n_hyper_col = n_hyper_col + sum(c.hyperactive(est_col));
    end
    for k = 1:nC
        mk = col == k;
        n_coul_tot(k) = n_coul_tot(k) + sum(mk);
        if isfield(c, 'hyperactive')
            hyper_obs(k) = hyper_obs(k) + sum(c.hyperactive(mk));
        end
    end
    for m = 1:nM
        v = double(c.(mesures{m})(:));
        for k = 1:nC
            mk = col == k & isfinite(v);
            ma = est_col & isfinite(v);
            if sum(mk) < nmin || sum(ma) == 0, continue; end
            if strcmp(mesures{m}, 'membre')
                % part de membres (pseudo-compte : une couleur sans membre reste finie)
                r = ((sum(v(mk)) + 0.5) / (sum(mk) + 1)) / ((sum(v(ma)) + 0.5) / (sum(ma) + 1));
            else
                r = mean(v(mk)) / mean(v(ma));
            end
            if r > 0 && isfinite(r)
                rapports(s, k, m) = r;
            end
        end
        lignes_det(end+1, :) = {nom{s}, animal{s}, champ{s}, mesures{m}, squeeze(rapports(s, :, m))}; %#ok<AGROW>
    end
end
Tdet = cell2table(lignes_det(:, 1:4), 'VariableNames', {'session', 'animal', 'champ', 'mesure'});
R6 = vertcat(lignes_det{:, 5});
for k = 1:nC
    Tdet.(['rapport_' noms_coul{k}]) = R6(:, k);
end
writetable(Tdet, fullfile(dossier, 'detection_couleurs.csv'));

det = struct();
det.mesures = mesures;
det.noms = noms_mesures;
det.animal = cell(nM, nC); det.animaux = cell(nM, nC);
det.rapport_median = nan(nM, nC); det.p = nan(nM, nC);
for m = 1:nM
    for k = 1:nC
        [det.animal{m, k}, det.animaux{m, k}] = par_animal(log2(rapports(:, k, m)), animal, champ);
        if ~isempty(det.animal{m, k})
            det.rapport_median(m, k) = 2 ^ median(det.animal{m, k});
            det.p(m, k) = test_animaux(det.animal{m, k});
        end
    end
end
part_coul = n_coul_tot / max(sum(n_coul_tot), 1);
det.hyper_obs = hyper_obs;
det.hyper_attendu = n_hyper_col * part_coul;
det.n_hyper_tot = n_hyper_tot;
det.n_cell_tot = n_cell_tot;

animaux_tous = unique(animal);
save(fullfile(dossier, 'controles_champ.mat'), 'Tsp', 'Tdet', 'esp', 'det', 'noms_coul', 'animaux_tous');

%% Texte du résumé
nA = numel(esp.animaux);
lignes{end+1} = sprintf('Contrôles du champ (sans assemblées ; %d animaux, %d sessions, chacune son champ de vue) :', ...
    nA, numel(esp.champs));
lignes{end+1} = sprintf(['   Couleurs regroupées dans l''espace : à moins de %g µm, %s %% de paires de même ' ...
    'couleur en plus que le hasard (couleurs mélangées dans le champ ; médiane des animaux ; Wilcoxon sur ' ...
    'les animaux : p = %s) ; %d/%d sessions à p < 0,05.'], esp.rayon_um, ...
    nb_txt(100 * (esp.indice_median - 1), '%+.0f'), pv_txt(esp.p_wilcoxon), esp.n_champs_sig, numel(esp.champs));
t = {};
for k = 1:nC
    if isnan(esp.c_indice_median(k)), continue; end
    t{end+1} = sprintf('%s %s %% (p = %s, %d sessions)', noms_coul{k}, ...
        nb_txt(100 * (esp.c_indice_median(k) - 1), '%+.0f'), pv_txt(esp.c_p(k)), esp.c_n_champs(k)); %#ok<AGROW>
end
if ~isempty(t)
    lignes{end+1} = sprintf('   Par couleur (paires c-c à moins de %g µm, en plus du hasard) : %s', esp.rayon_um, strjoin(t, ' ; '));
end
lignes{end+1} = sprintf(['   (Indice de distance, toutes distances : paires de même couleur %s %% plus proches que ' ...
    'les paires de couleurs différentes, p = %s.)'], nb_txt(100 * (1 - esp.dist_indice_median), '%.0f'), ...
    pv_txt(esp.dist_p_wilcoxon));
lignes{end+1} = ['   Détection selon la couleur (moyenne de la couleur / moyenne des cellules colorées ' ...
    'du champ ; Wilcoxon sur les animaux, non corrigé : pour un contrôle, c''est le choix prudent) :'];
for m = 1:nM
    sig = find(det.p(m, :) < 0.05);
    if isempty(sig)
        lignes{end+1} = sprintf('      %s : aucune couleur ne diffère', noms_mesures_fr{m}); %#ok<AGROW>
    else
        t = arrayfun(@(k) sprintf('%s x%s (p = %s)', noms_coul{k}, ...
            nb_txt(det.rapport_median(m, k), '%.2f'), pv_txt(det.p(m, k))), sig, 'UniformOutput', false);
        lignes{end+1} = sprintf('      %s : %s', noms_mesures_fr{m}, strjoin(t, ', ')); %#ok<AGROW>
    end
end
if n_hyper_tot > 0
    t = arrayfun(@(k) sprintf('%s %d (%s attendues)', noms_coul{k}, hyper_obs(k), ...
        nb_txt(det.hyper_attendu(k), '%.1f')), 1:nC, 'UniformOutput', false);
    lignes{end+1} = sprintf('   Cellules hyperactives : %d (%s %% des cellules) ; colorées : %s', ...
        n_hyper_tot, nb_txt(100 * n_hyper_tot / max(n_cell_tot, 1), '%.1f'), strjoin(t, ', '));
else
    lignes{end+1} = '   Cellules hyperactives : aucune (ou non enregistrées dans pretraitement.mat)';
end
fprintf('%s\n', lignes{:});

%% Figures
try
    figures_controles_champ(dossier);
catch err
    warning('figures_controles_champ non exécutée : %s', err.message);
end
end


function [v_animal, animaux, v_champ, champs, animal_du_champ] = par_animal(v, animal, champ)
% Sessions d'un même champ moyennées, puis champs moyennés par animal
v = v(:);
ok = isfinite(v);
v_animal = zeros(0, 1); animaux = cell(0, 1); v_champ = zeros(0, 1);
champs = cell(0, 1); animal_du_champ = cell(0, 1);
if ~any(ok), return; end
idx = find(ok);
[champs, ~, gc] = unique(champ(idx));
v_champ = accumarray(gc, v(idx), [], @mean);
premier = accumarray(gc, idx, [], @min);
animal_du_champ = animal(premier);
animal_du_champ = animal_du_champ(:);
[animaux, ~, ga] = unique(animal_du_champ);
v_animal = accumarray(ga, v_champ, [], @mean);
end


function p = test_animaux(x)
% Wilcoxon signé bilatéral sur les animaux (x : une valeur par animal, 0 = hasard)
x = x(isfinite(x));
if numel(x) < 2
    p = NaN;
else
    p = signrank(x);
end
end


function t = nb_txt(x, fmt)
if isnan(x)
    t = 'NA';
else
    t = strrep(sprintf(fmt, x), '.', ',');
end
end


function t = pv_txt(p)
if isnan(p)
    t = 'NA';
elseif p < 0.001
    t = '< 0,001';
else
    t = strrep(sprintf('%.3f', p), '.', ',');
end
end
