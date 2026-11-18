function test_pipeline_synthetique()
%% =========================================================================
%  TEST DU PIPELINE SUR DES DONNÉES SIMULÉES — à lancer AVANT les vraies données
%  -------------------------------------------------------------------------
%  Crée 3 animaux (1 session chacun) dans un dossier temporaire :
%    - 120 cellules (+ 10 ROI non-cellules dans Fall.mat), 5 min à 30 Hz,
%      4 épisodes de course
%    - 4 assemblées implantées de 15 cellules, 30 co-activations chacune au repos
%    - assemblées 1 et 2 « clonales » : 9 de leurs 15 cellules ont la même couleur
%  Lance les 4 étapes avec des paramètres allégés, puis vérifie :
%    1. que les assemblées retrouvées correspondent aux assemblées implantées
%    2. que les assemblées clonales ont un p_minp petit, les autres non
%  Durée : quelques minutes. Les anciens scripts ne sont pas lancés.
%  =========================================================================

racine = fullfile(tempdir, 'test_brainbow');
if isfolder(racine), rmdir(racine, 's'); end
mkdir(racine);
rng(0, 'twister');    % générateur imposé : mêmes données simulées à chaque lancement

animaux = {'900001', '900002', '900003'};
lignes = cell(numel(animaux), 9);
verite = cell(numel(animaux), 1);
for a = 1:numel(animaux)
    [chemin_data, dossier_s2p, verite{a}] = simuler_session(racine, animaux{a});
    lignes(a, :) = {1, animaux{a}, 'plane0', '260101', 'sans_indice', chemin_data, dossier_s2p, '', NaN};
end
T = cell2table(lignes, 'VariableNames', {'inclure', 'id_animal', 'id_champ', 'id_session', ...
    'id_condition', 'fichier_data', 'dossier_suite2p', 'fichier_colorcell', 'um_par_px'});

% Paramètres allégés
config = config_brainbow();
config.fichier_sessions = fullfile(racine, 'sessions_test.xlsx');
writetable(T, config.fichier_sessions);
config.PathSave              = fullfile(racine, 'resultats');
config.conditions            = {};
config.pre.motion_correction = false;
config.pre.sce_regle         = 'percentile';   % seuil calculé (120 cellules : un seuil fixe à 14 serait trop haut)
config.pre.sce_percentile    = 99;
config.pre.sce_n_cells_threshold = 5;
config.liste_K               = 3:6;
config.choixK.k_start        = 3;
config.asm.kmean_iter        = 10;
config.asm.kmeans_surrogate  = 20;
config.asm.NShuf             = 1000;
config.clonal.nPerm          = 2000;
config.agreg.min_cellules_asm = 6;
config.agreg.min_SCE_asm      = 3;
config.executer_scripts_anciens = false;

main_brainbow(config);

%% Vérifications
fprintf('\n================ VÉRIFICATIONS ================\n');
for a = 1:numel(animaux)
    nom = sprintf('%s_260101_plane0_sans_indice', animaux{a});
    dossier = fullfile(config.PathSave, nom);
    choix = load(fullfile(dossier, 'choix_K.mat'));
    if ~strcmp(choix.statut_choix, 'ok')
        fprintf('Animal %s : %s\n', animaux{a}, choix.statut_choix);
        continue
    end
    A = load(fullfile(dossier, sprintf('k%d', choix.best_K), 'assemblees.mat'), 'assemblyortho');
    C = load(fullfile(dossier, 'clonal', 'stats_clonales.mat'));
    fprintf('\nAnimal %s : K retenu = %d, %d assemblées retrouvées (4 implantées)\n', ...
        animaux{a}, choix.best_K, numel(A.assemblyortho));
    for i = 1:4
        vraie = verite{a}.assemblees{i};
        jac = cellfun(@(x) numel(intersect(x, vraie)) / numel(union(x, vraie)), A.assemblyortho);
        if isempty(jac)
            fprintf('  implantée %d : non retrouvée\n', i);
            continue
        end
        [jmax, meilleure] = max(jac);
        if verite{a}.clonale(i), type = 'clonale'; else, type = 'non clonale'; end
        fprintf('  implantée %d (%s) : Jaccard = %.2f avec l''assemblée %d, p_minp = %.4f\n', ...
            i, type, jmax, meilleure, C.res.p_minp(meilleure));
    end
end
fprintf(['\nAttendu : Jaccard proche de 1 ; p_minp petit (< 0.01) pour les clonales, ' ...
    'quelconque pour les autres.\nRésultats complets dans %s\n'], config.PathSave);

%% Étape 3b : synchronie des paires selon la couleur
% Les 9 cellules clonales des assemblées 1 (RG) et 2 (RB) s'activent ensemble
% mais sont dispersées dans le champ : attendu delta > 0 pour RG et RB à
% toutes les distances, et autour de 0 pour les autres couleurs.
config.paires.nPerm = 200;
main_paires(config);
end


%% =========================================================================
function [chemin_data, dossier_s2p, verite] = simuler_session(racine, id_animal)
fs = 30;  Nz = 300 * fs;  NCell = 120;  n_non_cellules = 10;
dossier_s2p = fullfile(racine, ['suite2p_' id_animal]);
mkdir(dossier_s2p);

% Vitesse : 4 épisodes de course de 10 s
speed = zeros(1, Nz);
for e = 1:4
    t0 = round((e * 70 - 40) * fs);
    speed(t0:t0 + 10 * fs) = 5;
end
% Instants de repos sûrs : loin de la course et des bords
loin_course = ~logical(movmax(double(speed > 0), [5 * fs 5 * fs]));
loin_course([1:10 * fs, end - 10 * fs + 1:end]) = false;
candidats = find(loin_course);
candidats = candidats(1:15:end);                       % >= 15 trames entre deux évènements
temps = candidats(randperm(numel(candidats), 120));    % 30 évènements x 4 assemblées

% Potentiels d'action : fond aléatoire (0.02 Hz) + assemblées
spikes = rand(NCell, Nz) < 0.02 / fs;
membres = randperm(NCell, 60);
assemblees = cell(1, 4);
for k = 1:4
    assemblees{k} = sort(membres((k - 1) * 15 + (1:15)));
    for t = temps((k - 1) * 30 + (1:30))
        actives = assemblees{k}(rand(1, 15) < 0.7);
        decalage = randi([0 1], 1, numel(actives));
        spikes(sub2ind([NCell Nz], actives, t + decalage)) = true;
    end
end

% Fluorescence : ligne de base lente + transients (décroissance 0.5 s) + bruit
tau = 0.5 * fs;
noyau = exp(-(0:round(5 * tau)) / tau);
calcium_signal = conv2(double(spikes), noyau);
calcium_signal = calcium_signal(:, 1:Nz);
t_s = (1:Nz) / fs;
base = 100 * (1 + 0.05 * sin(2 * pi * t_s / 200 + 2 * pi * rand(NCell, 1)));
F = base .* (1 + 0.6 * calcium_signal) + 4 * randn(NCell, Nz);

% Couleurs : tirées au hasard, puis 9 cellules d'un même clone dans les assemblées 1 et 2
p_couleurs = [.25 .2 .2 .1 .1 .05 .05 .05];
colorcell = randsample(8, NCell, true, p_couleurs);
couleur_clone = [4 5];                     % RG et RB
clonale = [true true false false];
for k = 1:2
    clone = assemblees{k}(randperm(15, 9));
    colorcell(clone) = couleur_clone(k);
end

% ROI suite2p (pixels comptés à partir de 0) ; les non-cellules sont mélangées
nROI = NCell + n_non_cellules;
est_cellule = [ones(NCell, 1); zeros(n_non_cellules, 1)];
est_cellule = est_cellule(randperm(nROI));
iscell = [est_cellule, est_cellule]; %#ok<NASGU>
stat = cell(1, nROI);
[gx, gy] = meshgrid(-5:5);
disque = gx.^2 + gy.^2 <= 25;
angle = linspace(0, 2 * pi, 25);  angle(end) = [];
img = 200 * ones(256);
for r = 1:nROI
    cx = randi([10 246]);  cy = randi([10 246]);
    roi = struct();
    roi.xpix = int32(cx + gx(disque))' - 1;
    roi.ypix = int32(cy + gy(disque))' - 1;
    roi.xext = int32(round(cx + 5 * cos(angle))) - 1;
    roi.yext = int32(round(cy + 5 * sin(angle))) - 1;
    stat{r} = roi;
    img(sub2ind([256 256], cy + gy(disque), cx + gx(disque))) = 1000;
end
save(fullfile(dossier_s2p, 'Fall.mat'), 'stat', 'iscell');
imwrite(uint16(img), fullfile(dossier_s2p, 'calcium.tif'));
save(fullfile(dossier_s2p, 'colorcell_GT.mat'), 'colorcell');

data_recording = struct('identifier', sprintf('%s_260101_plane0', id_animal), 'F', F, ...
    'sampling_rate', fs, 'speed', speed, 'path', dossier_s2p); %#ok<NASGU>
chemin_data = fullfile(racine, ['data_' id_animal '.mat']);
save(chemin_data, 'data_recording', '-v7.3');
verite = struct('assemblees', {assemblees}, 'clonale', clonale);
end
