function T = paires_brutes_session(session, config)
%PAIRES_BRUTES_SESSION  Similarité BRUTE des paires de même couleur et de
%   couleurs différentes, pour UNE session (étape 3b, 5 nov.).
%   « Brute » : la valeur de la mesure telle quelle (STTC, Pearson du dF/F,
%   ...), SANS résidus : ni l'effet de la distance ni le niveau de chaque
%   cellule ne sont retirés. Aucune permutation : seulement les moyennes
%   observées, donc rapide.
%   Mêmes paires que l'étape 3b : cellules colorées (1 à 6), ROI qui se
%   chevauchent exclues, mêmes distances (disques ou tranches), mêmes
%   mesures ; une distance avec moins de config.paires.min_paires paires de
%   même couleur ou de couleurs différentes : NaN.
%   Écrit <session>\paires\paires_brutes.mat (T, bornes, distances) et
%   renvoie T : une ligne par mesure x distance (toutes couleurs).

opts = config.paires;
dossier = fullfile(session.dossier, 'paires');
if ~isfolder(dossier), mkdir(dossier); end

% Activité au repos (aucune couleur)
pre = load(fullfile(session.dossier_source, 'pretraitement.mat'), 'Tr1b', 'Raster', 'WinRest', 'sampling_rate', 'Race');
NCell = size(pre.Raster, 1);
S = similarite_paires(pre.Tr1b, pre.Raster, pre.WinRest, pre.sampling_rate, opts, pre.Race);
clear pre

% Géométrie : distances en µm, ROI qui se chevauchent exclues
path = avec_separateur(session.dossier_suite2p);
[~, ~, gcamp_centroid, ~, ~, Coomasqx, Coomasqy] = charger_geometrie(path);
um = session.um_par_px;
if isnan(um), um = config.um_par_px_defaut; end
xy = gcamp_centroid * um;
dist_um = sqrt((xy(:, 1) - xy(:, 1)') .^ 2 + (xy(:, 2) - xy(:, 2)') .^ 2);
V0 = ~masques_chevauchent(Coomasqx, Coomasqy) & isfinite(dist_um);
V0(logical(eye(NCell))) = false;
haut = triu(true(NCell), 1);                 % chaque paire une seule fois

% Couleurs
colorcell = charger_couleurs(path, session.fichier_colorcell);
colorcell = couleurs_sans_doublons(colorcell, session);    % doublons de l'étape 1 : non colorés
col = colorcell(:);
ok = col >= 1 & col <= 6;
meme = ok & ok' & (col == col');
differentes = ok & ok' & (col ~= col');

[bornes, cumulatif] = distances_paires(opts);
mesures = liste_mesures_paires(opts);
lignes = cell(0, 11);
for im = 1:numel(mesures)
    C = S.(mesures{im});
    V = V0 & isfinite(C) & haut;
    for b = 1:numel(bornes) - 1
        if cumulatif
            m = V & dist_um <= bornes(b + 1);
            dmin = 0;
        else
            m = V & dist_um >= bornes(b) & dist_um < bornes(b + 1);
            dmin = bornes(b);
        end
        cm = C(m & meme);
        cd = C(m & differentes);
        b_meme = NaN; b_diff = NaN;
        if numel(cm) >= opts.min_paires && numel(cd) >= opts.min_paires
            b_meme = mean(cm);
            b_diff = mean(cd);
        end
        lignes(end+1, :) = {session.id_animal, session.id_champ, session.id_session, session.id_condition, ...
            mesures{im}, dmin, bornes(b + 1), b_meme, b_diff, numel(cm), numel(cd)}; %#ok<AGROW>
    end
end
T = cell2table(lignes, 'VariableNames', {'animal', 'champ', 'session', 'condition', 'mesure', ...
    'dist_min_um', 'dist_max_um', 'brut_meme', 'brut_diff', 'n_paires_meme', 'n_paires_diff'});
distances = 'tranches';
if isfield(opts, 'distances'), distances = opts.distances; end
save(fullfile(dossier, 'paires_brutes.mat'), 'T', 'bornes', 'distances');
end
