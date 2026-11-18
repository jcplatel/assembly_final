function T = diagnostic_doublons(config)
%% =========================================================================
%  DIAGNOSTIC_DOUBLONS : ROI qui sont probablement la même cellule comptée
%  deux fois, et ROI trop petites ou trop grandes, dans toutes les sessions
%  déjà prétraitées (7 nov.)
%  -------------------------------------------------------------------------
%  NE RECALCULE ET NE MODIFIE RIEN : lit pretraitement.mat (étape 1), les
%  masques suite2p (Fall.mat) et, pour information, les couleurs et les
%  assemblées au K retenu. Même critère que config.pre.retirer_doublons
%  (detecter_doublons : traces au repos corrélées à plus de doublons_r_min
%  et masques qui se touchent ou centres à moins d'un diamètre).
%  Tailles : même règle que config.pre.retirer_tailles (tailles_roi : moins
%  de roi_min_px pixels, ou plus de roi_max_facteur x la surface médiane).
%  Sert à voir l'ampleur du problème avant de refaire l'étape 1.
%
%  Lancer :   diagnostic_doublons
%             config = config_brainbow(); config.seulement = {'444175'}; diagnostic_doublons(config)
%  Écrit <PathSave>\doublons_paires.csv (une ligne par paire),
%  roi_hors_taille.csv (une ligne par ROI trop petite ou trop grande) et
%  doublons_sessions.csv (une ligne par session, avec les tailles des ROI).
%  =========================================================================

if nargin < 1, config = config_brainbow(); end
addpath(genpath(fileparts(mfilename('fullpath'))));
config = dossiers_analyse(config);
sessions = lire_table_sessions(config);
a_traiter = 1:height(sessions);
if isfield(config, 'seulement') && ~isempty(config.seulement)
    seulement = cellstr(config.seulement);
    a_traiter = find(ismember(sessions.nom, seulement) | ismember(sessions.id_animal, seulement))';
end
r_min = 0.7;
if isfield(config.pre, 'doublons_r_min') && ~isempty(config.pre.doublons_r_min)
    r_min = config.pre.doublons_r_min;
end
min_px = 4; max_facteur = 3;
if isfield(config.pre, 'roi_min_px'), min_px = config.pre.roi_min_px; end
if isfield(config.pre, 'roi_max_facteur'), max_facteur = config.pre.roi_max_facteur; end
noms_c = {'R', 'G', 'B', 'RG', 'RB', 'GB', 'all', 'none'};
lignes = cell(0, 12);
lignes_t = cell(0, 7);
resume = cell(0, 11);
for k = a_traiter
    session = table2struct(sessions(k, :));
    f = fullfile(session.dossier_source, 'pretraitement.mat');
    if ~isfile(f), continue; end
    try
        P = load(f, 'Tr1b', 'WinRest');
        NCell = size(P.Tr1b, 1);
        path = avec_separateur(session.dossier_suite2p);
        [~, ~, ~, ~, ~, Cx, Cy] = charger_geometrie(path);
        Ta = tailles_roi(Cx, min_px, max_facteur);
        Tr_d = P.Tr1b;
        Tr_d(Ta.retirees, :) = 0;
        D = detecter_doublons(Tr_d, P.WinRest, Cx, Cy, r_min);
        col = nan(NCell, 1);
        try
            c = charger_couleurs(path, session.fichier_colorcell);
            if numel(c) == NCell, col = c(:); end
        catch
        end
        asm = assemblees_retenues(session, NCell);
        um = session.um_par_px;
        if isnan(um), um = config.um_par_px_defaut; end
        for m = 1:height(D.paires)
            i = D.paires.roi_1(m);
            j = D.paires.roi_2(m);
            ensemble = any(cellfun(@(a) ismember(i, a) && ismember(j, a), asm));
            une_seule = any(cellfun(@(a) xor(ismember(i, a), ismember(j, a)), asm));
            lignes(end+1, :) = {session.id_animal, session.nom, i, j, D.paires.r(m), ...
                D.paires.distance_px(m) * um, D.paires.se_touchent(m), D.paires.retiree(m), ...
                nom_couleur(col(i), noms_c), nom_couleur(col(j), noms_c), ensemble, une_seule}; %#ok<AGROW>
        end
        for m = Ta.retirees
            if ismember(m, Ta.petites), quoi = 'trop petite'; else, quoi = 'trop grande'; end
            lignes_t(end+1, :) = {session.id_animal, session.nom, m, Ta.npix(m), Ta.npix(m) / Ta.mediane, ...
                quoi, nom_couleur(col(m), noms_c)}; %#ok<AGROW>
        end
        resume(end+1, :) = {session.id_animal, session.nom, NCell, height(D.paires), numel(D.retirees), ...
            min(Ta.npix), Ta.mediane, max(Ta.npix), numel(Ta.petites), numel(Ta.grandes), ...
            Ta.mediane * um ^ 2}; %#ok<AGROW>
        fprintf('[%s] %d cellules : %d doublon(s) %s ; ROI de %d à %d pixels (médiane %g) : %d trop petite(s), %d trop grande(s)\n', ...
            session.nom, NCell, numel(D.retirees), mat2str(D.retirees), min(Ta.npix), max(Ta.npix), Ta.mediane, ...
            numel(Ta.petites), numel(Ta.grandes));
    catch err
        fprintf(2, '[%s] : ÉCHEC : %s\n', session.nom, err.message);
    end
end
T = cell2table(lignes, 'VariableNames', {'animal', 'session', 'roi_1', 'roi_2', 'r', 'distance_um', ...
    'se_touchent', 'roi_retiree', 'couleur_1', 'couleur_2', 'meme_assemblee', 'une_seule_dans_assemblee'});
R = cell2table(resume, 'VariableNames', {'animal', 'session', 'n_cellules', 'n_paires', 'n_retirees', ...
    'pixels_min', 'pixels_mediane', 'pixels_max', 'n_trop_petites', 'n_trop_grandes', 'surface_mediane_um2'});
Tt = cell2table(lignes_t, 'VariableNames', {'animal', 'session', 'roi', 'pixels', 'fois_la_mediane', ...
    'raison', 'couleur'});
writetable(T, fullfile(config.PathSave, 'doublons_paires.csv'));
writetable(R, fullfile(config.PathSave, 'doublons_sessions.csv'));
writetable(Tt, fullfile(config.PathSave, 'roi_hors_taille.csv'));
if isempty(R), return; end
fprintf(['\n%d sessions : %d doublons à retirer sur %d ROI (%.1f %%), de 0 à %d par session ; ' ...
    '%d paires de même couleur sur %d ; %d paires toutes les deux dans une même assemblée.\n'], ...
    height(R), sum(R.n_retirees), sum(R.n_cellules), 100 * sum(R.n_retirees) / sum(R.n_cellules), ...
    max(R.n_retirees), sum(strcmp(T.couleur_1, T.couleur_2)), height(T), sum(T.meme_assemblee));
fprintf(['Tailles : ROI de %d à %d pixels ; %d trop petites (< %d px), %d trop grandes ' ...
    '(> %g x la médiane de leur session).\n'], min(R.pixels_min), max(R.pixels_max), ...
    sum(R.n_trop_petites), min_px, sum(R.n_trop_grandes), max_facteur);
end


function t = nom_couleur(c, noms_c)
if isnan(c) || c < 1 || c > numel(noms_c), t = '?'; else, t = noms_c{c}; end
end
