function bilan = reextraire_fluo_masques(dossier, enregistrer)
%% =========================================================================
%  REEXTRAIRE_FLUO_MASQUES (18 nov. 2026, demandée par l'utilisateur)
%  Animaux 446178 et 446198 : les fluorescences de data_recording.mat
%  viennent des masques suite2p trop gros (halo de pixels faibles). On les
%  recalcule avec les masques corrigés (corriger_masques, color_maxwellV2 :
%  mêmes masques que pour les couleurs), à partir du film recalé contenu dans
%  le NWB de la session (/acquisition/motion_corrected_ci_movie, uint16).
%
%  Le film du NWB est une copie remise à l'échelle du film de suite2p
%  (x_suite2p = a * x_nwb + b, a et b estimés sur la session) et quelques
%  pixels très brillants y sont saturés (65535). Pour ne pas en souffrir, on
%  part de l'ancienne trace (exacte) et on retire seulement la part des
%  pixels enlevés (halo, faibles, rarement saturés) :
%     F_nouv = (F_anc * S_tous - sum_{enlevés} lam * (a * x + b)) / S_gardés
%  (F de suite2p = somme des pixels pondérés par lam / somme des lam,
%  allow_overlap = 1). Le neuropile n'est pas touché (le pipeline utilise F
%  brut).
%
%  enregistrer = true : data_recording.mat d'origine copié en
%  data_recording_avant_masquescorr.mat (une seule fois), puis
%  data_recording.F remplacé ; champ data_recording.F_source ajouté.
%  Si la copie existe déjà, F_anc est relu dans la copie (relance sans risque).
%  =========================================================================
if nargin < 2, enregistrer = false; end
d = regexprep(char(dossier), '[\\/]+$', '');
f_dr = fullfile(d, 'data_recording.mat');
f_bak = fullfile(d, 'data_recording_avant_masquescorr.mat');
if isfile(f_bak), R = load(f_bak); else, R = load(f_dr); end
R = R.data_recording;
F_anc = double(R.F);
Fa = load(fullfile(d, 'Fall.mat'), 'stat', 'iscell', 'ops');
roi = find(Fa.iscell(:, 1) == 1);
if size(F_anc, 1) ~= numel(roi)
    error('%s : data_recording a %d cellules, Fall.mat %d (iscell).', d, size(F_anc, 1), numel(roi));
end
[~, stat_n] = corriger_masques(d, [], false);
Ly = double(Fa.ops.Ly); Lx = double(Fa.ops.Lx);

% Film : d'abord un tif recalé du dossier (> 1 Go, avec le bon nombre
% d'images ; 18 nov. : films complets transférés pour 2 sessions), sinon le
% NWB qui contient le film avec le bon nombre d'images
Nt = size(F_anc, 2);
source = ''; f_film = ''; transposer = true;
% 18 nov. : film refait à partir des tifs bruts (recaler_avec_decalages.py)
% ou data.bin d'origine de suite2p, copié dans le dossier de la session
for nom_bin = {'data_recale.bin', 'data.bin'}
    fb = fullfile(d, nom_bin{1});
    if isempty(source) && isfile(fb)
        Lb = dir(fb);
        if Lb.bytes == 2 * Lx * Ly * Nt, source = 'bin'; f_film = fb; end
    end
end
L = dir(fullfile(d, '*.tif'));
if ~isempty(source), L = L([]); end
L = L([L.bytes] > 1e9);
for k = 1:numel(L)
    f = fullfile(d, L(k).name);
    info_tif = imfinfo(f);
    if numel(info_tif) == Nt && info_tif(1).Height == Ly && info_tif(1).Width == Lx
        source = 'tif'; f_film = f; break;
    end
end
if isempty(source)
    L = dir(fullfile(d, '*.nwb'));
    for k = 1:numel(L)
        f = fullfile(d, L(k).name);
        try
            I = h5info(f, '/acquisition/motion_corrected_ci_movie/data');
            if I.Dataspace.Size(3) == Nt, source = 'nwb'; f_film = f; break; end
        catch
        end
    end
    if isempty(source), error('%s : pas de NWB ni de tif avec le film (%d images).', d, Nt); end
    taille = double(I.Dataspace.Size);             % [x y t] vu de MATLAB
    if taille(1) ~= Lx || taille(2) ~= Ly, error('%s : taille du film différente de ops.', d); end
end
taille = [Lx Ly Nt];
lire = @(t, nt) lire_film(source, f_film, t, nt, Lx, Ly, transposer);

% poids : tous les pixels (pour estimer a et b) et pixels enlevés
n = numel(roi);
[i_tous, j_tous, v_tous, i_rem, j_rem, v_rem] = deal([]);
S_tous = zeros(n, 1); S_gard = zeros(n, 1); n_avant = zeros(n, 1); n_apres = zeros(n, 1);
for m = 1:n
    a0 = Fa.stat{roi(m)}; a1 = stat_n{roi(m)};
    x = double(a0.xpix(:)); y = double(a0.ypix(:)); lam = double(a0.lam(:));
    idx = sub2ind([Lx Ly], x + 1, y + 1);
    idx1 = sub2ind([Lx Ly], double(a1.xpix(:)) + 1, double(a1.ypix(:)) + 1);
    rem = ~ismember(idx, idx1);
    S_tous(m) = sum(lam); S_gard(m) = sum(lam(~rem));
    n_avant(m) = numel(idx); n_apres(m) = sum(~rem);
    i_tous = [i_tous; m * ones(numel(idx), 1)]; j_tous = [j_tous; idx]; v_tous = [v_tous; lam / sum(lam)]; %#ok<AGROW>
    i_rem = [i_rem; m * ones(sum(rem), 1)]; j_rem = [j_rem; idx(rem)]; v_rem = [v_rem; lam(rem)]; %#ok<AGROW>
end
W_tous = sparse(i_tous, j_tous, v_tous, n, Lx * Ly);
W_rem = sparse(i_rem, j_rem, v_rem, n, Lx * Ly);
P_rem = full(sum(W_rem, 2));                   % somme des lam enlevés

% a et b : sur les 1000 premières images, cellules sans pixel saturé
T0 = min(1000, taille(3));
[a, b, err_rel, ab] = estimer_echelle(lire(1, T0), W_tous, F_anc(:, 1:T0));
if strcmp(source, 'tif') && ~(err_rel < 1e-2)
    % tif dans l'autre orientation ?
    transposer = false; lire = @(t, nt) lire_film(source, f_film, t, nt, Lx, Ly, transposer);
    [a, b, err_rel, ab] = estimer_echelle(lire(1, T0), W_tous, F_anc(:, 1:T0));
end
if ~(err_rel < 1e-2)
    error('%s : le film ne correspond pas aux traces de suite2p (erreur relative %.3g).', d, err_rel);
end

% part des pixels enlevés, sur tout le film (par blocs)
Rr = zeros(n, taille(3)); err_blocs = [];
bloc = 1000;
for t = 1:bloc:taille(3)
    nt = min(bloc, taille(3) - t + 1);
    M = lire(t, nt);
    Rr(:, t:t + nt - 1) = a * (W_rem * M) + b * P_rem;
    % contrôle sur tout le film : l'ancienne trace est-elle retrouvée ? (images décalées, manquantes...)
    ecart = abs(F_anc(:, t:t + nt - 1) - (a * (W_tous * M) + b)) ./ mean(F_anc(:, t:t + nt - 1), 2);
    err_blocs(end + 1) = median(max(ecart, [], 2)); %#ok<AGROW>
end
if max(err_blocs) > 1e-2
    error('%s : le film ne retrouve pas les anciennes traces partout (écart max par bloc %.3g).', d, max(err_blocs));
end
F_nouv = (F_anc .* S_tous - Rr) ./ S_gard;

c = arrayfun(@(m) corr(F_anc(m, :)', F_nouv(m, :)'), 1:n);
bilan = struct('dossier', d, 'film', f_film, 'source', source, 'n_cellules', n, 'a', a, 'b', b, ...
    'n_cellules_fit', sum(~isnan(ab(:, 1))), 'erreur_rel_fit', err_rel, 'erreur_max_blocs', max(err_blocs), ...
    'npix_median_avant', median(n_avant), 'npix_median_apres', median(n_apres), ...
    'corr_median_anc_nouv', median(c), 'corr_min_anc_nouv', min(c), ...
    'rapport_moyenne_nouv_anc', median(mean(F_nouv, 2) ./ mean(F_anc, 2)));
if enregistrer
    if ~isfile(f_bak), copyfile(f_dr, f_bak); end
    data_recording = R;
    data_recording.F = cast(F_nouv, 'like', R.F);
    data_recording.F_source = sprintf(['F recalculé le %s avec les masques corrigés (corriger_masques, ' ...
        'pixels de plus fort lam) à partir du film recalé (%s) ; ancien F dans data_recording_avant_masquescorr.mat'], ...
        datestr(now, 'yyyy-mm-dd'), source); %#ok<DATST,TNOW1>
    save(f_dr, 'data_recording', '-v7.3');
end
end


function M = lire_film(source, f, t, nt, Lx, Ly, transposer)
% images t à t+nt-1, en colonnes de Lx*Ly pixels dans l'ordre (x, y)
if strcmp(source, 'nwb')
    M = double(reshape(h5read(f, '/acquisition/motion_corrected_ci_movie/data', [1 1 t], [Lx Ly nt]), Lx * Ly, nt));
    return
end
if strcmp(source, 'bin')
    % int16, images Ly x Lx écrites ligne par ligne (Python) : x varie le plus vite
    fid = fopen(f, 'r'); fermer = onCleanup(@() fclose(fid));
    fseek(fid, 2 * Lx * Ly * (t - 1), 'bof');
    M = double(reshape(fread(fid, Lx * Ly * nt, 'int16=>int16'), Lx * Ly, nt));
    return
end
M = zeros(Lx * Ly, nt);
T = Tiff(f, 'r'); fermer = onCleanup(@() close(T));
setDirectory(T, t);
for k = 1:nt
    im = double(read(T));                      % Ly x Lx (lignes = y)
    if transposer, im = im'; end
    M(:, k) = im(:);
    if k < nt, nextDirectory(T); end
end
end


function [a, b, err_rel, ab] = estimer_echelle(M, W_tous, F0)
% x_suite2p = a * x_film + b, estimé sur les cellules sans pixel saturé
G = W_tous * M;
hi = max(M(:)); lo = min(M(:));
sat = (W_tous > 0) * double(M >= hi | M <= lo) > 0;
n = size(W_tous, 1); ab = nan(n, 2);
for m = 1:n
    if any(sat(m, :)), continue; end
    ab(m, :) = polyfit(G(m, :), F0(m, :), 1);
end
a = median(ab(:, 1), 'omitnan'); b = median(ab(:, 2), 'omitnan');
res = F0 - (a * G + b);
ok = ~any(sat, 2) & ~isnan(ab(:, 1));
err_rel = median(max(abs(res(ok, :)), [], 2) ./ mean(F0(ok, :), 2));
end
