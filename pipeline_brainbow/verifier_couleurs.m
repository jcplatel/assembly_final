function [Ts, Tc] = verifier_couleurs(config)
%% =========================================================================
%  VERIFIER_COULEURS : contrôle des fichiers de couleurs (colorcell) de
%  toutes les sessions (10 nov., demande de l'utilisateur : « c'est possible
%  que des erreurs se soient glissées »)
%  -------------------------------------------------------------------------
%  NE MODIFIE RIEN. Lit, pour chaque session de la table : le fichier de
%  couleurs que le pipeline utilise (charger_couleurs : colonne
%  fichier_colorcell de la table, sinon colorcell_GT, _Maxwellnew,
%  _Maxwell), Fall.mat (stat, iscell) et, s'il existe, pretraitement.mat
%  (nombre de lignes du Raster). Les couleurs sont lues BRUTES (doublons et
%  ROI hors taille non retirés).
%
%  Pour chaque session :
%   1. TAILLE : nombre de valeurs de colorcell = nombre de cellules
%      iscell = 1 de Fall.mat = nombre de lignes du Raster ? Si colorcell a
%      autant de valeurs que TOUTES les ROI de Fall.mat, il a été calculé
%      avant la curation (iscell).
%   2. VALEURS : uniquement 1 à 8 (1 R, 2 G, 3 B, 4 RG, 5 RB, 6 GB, 7 blanc,
%      8 sans couleur) ? NaN, 0, valeurs non entières ou > 8 sont signalés.
%   3. AUTRES FICHIERS : si plusieurs fichiers de couleurs existent dans le
%      dossier (GT, Maxwellnew, Maxwell), accord avec celui utilisé.
%   4. ALIGNEMENT couleurs / ROI : les clones rendent les voisines plus
%      souvent de même couleur. Part des 3 plus proches voisines colorées
%      qui ont la même couleur, divisée par la même part après mélange des
%      couleurs (500 mélanges) : indice 1 = aucun regroupement. Le même
%      indice avec colorcell décalé d'une position (+1 et -1) doit retomber
%      à 1. Un indice proche de 1, ou un indice décalé plus fort que le
%      vrai, fait soupçonner des couleurs mal alignées sur les ROI (ordre
%      différent). Ce n'est pas une preuve : un champ sans clone visible
%      donne aussi 1.
%   6. IMAGE SUR LES MASQUES (11 nov., après 443669_220615_plane0) : dans
%      compositecalcium.tif (rouge, vert, bleu, calcium), cherche le
%      déplacement de l'image qui met le plus de lumière dans les masques
%      suite2p (alignement_composite) ; bien posée = meilleur décalage
%      proche de 0. Voit un recalage raté, que le contrôle 4 ne voit pas
%      (une image décalée d'un bloc garde des couleurs regroupées).
%   5. COPIES : deux sessions avec exactement les mêmes couleurs, ou le
%      même fichier, sont signalées.
%  Composition (% de chaque couleur parmi les cellules colorées) :
%   - dans un même animal : écart de chaque session à la moyenne des autres
%     sessions de l'animal (distance de variation totale, en points de %) ;
%     sessions au-delà de médiane + 3 écarts absolus médians signalées ;
%   - entre animaux : composition moyenne de chaque animal.
%  Les champs d'un même animal ne sont pas les mêmes cellules (autre jour,
%  autre profondeur, autres clones) : une composition différente n'est pas
%  forcément une erreur, c'est une session à regarder.
%
%  Lancer :   verifier_couleurs
%             config = config_brainbow(); config.seulement = {'444175'}; verifier_couleurs(config)
%  Écrit dans <dossier de l'étape 1>\verification_couleurs\ :
%  verification_couleurs_sessions.csv (une ligne par session),
%  composition_animaux.csv, resume_verification_couleurs.txt (aussi affiché),
%  verification_couleurs.png.
%  =========================================================================

if nargin < 1, config = config_brainbow(); end
addpath(genpath(fileparts(mfilename('fullpath'))));
etat_rng = rng;
remettre_rng = onCleanup(@() rng(etat_rng)); %#ok<NASGU>
config = dossiers_analyse(config);
sessions = lire_table_sessions(config);
a_traiter = 1:height(sessions);
if isfield(config, 'seulement') && ~isempty(config.seulement)
    seulement = cellstr(config.seulement);
    a_traiter = find(ismember(sessions.nom, seulement) | ismember(sessions.id_animal, seulement))';
end
noms_c = {'R', 'G', 'B', 'RG', 'RB', 'GB'};
dossier = fullfile(config.PathEtape1, 'verification_couleurs');
if ~isfolder(dossier), mkdir(dossier); end
fprintf('=== VÉRIFICATION DES COULEURS : %d session(s) ===\n', numel(a_traiter));

n = numel(a_traiter);
animal = cell(n, 1); nom = cell(n, 1); fichier = cell(n, 1); source = cell(n, 1);
n_colorcell = nan(n, 1); n_iscell = nan(n, 1); n_roi = nan(n, 1); n_raster = nan(n, 1);
n_invalides = nan(n, 1); valeurs_invalides = repmat({''}, n, 1);
comp = nan(n, 6); n_col = nan(n, 1); pct_blanc = nan(n, 1); pct_sans = nan(n, 1);
autres = repmat({''}, n, 1);
indice = nan(n, 1); p_indice = nan(n, 1); indice_plus1 = nan(n, 1); indice_moins1 = nan(n, 1);
probleme = repmat({''}, n, 1);
image_composite = repmat({''}, n, 1); dec_bb = nan(n, 2); c_bb0 = nan(n, 1); c_bbmax = nan(n, 1);
dec_ca = nan(n, 2); c_ca0 = nan(n, 1); bord_vide = nan(n, 1);
vecteurs = cell(n, 1);
for s = 1:n
    session = table2struct(sessions(a_traiter(s), :));
    animal{s} = char(session.id_animal);
    nom{s} = char(session.nom);
    pb = {};
    try
        path = avec_separateur(session.dossier_suite2p);
        [cc, fichier{s}, source{s}] = charger_couleurs(path, session.fichier_colorcell);
        vecteurs{s} = cc;
        n_colorcell(s) = numel(cc);

        % 1. tailles
        if isfile([path 'Fall.mat'])
            F = load([path 'Fall.mat'], 'iscell');
            n_roi(s) = size(F.iscell, 1);
            n_iscell(s) = sum(F.iscell(:, 1) == 1);
        else
            pb{end+1} = 'Fall.mat introuvable'; %#ok<AGROW>
        end
        dsrc = session.dossier;
        if isfield(session, 'dossier_source') && ~isempty(session.dossier_source), dsrc = session.dossier_source; end
        fp = fullfile(dsrc, 'pretraitement.mat');
        if isfile(fp)
            info = whos('-file', fp, 'Raster');
            if ~isempty(info), n_raster(s) = info.size(1); end
        end
        if ~isnan(n_iscell(s)) && n_colorcell(s) ~= n_iscell(s)
            if n_colorcell(s) == n_roi(s)
                pb{end+1} = sprintf(['colorcell a %d valeurs = TOUTES les ROI de Fall.mat (%d cellules ' ...
                    'iscell) : calculé avant la curation'], n_colorcell(s), n_iscell(s)); %#ok<AGROW>
            else
                pb{end+1} = sprintf('colorcell a %d valeurs, Fall.mat %d cellules iscell', ...
                    n_colorcell(s), n_iscell(s)); %#ok<AGROW>
            end
        end
        if ~isnan(n_raster(s)) && ~isnan(n_iscell(s)) && n_raster(s) ~= n_iscell(s)
            pb{end+1} = sprintf('Raster de l''étape 1 : %d cellules, Fall.mat : %d (curation changée depuis ?)', ...
                n_raster(s), n_iscell(s)); %#ok<AGROW>
        end

        % 2. valeurs
        inv = ~(isfinite(cc) & cc == round(cc) & cc >= 1 & cc <= 8);
        n_invalides(s) = sum(inv);
        if any(inv)
            u = unique(cc(inv));
            valeurs_invalides{s} = strjoin(arrayfun(@(v) sprintf('%g', v), u(1:min(end, 5))', ...
                'UniformOutput', false), ' ');
            pb{end+1} = sprintf('%d valeur(s) hors 1-8 (%s)', sum(inv), valeurs_invalides{s}); %#ok<AGROW>
        end
        ok = ~inv;
        colo = ok & cc >= 1 & cc <= 6;
        n_col(s) = sum(colo);
        comp(s, :) = 100 * accumarray(cc(colo), 1, [6 1])' / max(n_col(s), 1);
        pct_blanc(s) = 100 * mean(cc == 7);
        pct_sans(s) = 100 * mean(cc == 8);

        % 3. autres fichiers de couleurs dans le dossier
        cand = {'colorcell_GT.mat', 'colorcell_Maxwellnew.mat', 'colorcell_Maxwell.mat'};
        txt = {};
        for f = 1:numel(cand)
            fa = [path cand{f}];
            if ~isfile(fa) || strcmpi(fa, fichier{s}), continue; end
            A = load(fa, 'colorcell');
            if ~isfield(A, 'colorcell'), continue; end
            a = double(A.colorcell(:));
            if numel(a) ~= numel(cc)
                txt{end+1} = sprintf('%s : %d valeurs', cand{f}, numel(a)); %#ok<AGROW>
            else
                txt{end+1} = sprintf('%s : %.0f %% identique', cand{f}, 100 * mean(a == cc)); %#ok<AGROW>
            end
        end
        autres{s} = strjoin(txt, ' ; ');

        % 4. alignement : regroupement des couleurs entre voisines
        if ~isnan(n_iscell(s)) && n_colorcell(s) == n_iscell(s) && n_col(s) >= 20
            [~, ~, xy] = charger_geometrie(path);
            fixer_flux_aleatoire(session.graine, 40, 1);
            [indice(s), p_indice(s), indice_plus1(s), indice_moins1(s)] = regroupement_voisines(cc, xy, 3, 500);
            if indice(s) < 1.05 && p_indice(s) > 0.05
                pb{end+1} = sprintf('aucun regroupement des couleurs entre voisines (indice %.2f)', indice(s)); %#ok<AGROW>
            end
            if max(indice_plus1(s), indice_moins1(s)) > indice(s)
                pb{end+1} = sprintf(['couleurs décalées d''une position plus regroupées que les vraies ' ...
                    '(%.2f / %.2f contre %.2f) : alignement à vérifier'], indice_plus1(s), indice_moins1(s), ...
                    indice(s)); %#ok<AGROW>
            end
        end

        % 6. image Brainbow posée sur les masques ? (11 nov., compositecalcium.tif)
        fc = trouver_composite(path);
        if isempty(fc)
            image_composite{s} = '(aucune)';
        else
            [~, nf, ef] = fileparts(fc);
            image_composite{s} = [nf ef];
            Al = alignement_composite(fc, path, 80);
            dec_bb(s, :) = Al.decalage_bb; c_bb0(s) = Al.contraste_bb_0; c_bbmax(s) = Al.contraste_bb_max;
            dec_ca(s, :) = Al.decalage_ca; c_ca0(s) = Al.contraste_ca_0; bord_vide(s) = Al.bord_vide_px;
            if any(abs(Al.decalage_bb) > 4) && Al.contraste_bb_max > 1.05 * Al.contraste_bb_0
                pb{end+1} = sprintf(['IMAGE BRAINBOW DÉCALÉE : les masques tombent mieux sur les cellules avec ' ...
                    'l''image déplacée de %d px en x et %d px en y (contraste %.2f au lieu de %.2f) : couleurs ' ...
                    'probablement fausses'], Al.decalage_bb(1), Al.decalage_bb(2), Al.contraste_bb_max, ...
                    Al.contraste_bb_0); %#ok<AGROW>
            elseif Al.contraste_bb_0 < 1.05
                pb{end+1} = sprintf(['les masques ne tombent pas sur des cellules brillantes de l''image Brainbow ' ...
                    '(contraste %.2f)'], Al.contraste_bb_0); %#ok<AGROW>
            end
            if Al.bord_vide_px > 10
                pb{end+1} = sprintf('bord vide de %d px dans un canal Brainbow (recalage par grand décalage)', ...
                    Al.bord_vide_px); %#ok<AGROW>
            end
            if any(abs(Al.decalage_ca) > 2)
                pb{end+1} = sprintf('canal calcium du composite décalé de [%d %d] px par rapport aux masques', ...
                    Al.decalage_ca(1), Al.decalage_ca(2)); %#ok<AGROW>
            end
        end
    catch err
        pb{end+1} = ['ÉCHEC : ' err.message]; %#ok<AGROW>
    end
    probleme{s} = strjoin(pb, ' | ');
    if isempty(pb)
        fprintf('[%s] ok (%d cellules, %d colorées)\n', nom{s}, n_colorcell(s), n_col(s));
    else
        fprintf(2, '[%s] %s\n', nom{s}, probleme{s});
    end
end

% 5. copies entre sessions
for s = 1:n
    for t = s+1:n
        if isempty(vecteurs{s}) || isempty(vecteurs{t}), continue; end
        msg = '';
        if strcmpi(fichier{s}, fichier{t})
            msg = 'même fichier de couleurs que ';
        elseif numel(vecteurs{s}) == numel(vecteurs{t}) && isequaln(vecteurs{s}, vecteurs{t})
            msg = 'couleurs identiques à ';
        end
        if ~isempty(msg)
            probleme{s} = strjoin([probleme(s), {[msg nom{t}]}], ' | ');
            probleme{t} = strjoin([probleme(t), {[msg nom{s}]}], ' | ');
        end
    end
end
probleme = regexprep(probleme, '^ \| ', '');

% Composition : écart à la moyenne des autres sessions de l'animal
ecart_animal = nan(n, 1);
for s = 1:n
    autres_s = strcmp(animal, animal{s}) & (1:n)' ~= s & isfinite(comp(:, 1));
    if any(autres_s) && isfinite(comp(s, 1))
        ecart_animal(s) = 0.5 * sum(abs(comp(s, :) - mean(comp(autres_s, :), 1)));
    end
end
e = ecart_animal(isfinite(ecart_animal));
seuil_ecart = median(e) + 3 * 1.4826 * median(abs(e - median(e)));
atypique = ecart_animal > seuil_ecart;

Ts = table(animal, nom, source, fichier, n_colorcell, n_iscell, n_roi, n_raster, n_invalides, n_col, ...
    pct_blanc, pct_sans, comp(:, 1), comp(:, 2), comp(:, 3), comp(:, 4), comp(:, 5), comp(:, 6), ...
    ecart_animal, atypique, indice, p_indice, indice_plus1, indice_moins1, image_composite, dec_bb(:, 1), ...
    dec_bb(:, 2), c_bb0, c_bbmax, dec_ca(:, 1), dec_ca(:, 2), c_ca0, bord_vide, autres, probleme, ...
    'VariableNames', {'animal', 'session', 'source', 'fichier', 'n_colorcell', 'n_iscell', 'n_roi_fall', ...
    'n_raster', 'n_invalides', 'n_colorees', 'pct_blanc', 'pct_sans_couleur', 'pct_R', 'pct_G', 'pct_B', ...
    'pct_RG', 'pct_RB', 'pct_GB', 'ecart_animal_points', 'composition_atypique', 'indice_voisines', ...
    'p_voisines', 'indice_decale_plus1', 'indice_decale_moins1', 'image_composite', 'decalage_image_x_px', ...
    'decalage_image_y_px', 'contraste_masques', 'contraste_masques_max', 'decalage_calcium_x_px', ...
    'decalage_calcium_y_px', 'contraste_calcium', 'bord_vide_px', 'autres_fichiers', 'probleme'});
writetable(Ts, fullfile(dossier, 'verification_couleurs_sessions.csv'));

% Animaux
[G, anim] = findgroups(animal);
Tc = table(anim, splitapply(@numel, n_col, G), splitapply(@(x) sum(x, 'omitnan'), n_col, G), ...
    'VariableNames', {'animal', 'n_sessions', 'n_colorees'});
for c = 1:6
    Tc.(['pct_' noms_c{c}]) = splitapply(@(x) mean(x, 'omitnan'), comp(:, c), G);
end
Tc.pct_sans_couleur = splitapply(@(x) mean(x, 'omitnan'), pct_sans, G);
writetable(Tc, fullfile(dossier, 'composition_animaux.csv'));

% Résumé
L = {};
L{end+1} = sprintf('VÉRIFICATION DES COULEURS (verifier_couleurs, %d sessions, %d animaux) :', n, numel(anim));
taille_ok = n_colorcell == n_iscell;
L{end+1} = sprintf(['   Taille : %d sessions où colorcell a le même nombre de valeurs que les cellules iscell ' ...
    'de Fall.mat, %d où il diffère, %d non vérifiables.'], sum(taille_ok), ...
    sum(~taille_ok & isfinite(n_iscell) & isfinite(n_colorcell)), sum(~isfinite(n_iscell) | ~isfinite(n_colorcell)));
L{end+1} = sprintf('   Valeurs hors 1-8 : %d session(s).', sum(n_invalides > 0));
dd = max(abs(dec_bb), [], 2);
L{end+1} = sprintf(['   Image Brainbow sur les masques (compositecalcium.tif, %d session(s) avec image) : %d bien posée(s) ' ...
    '(meilleur décalage <= 4 px), %d décalée(s) ; contraste des masques médian %.2f.'], sum(isfinite(dd)), ...
    sum(dd <= 4), sum(dd > 4), median(c_bb0, 'omitnan'));
L{end+1} = sprintf(['   Regroupement des couleurs entre voisines (indice, 1 = hasard) : médiane %.2f [%.2f-%.2f] ; ' ...
    '%d session(s) sans regroupement net ; avec les couleurs décalées d''une position : médiane %.2f.'], ...
    median(indice, 'omitnan'), min(indice), max(indice), sum(indice < 1.05 & p_indice > 0.05), ...
    median([indice_plus1; indice_moins1], 'omitnan'));
L{end+1} = sprintf(['   Composition : écart de chaque session aux autres sessions de son animal, médiane %.0f ' ...
    'points ; %d session(s) atypiques (> %.0f points).'], median(e), sum(atypique), seuil_ecart);
L{end+1} = '   Composition moyenne par animal (% des cellules colorées ; sans couleur en % de toutes les cellules) :';
for a = 1:height(Tc)
    L{end+1} = sprintf('      %-8s %2d sessions : R %4.1f  G %4.1f  B %4.1f  RG %4.1f  RB %4.1f  GB %4.1f | sans %4.1f', ...
        Tc.animal{a}, Tc.n_sessions(a), Tc.pct_R(a), Tc.pct_G(a), Tc.pct_B(a), Tc.pct_RG(a), Tc.pct_RB(a), ...
        Tc.pct_GB(a), Tc.pct_sans_couleur(a)); %#ok<AGROW>
end
L{end+1} = '   Sessions à regarder :';
a_voir = find(~cellfun(@isempty, probleme) | atypique);
if isempty(a_voir), L{end+1} = '      aucune.'; end
for s = a_voir'
    t = probleme{s};
    if atypique(s)
        t = strjoin([{t}, {sprintf('composition atypique pour %s (écart %.0f points)', animal{s}, ecart_animal(s))}], ' | ');
        t = regexprep(t, '^ \| ', '');
    end
    L{end+1} = sprintf('      %s : %s', nom{s}, t); %#ok<AGROW>
end
fid = fopen(fullfile(dossier, 'resume_verification_couleurs.txt'), 'w', 'n', 'UTF-8');
if fid >= 0
    fprintf(fid, '%s\n', L{:});
    fclose(fid);
end
fprintf('\n%s\n', L{:});
try
    figure_verification(Ts, noms_c, dossier);
catch err
    warning('Figure non faite : %s', err.message);
end
end


function [ind, p, ind_p1, ind_m1] = regroupement_voisines(cc, xy, k, nPerm)
% Part des k plus proches voisines colorées de même couleur, divisée par la
% moyenne après mélange des couleurs (1 = hasard) ; aussi pour colorcell
% décalé de +1 et -1 (doit retomber à 1 si l'alignement est bon)
colo = find(cc >= 1 & cc <= 6);
xc = xy(colo, :);
D = sqrt((xc(:, 1) - xc(:, 1)').^2 + (xc(:, 2) - xc(:, 2)').^2);
D(1:size(D, 1) + 1:end) = Inf;
[~, ordre] = sort(D, 2);
V = ordre(:, 1:k);
part = @(c) mean(mean(c(V) == c(:, ones(1, k)), 2));
c0 = cc(colo);
obs = part(c0);
nul = zeros(nPerm, 1);
for r = 1:nPerm
    nul(r) = part(c0(randperm(numel(c0))));
end
ind = obs / mean(nul);
p = (1 + sum(nul >= obs)) / (nPerm + 1);
ind_p1 = decale(cc, colo, 1, part) / mean(nul);
ind_m1 = decale(cc, colo, -1, part) / mean(nul);
end


function v = decale(cc, colo, pas, part)
% couleurs décalées d'une position ; les cellules qui deviennent non
% colorées gardent une couleur tirée au hasard parmi les colorées (même urne)
c2 = circshift(cc, pas);
c2 = c2(colo);
mauvais = ~(c2 >= 1 & c2 <= 6);
c0 = cc(colo);
c2(mauvais) = c0(randi(numel(c0), sum(mauvais), 1));
v = part(c2);
end


function figure_verification(Ts, noms_c, dossier)
% A. composition de chaque session, rangée par animal ; B. indice de
% regroupement entre voisines, vrai et décalé
rgb = [1 0.15 0.15; 0.15 0.8 0.15; 0.25 0.5 1; 0.95 0.85 0.1; 1 0.2 1; 0.1 0.85 0.85];
[~, ordre] = sortrows(Ts(:, {'animal', 'session'}));
T = Ts(ordre, :);
n = height(T);
fig = figure('Color', 'w', 'Position', [40 60 1500 820]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Colour files check (%d sessions)', n), 'FontWeight', 'bold', 'FontSize', 12);
ax = nexttile(tl); hold(ax, 'on');
C = [T.pct_R, T.pct_G, T.pct_B, T.pct_RG, T.pct_RB, T.pct_GB];
C(~isfinite(C)) = 0;
b = bar(ax, 1:n, C, 0.85, 'stacked', 'EdgeColor', 'none');
for c = 1:6, b(c).FaceColor = rgb(c, :); end
for s = find(T.composition_atypique)'
    text(ax, s, 102, '*', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
end
for s = find(~cellfun(@isempty, T.probleme))'
    text(ax, s, 108, '!', 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', [0.8 0 0], 'FontWeight', 'bold');
end
separateurs(ax, T.animal, [0 112]);
ylim(ax, [0 112]); xlim(ax, [0.4, n + 0.6]);
set(ax, 'XTick', 1:n, 'XTickLabel', T.session, 'XTickLabelRotation', 60, 'FontSize', 7, 'TickDir', 'out', ...
    'TickLabelInterpreter', 'none');
ylabel(ax, '% of coloured cells');
legend(ax, b, noms_c, 'Location', 'eastoutside');
title(ax, 'A. Composition of each session, grouped by animal (* atypical for its animal, ! problem: see summary)', ...
    'FontSize', 10);
box(ax, 'off');
ax = nexttile(tl); hold(ax, 'on');
plot(ax, 1:n, T.indice_decale_plus1, 'x', 'Color', [0.6 0.6 0.6], 'MarkerSize', 7);
plot(ax, 1:n, T.indice_decale_moins1, '+', 'Color', [0.6 0.6 0.6], 'MarkerSize', 7);
h = plot(ax, 1:n, T.indice_voisines, 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
yline(ax, 1, '--', 'Color', [0.8 0.1 0.1]);
separateurs(ax, T.animal, ylim(ax));
xlim(ax, [0.4, n + 0.6]);
set(ax, 'XTick', 1:n, 'XTickLabel', T.session, 'XTickLabelRotation', 60, 'FontSize', 7, 'TickDir', 'out', ...
    'TickLabelInterpreter', 'none');
ylabel(ax, 'same colour among 3 nearest / chance');
legend(ax, h, {'true colours (grey: shifted by one position)'}, 'Location', 'northeast');
title(ax, 'B. Alignment check: neighbours share colours more than chance (should fall to 1 when shifted)', 'FontSize', 10);
box(ax, 'off');
exportgraphics(fig, fullfile(dossier, 'verification_couleurs.png'), 'Resolution', 200);
end


function separateurs(ax, animaux, yl)
% traits verticaux entre animaux, nom de l'animal en haut
[~, debut] = unique(animaux, 'stable');
debut = sort(debut);
for k = 1:numel(debut)
    if k > 1, plot(ax, (debut(k) - 0.5) * [1 1], yl, '-', 'Color', [0.4 0.4 0.4]); end
end
end
