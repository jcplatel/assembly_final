function statut = etape_soeurs_session(session, config)
%% =========================================================================
%  ÉTAPE 3d, 1er passage, UNE session (7 nov.)
%  -------------------------------------------------------------------------
%  Question : la similarité d'activité d'une paire augmente-t-elle avec sa
%  probabilité d'être sœurs ?
%  Ce passage prépare deux choses indépendantes :
%  1. COULEURS ET POSITIONS SEULEMENT (aucune activité) : comptes de paires
%     de cellules colorées par tranche de distance (et en 2D : écart en x,
%     écart en y), toutes et de même couleur, couleur par couleur. Ils
%     serviront à estimer la probabilité d'être sœurs (carte_soeurs, sur
%     toutes les sessions).
%  2. ACTIVITÉ SEULEMENT (aucune couleur) : similarité de chaque paire,
%     puis on retire la distance et le niveau propre de chaque cellule
%     (residus_paires, comme l'étape 3b) :
%       'sttc'        : STTC au repos (fenêtre de config.paires) ;
%       'sce'         : co-participation aux SCE (phi, coparticipation_sce,
%                       8 nov.) ;
%       'coassemblee' : 1 si les deux cellules sont dans la même assemblée
%                       au K retenu (étapes 2-3), sinon 0.
%  Paires dont les ROI se chevauchent exclues. Écrit <session>\soeurs\soeurs.mat.
%  =========================================================================

opts = options_soeurs(config);
dossier = fullfile(session.dossier, 'soeurs');
if ~isfolder(dossier), mkdir(dossier); end
fichier = fullfile(dossier, 'soeurs.mat');
if isfile(fichier), delete(fichier); end

pre = load(fullfile(session.dossier_source, 'pretraitement.mat'), 'Tr1b', 'Raster', 'WinRest', 'sampling_rate', 'Race');
NCell = size(pre.Raster, 1);
path = avec_separateur(session.dossier_suite2p);
[outline_gcampx, ~, gcamp_centroid, ~, ~, Coomasqx, Coomasqy] = charger_geometrie(path);
um = session.um_par_px;
if isnan(um), um = config.um_par_px_defaut; end
xy_um = gcamp_centroid * um;
dist_um = sqrt((xy_um(:, 1) - xy_um(:, 1)') .^ 2 + (xy_um(:, 2) - xy_um(:, 2)') .^ 2);
V0 = ~masques_chevauchent(Coomasqx, Coomasqy) & isfinite(dist_um);
V0(logical(eye(NCell))) = false;
haut = triu(true(NCell), 1);

%% 1. Couleurs et positions (aucune activité)
colorcell = charger_couleurs(path, session.fichier_colorcell);
colorcell = couleurs_sans_doublons(colorcell, session);    % doublons de l'étape 1 : non colorés
verifier_alignement(NCell, colorcell, outline_gcampx, session.nom);
col = colorcell(:);
colorees = col >= 1 & col <= 6;
[i, j] = find(V0 & haut & colorees & colorees');
d = dist_um(sub2ind([NCell NCell], i, j));
ci = col(i);
meme = ci == col(j);
nb = numel(opts.bornes_um) - 1;
b = discretize(d, opts.bornes_um);
C = struct();
C.N_b = accumarray(b, 1, [nb 1]);                                % paires colorées par tranche
C.S_cb = accumarray([b(meme), ci(meme)], 1, [nb 6]);            % paires de même couleur c
C.q_c = sum(C.S_cb, 1) / max(1, numel(d));                       % part des paires de même couleur c dans le champ
C.q = sum(meme) / max(1, numel(d));
% direction des bandes de couleur (7 nov., choix de l'utilisateur : les
% bandes n'ont pas la même direction d'une session à l'autre)
[C.angle, theta_sans] = direction_bandes(xy_um, col, i, j, d, ci, meme, opts);
% carte 2D alignée sur les bandes : les paires de la couleur c sont tournées
% de l'angle trouvé avec les 5 AUTRES couleurs (sinon on s'aligne en partie
% sur le bruit de ces mêmes paires et on voit toujours un faux allongement)
nk = ceil(opts.carte_max_um / opts.carte_pas_um);
dxs = xy_um(j, 1) - xy_um(i, 1);
dys = xy_um(j, 2) - xy_um(i, 2);
C.S_uv = zeros(nk);                                              % lignes = écart en travers, colonnes = le long
C.E_uv = zeros(nk);
for c = 1:6
    t = theta_sans(c);
    u = abs(dxs * cos(t) + dys * sin(t));
    v = abs(-dxs * sin(t) + dys * cos(t));
    ok2 = u < opts.carte_max_um & v < opts.carte_max_um;
    bu = floor(u / opts.carte_pas_um) + 1;
    bv = floor(v / opts.carte_pas_um) + 1;
    mc = ok2 & meme & ci == c;
    C.S_uv = C.S_uv + accumarray([bv(mc), bu(mc)], 1, [nk nk]);
    C.E_uv = C.E_uv + C.q_c(c) * accumarray([bv(ok2), bu(ok2)], 1, [nk nk]);
end

%% 2. Activité (aucune couleur) : résidus
popts = config.paires;
M = struct();
for k = 1:numel(opts.mesures)
    nom = opts.mesures{k};
    switch nom
        case 'sttc'
            if ~exist('S', 'var')
                S = similarite_paires(pre.Tr1b, pre.Raster, pre.WinRest, pre.sampling_rate, popts);
            end
            brut = S.sttc;
        case 'sce'
            brut = coparticipation_sce(pre.Race);
        case 'coassemblee'
            brut = coassemblee(session, NCell);
            if isempty(brut)
                M.(nom) = struct('E', [], 'V', [], 'brut_100', NaN, 'statut', 'pas d''assemblées au K retenu');
                continue
            end
        otherwise
            error('config.soeurs.mesures : ''sttc'', ''sce'' ou ''coassemblee'' (reçu ''%s'').', nom);
    end
    V = V0 & isfinite(brut);
    E = residus_paires(brut, V, dist_um, popts.pas_distance_um, false);
    M.(nom) = struct('E', single(E), 'V', V, 'brut_100', mean(brut(V & haut & dist_um <= 100)), 'statut', 'ok');
end
clear pre S

R = struct('nom', session.nom, 'animal', session.id_animal, 'opts', opts, 'col', col, 'xy_um', xy_um, ...
    'comptes', C, 'mesures', M);
save(fichier, '-struct', 'R', '-v7.3');
statut = sprintf('%d cellules colorées, %d paires ; mesures : %s', sum(colorees), numel(d), ...
    strjoin(opts.mesures, ', '));
end


function [A, theta_sans] = direction_bandes(xy, col, i, j, d, ci, meme, opts)
% Direction où les paires de même couleur sont en excès (COULEURS ET
% POSITIONS SEULEMENT). Chaque paire de cellules colorées à une distance
% dans opts.angle_um a une orientation phi (0-180°) ; on représente phi par
% le vecteur unité d'angle 2*phi (une orientation, pas un sens). Excès =
% somme des vecteurs des paires de même couleur - part attendue x somme des
% vecteurs de toutes les paires colorées (ce second terme retire la forme
% du champ). Angle des bandes = angle de l'excès / 2.
% theta_sans(c) : même calcul sans AUCUNE paire touchant une cellule de
% couleur c (radians) ; sert à aligner les paires de la couleur c sans
% s'aligner sur leur propre bruit.
% A.theta_deg : 0° = horizontale de l'image ; l'axe y de l'image descend.
% A.long, A.travers : paires de même couleur observées / attendues à moins
% de 30° de la direction, et à plus de 60° ; chaque couleur mesurée sur
% l'angle trouvé sans elle (contrôle croisé : sans direction commune aux
% clones, long = travers en moyenne).
% Simulation (7 nov., 40 champs) : clones allongés à 30°, erreur d'angle
% médiane 7° ; clones ronds, carte le long - travers = -0,01 +- 0,01. Un p
% par mélange des couleurs a été essayé et écarté (70 % de « direction
% nette » avec des clones ronds : les paires d'un même clone ne sont pas
% indépendantes).
rg = d >= opts.angle_um(1) & d <= opts.angle_um(2);
phi = atan2(xy(j, 2) - xy(i, 2), xy(j, 1) - xy(i, 1));
z = exp(2i * phi);
cj = col(j);
n_tot = numel(d);
V = sum(z(rg & meme)) - sum(meme) / max(1, n_tot) * sum(z(rg));
A = struct();
A.theta_deg = mod(angle(V) / 2 * 180 / pi, 180);
theta_sans = zeros(1, 6);
obs = [0 0];
att = [0 0];
for c = 1:6
    sans = ci ~= c & cj ~= c;                                       % paires sans cellule de couleur c
    ms = meme & sans;
    Vs = sum(z(rg & ms)) - sum(ms) / max(1, sum(sans)) * sum(z(rg & sans));
    theta_sans(c) = angle(Vs) / 2;
    delta = abs(mod(phi - theta_sans(c) + pi / 2, pi) - pi / 2);    % écart à la direction, 0-90°
    mc = meme & ci == c;
    q_c = sum(mc) / max(1, n_tot);
    le_long = rg & delta < pi / 6;
    travers = rg & delta > pi / 3;
    obs = obs + [sum(mc & le_long), sum(mc & travers)];
    att = att + q_c * [sum(le_long), sum(travers)];
end
A.long = obs(1) / max(eps, att(1));
A.travers = obs(2) / max(eps, att(2));
A.n_paires_meme = sum(rg & meme);
end


function A = coassemblee(session, NCell)
% 1 si deux cellules sont dans la même assemblée au K retenu, sinon 0 ; [] si rien
A = [];
asm = assemblees_retenues(session, NCell);
if isempty(asm), return; end
A = zeros(NCell);
for a = 1:numel(asm)
    A(asm{a}, asm{a}) = 1;
end
A(logical(eye(NCell))) = NaN;
end
