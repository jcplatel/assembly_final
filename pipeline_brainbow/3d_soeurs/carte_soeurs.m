function K = carte_soeurs(sessions, config)
%CARTE_SOEURS  Probabilité qu'une paire de même couleur soit une paire de
%   sœurs, selon sa couleur et sa distance. COULEURS ET POSITIONS SEULEMENT
%   (aucune activité), toutes les sessions ensemble.
%   Pour la couleur c et la tranche de distance b :
%     observé = paires de même couleur c (somme des sessions) ;
%     attendu = paires de couleurs différentes de la tranche x (paires de
%               même couleur c / paires de couleurs différentes au-delà de
%               opts.pi_max_um), session par session, puis somme ;
%     rapport r = observé / attendu ; pi = 1 - 1 / r (0 si r <= 1).
%   Idée : les paires de couleurs différentes ne sont jamais sœurs (un clone
%   a une seule couleur) et se répartissent dans l'espace comme les paires de
%   même couleur non sœurs (clones différents) ; au-delà de 200 µm il n'y a
%   pas de sœurs (choix de l'utilisateur) : ce rapport y donne le niveau du
%   hasard. L'excès au-dessus de ce niveau = les sœurs.
%   (Remplace le 7 nov., à la demande de l'utilisateur, l'ancien attendu =
%   paires colorées x part des paires de même couleur c dans tout le champ :
%   cette part contient les sœurs elles-mêmes, et pi était très sous-estimé,
%   surtout pour les couleurs rares. Simulation : couleur rare en un seul
%   clone, 100 % de sœurs : ancien 86 % à 0-25 µm et 26 % à 150-200 µm,
%   nouveau 98 % et 87 % ; couleur fréquente, 55 % de sœurs à 0-25 µm :
%   ancien 31 %, nouveau 52 %.)
%   Hypothèses : un clone a une seule couleur (les erreurs de classement
%   font baisser pi) ; pas de sœurs au-delà de 200 µm (les clones de 500 µm
%   de long font un peu baisser pi) ; tout l'excès est clonal (un gradient
%   technique de couleur dans le champ le gonflerait).
%   Couleur avec moins de opts.min_observees paires de même couleur dans une
%   tranche : pi de toutes les couleurs (trop peu de paires pour estimer).
%   Une couleur rare a souvent un attendu < 5 à courte distance : c'est
%   normal (presque toutes ses paires proches sont sœurs), pas un manque de
%   données.
%   pi = 0 au-delà de opts.pi_max_um (200 µm, choix de l'utilisateur le
%   7 nov. : taille maximale d'un clone ; l'excès de paires de même couleur
%   très éloignées, ex. G à > 400 µm sur 444119, ne peut pas venir de sœurs).
%   pi_brut_b, pi_brut_cb : avant cette limite (figure seulement).
%   Aussi la carte 2D alignée sur les bandes de couleur de chaque session
%   (écart le long des bandes, écart en travers ; chaque couleur tournée de
%   l'angle trouvé avec les autres couleurs, etape_soeurs_session), pour
%   voir la forme des clones (descriptif), et le tableau des directions
%   (angles_bandes.csv).
%   Écrit <PathSave>\soeurs_agregation\carte_soeurs.mat.

opts = options_soeurs(config);
nb = numel(opts.bornes_um) - 1;
nk = ceil(opts.carte_max_um / opts.carte_pas_um);
S_cb = zeros(nb, 6); E_cb = zeros(nb, 6);
S_uv = zeros(nk); E_uv = zeros(nk);
n_sessions = 0;
ang = cell(0, 6);
loin = opts.bornes_um(1:nb) >= opts.pi_max_um;     % tranches qui commencent au-delà de la limite
sans_loin = {};
for s = 1:height(sessions)
    f = fullfile(sessions.dossier{s}, 'soeurs', 'soeurs.mat');
    if ~isfile(f), continue; end
    R = load(f, 'comptes');
    C = R.comptes;
    D_b = C.N_b - sum(C.S_cb, 2);                      % paires de couleurs différentes par tranche
    D_loin = sum(D_b(loin));
    if D_loin > 0
        k_c = sum(C.S_cb(loin, :), 1) / D_loin;        % niveau du hasard de chaque couleur dans ce champ
        S_cb = S_cb + C.S_cb;
        E_cb = E_cb + D_b * k_c;
    else
        sans_loin{end+1} = sessions.nom{s}; %#ok<AGROW>
    end
    S_uv = S_uv + C.S_uv;
    E_uv = E_uv + C.E_uv;
    a = C.angle;
    ang(end+1, :) = {sessions.id_animal{s}, sessions.nom{s}, a.theta_deg, a.long, a.travers, ...
        a.n_paires_meme}; %#ok<AGROW>
    n_sessions = n_sessions + 1;
end
if ~isempty(sans_loin)
    warning('pi : %d session(s) sans paire au-delà de %g µm, non utilisées pour pi : %s', numel(sans_loin), ...
        opts.pi_max_um, strjoin(sans_loin, ', '));
end
K = struct('opts', opts, 'n_sessions', n_sessions, 'S_cb', S_cb, 'E_cb', E_cb, 'S_uv', S_uv, 'E_uv', E_uv);
K.sessions_sans_loin = sans_loin;
K.angles = cell2table(ang, 'VariableNames', {'animal', 'session', 'angle_deg', 'meme_couleur_le_long', ...
    'meme_couleur_en_travers', 'n_paires_meme'});
K.r_b = sum(S_cb, 2) ./ sum(E_cb, 2);                 % toutes couleurs
K.r_cb = S_cb ./ E_cb;                                % Inf si attendu 0 et observé > 0 : pi = 1
peu = S_cb < opts.min_observees;
r_tous = repmat(K.r_b, 1, 6);
K.r_cb(peu) = r_tous(peu);
K.pi_b = max(0, 1 - 1 ./ K.r_b);
K.pi_cb = max(0, 1 - 1 ./ K.r_cb);
K.pi_b(~isfinite(K.pi_b)) = 0;
K.pi_cb(~isfinite(K.pi_cb)) = 0;
K.pi_brut_b = K.pi_b;
K.pi_brut_cb = K.pi_cb;
K.pi_b(loin) = 0;
K.pi_cb(loin, :) = 0;
K.repli = peu;
K.r_uv = S_uv ./ E_uv;
dossier = fullfile(config.PathSave, 'soeurs_agregation');
if ~isfolder(dossier), mkdir(dossier); end
save(fullfile(dossier, 'carte_soeurs.mat'), '-struct', 'K');
writetable(K.angles, fullfile(dossier, 'angles_bandes.csv'));
end
