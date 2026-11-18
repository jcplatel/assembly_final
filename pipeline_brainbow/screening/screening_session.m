function R = screening_session(dossier_s, sess, G, asm)
%SCREENING_SESSION  Toutes les variantes d'assemblées d'UNE session, à partir
%   de l'étape 1 déjà faite (13 nov., screening_reglages). AUCUNE COULEUR.
%
%   dossier_s : <étape 1>\<session> (pretraitement.mat, k<K>\assemblees.mat)
%   sess      : struct de la session (nom, id_animal, graine)
%   G         : grille (grille_screening)
%   asm       : config.asm de base (avec NShuf)
%
%   Le k-means et ses surrogates ne sont PAS refaits : ils sont lus dans
%   k<K>\assemblees.mat (IDX2, sCl, sCl_rnd). Sont recalculés, comme le
%   ferait l'étape 1 avec ces réglages :
%    - les clusters de SCE retenus selon chaque règle de G.clusters ;
%    - les p des cellules pour TOUS les clusters (même tirage que
%      RACEORTHO_FDR2 : mêmes graines, donc mêmes valeurs) ;
%    - les membres selon chaque règle de G.membres, en version 1 (une
%      assemblée par cellule) et 2 (cellules partagées) ;
%    - les activations des assemblées dans les SCE (un seul tirage par SCE
%      pour toutes les variantes : mêmes graines que l'étape 1) ;
%    - les consensus entre K (consensus_assemblees), et leurs activations
%      (mêmes graines que consensus_session).
%   Avec les réglages de base, les assemblées en version 1 sont exactement
%   celles de l'étape 1 et de l'étape 2.
%
%   R.ens(e)   un jeu d'assemblées : cellules {1 x nA}, n_SCE (1 x nA :
%              activations significatives), type ('K' ou 'consensus'), K
%   R.map_K(iK, iRegle, iMembres, iVersion)      numéro du jeu (0 = aucune assemblée)
%   R.map_cons(iRegle, iMembres, iVersion, iChoix) idem pour les consensus
%   R.NClOK(iK, iRegle), R.exces(iK, iRegle)     clusters retenus, silhouette en excès

P = load(fullfile(dossier_s, 'pretraitement.mat'), 'Race');
Race = double(P.Race);
[NCell, NRace] = size(Race);
nK = numel(G.liste_K);
nR = size(G.clusters, 1);
nM = numel(G.membres);
nV = numel(G.versions);
nC = numel(G.choix_K);

R = struct();
R.nom_session = sess.nom;
R.id_animal   = sess.id_animal;
R.NCell       = NCell;
R.NRace       = NRace;
R.liste_K     = G.liste_K;
R.K_present   = false(1, nK);
R.NClOK       = nan(nK, nR);
R.exces       = nan(nK, nR);
R.map_K       = zeros(nK, nR, nM, nV);
R.map_cons    = zeros(nR, nM, nV, nC);
ens = struct('type', {}, 'K', {}, 'cellules', {}, 'n_SCE', {});

% réglages des membres de chaque variante
asm_m = cell(1, nM);
for m = 1:nM
    asm_m{m} = fusion_struct(asm, G.membres(m).modif);
end

%% ---------- 1. Chaque K -----------------------------------------------------
for iK = 1:nK
    K = G.liste_K(iK);
    fK = fullfile(dossier_s, sprintf('k%d', K), 'assemblees.mat');
    if ~isfile(fK), continue; end
    R.K_present(iK) = true;
    noms = intersect({'IDX2', 'sCl', 'sCl_rnd'}, who('-file', fK));
    if numel(noms) < 3, continue; end            % trop peu de SCE pour ce K : aucun cluster
    A = load(fK, 'IDX2', 'sCl', 'sCl_rnd');
    if isempty(A.sCl) || isempty(A.sCl_rnd), continue; end
    IDX2 = A.IDX2(:)';
    sCl = A.sCl(:)';
    nCl = numel(sCl);

    % clusters retenus selon chaque règle
    for r = 1:nR
        [n_ok, seuils] = clusters_retenus(sCl, A.sCl_rnd, G.clusters{r, 1}, G.clusters{r, 2});
        R.NClOK(iK, r) = n_ok;
        R.exces(iK, r) = sum(sCl(1:n_ok) - seuils(1:n_ok));
    end
    n_u = unique(R.NClOK(iK, :));
    n_u = n_u(n_u >= 1);
    if isempty(n_u), continue; end

    % p des cellules pour tous les clusters (graine de l'étape 1 :
    % SCE_clustering2(..., session.graine + K), RACEORTHO_FDR2(..., graine + 2))
    graine_K = sess.graine + K + 2;
    [pv, CellP, CellR, Nact] = p_cellules_clusters(Race, IDX2, max(n_u), asm, graine_K);

    % candidats de chaque (nombre de clusters, règle des membres, version)
    cand = struct('n', {}, 'm', {}, 'v', {}, 'M', {});
    for n = n_u
        for m = 1:nM
            for v = 1:nV
                M = candidats(pv(:, 1:n), CellP(:, 1:n), CellR(:, 1:n), Nact, asm_m{m}, G.versions{v});
                cand(end+1) = struct('n', n, 'm', m, 'v', v, 'M', M); %#ok<AGROW>
            end
        end
    end
    % activations : un seul tirage par SCE pour tous les candidats
    finals = activations_jointes(Race, cand, asm, asm_m, graine_K);
    for i = 1:numel(cand)
        if isempty(finals(i).cellules), continue; end
        ens(end+1) = struct('type', 'K', 'K', K, 'cellules', {finals(i).cellules}, ...
            'n_SCE', finals(i).n_SCE); %#ok<AGROW>
        for r = find(R.NClOK(iK, :) == cand(i).n)
            R.map_K(iK, r, cand(i).m, cand(i).v) = numel(ens);
        end
    end
end

%% ---------- 2. Consensus entre K --------------------------------------------
i_cons = find(strcmp({G.choix_K.type}, 'consensus'));
Rb = double(Race > 0);                            % comme consensus_session
cand = struct('n', {}, 'm', {}, 'v', {}, 'M', {});
cible = zeros(0, 4);                              % (r, m, v, c) de chaque candidat
cles = containers.Map('KeyType', 'char', 'ValueType', 'double');
ou = [];                                          % candidat de chaque ligne de cible
for c = i_cons
    ch = G.choix_K(c);
    [dans, iKs] = ismember(ch.liste_K, G.liste_K);
    iKs = iKs(dans);
    for r = 1:nR
        for m = 1:nM
            for v = 1:nV
                jeux = reshape(R.map_K(iKs, r, m, v), 1, []);
                % même entrée = même consensus : calculé une seule fois
                cle = sprintf('%s|%g|%g|%s|%d|%s', sprintf('%d,', jeux), ch.seuil, ch.min_K, ...
                    ch.normalisation, asm_m{m}.MIN_CELL_SIZE, sprintf('%d', R.K_present(iKs)));
                cle = sprintf('%s|%g|%g|%d', cle, valeur(asm_m{m}, 'FDR_Q_activations', asm_m{m}.FDR_Q), ...
                    asm_m{m}.FDR_Q, valeur(asm_m{m}, 'min_activations_asm', 1));
                if isKey(cles, cle)
                    j = cles(cle);
                else
                    asm_K = cell(1, numel(iKs));
                    for q = 1:numel(iKs)
                        if ~R.K_present(iKs(q))
                            asm_K{q} = [];                       % K non calculé
                        elseif jeux(q) == 0
                            asm_K{q} = cell(1, 0);               % K calculé, sans assemblée
                        else
                            asm_K{q} = ens(jeux(q)).cellules;
                        end
                    end
                    o = struct('seuil', ch.seuil, 'min_K', ch.min_K, 'normalisation', ch.normalisation, ...
                        'taille_min', asm_m{m}.MIN_CELL_SIZE + 1);
                    o.assemblees_K = asm_K;
                    C = consensus_assemblees('', G.liste_K(iKs), NCell, o);
                    M = false(NCell, numel(C.assemblees));
                    if numel(C.K_utilises) < 2, M = false(NCell, 0); end    % il faut au moins 2 K
                    for a = 1:size(M, 2)
                        M(C.assemblees{a}, a) = true;
                    end
                    cand(end+1) = struct('n', 0, 'm', m, 'v', v, 'M', M); %#ok<AGROW>
                    j = numel(cand);
                    cles(cle) = j;
                end
                cible(end+1, :) = [r, m, v, c]; %#ok<AGROW>
                ou(end+1) = j; %#ok<AGROW>
            end
        end
    end
end
if ~isempty(cand)
    % graine de consensus_session : activations_assemblees(..., session.graine + 7)
    finals = activations_jointes(Rb, cand, asm, asm_m, sess.graine + 7);
    num = zeros(1, numel(cand));
    for i = 1:numel(cand)
        if isempty(finals(i).cellules), continue; end
        ens(end+1) = struct('type', 'consensus', 'K', 0, 'cellules', {finals(i).cellules}, ...
            'n_SCE', finals(i).n_SCE); %#ok<AGROW>
        num(i) = numel(ens);
    end
    for i = 1:size(cible, 1)
        R.map_cons(cible(i, 1), cible(i, 2), cible(i, 3), cible(i, 4)) = num(ou(i));
    end
end
R.ens = ens;
end


%% =========================================================================
function [NClOK, seuils_rang] = clusters_retenus(sCl, sCl_rnd, regle, p_seuil)
% Copie de la partie 2 de SCE_clustering2 (clusters significatifs).
NCl = numel(sCl);
switch regle
    case 'rang'
        seuils_rang = prctile(sCl_rnd, p_seuil, 1);
        ok = sCl(:)' > seuils_rang(:)';
        NClOK = find(~ok, 1) - 1;
        if isempty(NClOK), NClOK = NCl; end
    case 'max'
        threshold = prctile(sCl_rnd(:, 1)', p_seuil);
        NClOK = sum(sCl > threshold);
        seuils_rang = repmat(threshold, 1, NCl);
    case 'tous'
        threshold = prctile(sCl_rnd(:), p_seuil);
        NClOK = sum(sCl > threshold);
        seuils_rang = repmat(threshold, 1, NCl);
    otherwise
        error('Règle des clusters ''%s'' inconnue (''rang'', ''max'' ou ''tous'').', regle);
end
seuils_rang = seuils_rang(:)';
end


function [pv, CellP, CellR, Nact] = p_cellules_clusters(Race, IDX2, nCl, asm, graine)
% p de chaque cellule pour les nCl premiers clusters de SCE : copie de
% l'étape 1 de RACEORTHO_FDR2 (mêmes tirages : la p d'un cluster ne dépend
% pas du nombre de clusters comptés).
[NCell, NRace] = size(Race);
CellP   = zeros(NCell, nCl);
CellR   = zeros(NCell, nCl);
nSCE_cl = zeros(1, nCl);
for i = 1:nCl
    CellP(:, i) = sum(Race(:, IDX2 == i), 2);
    nSCE_cl(i)  = sum(IDX2 == i);
    if nSCE_cl(i) > 0
        CellR(:, i) = CellP(:, i) / nSCE_cl(i);
    end
end
Nact = sum(Race, 2);
switch asm.null_cellules
    case 'uniforme'
        pv = hygecdf(CellP - 1, NRace, repmat(nSCE_cl, NCell, 1), repmat(Nact, 1, nCl), 'upper');
    case 'pondere'
        taille_SCE = sum(Race, 1);
        NShuf = asm.NShuf;
        pv = ones(NCell, nCl);
        parfor j = 1:NCell
            ligne = ones(1, nCl);
            if Nact(j) > 0
                fixer_flux_aleatoire(graine, j, NCell);
                RClr = tirages_ponderes(taille_SCE, Nact(j), NShuf, IDX2, nCl);
                ligne = ((sum(RClr >= CellP(j, :)', 2) + 1) / (NShuf + 1))';
            end
            pv(j, :) = ligne;
        end
    otherwise
        error('null_cellules doit valoir ''pondere'' ou ''uniforme''.');
end
pv(Nact == 0, :) = 1;
end


function M = candidats(pv, CellP, CellR, Nact, asm, version)
% Assemblées candidates (NCell x nCand, logique) : étapes 1.3 et 2 de
% RACEORTHO_FDR2 ('v1') ou de RACEORTHO_V2 ('v2').
[NCell, n] = size(pv);
famille = Nact > 0;
L = false(NCell, n);
switch valeur(asm, 'fdr_famille', 'globale')
    case 'globale'
        p_fam = pv(famille, :);
        sig = bh_fdr_mask(p_fam(:)', asm.FDR_Q);
        L(famille, :) = reshape(sig, size(p_fam));
    case 'par_cellule'
        for j = find(famille)'
            L(j, :) = bh_fdr_mask(pv(j, :), asm.FDR_Q);
        end
    otherwise
        error('fdr_famille doit valoir ''globale'' ou ''par_cellule''.');
end
L = L & (CellP >= valeur(asm, 'min_participations', 1));
CellCl = double(L');                         % [n x NCell]
MIN_CELL_SIZE = asm.MIN_CELL_SIZE;
taille_cl = sum(CellCl, 2);
switch version
    case 'v1'
        CellCl(taille_cl <= MIN_CELL_SIZE, :) = 0;
        A2 = find(sum(CellCl, 1) >= 2);      % cellules de plusieurs clusters : celui où elles participent le plus
        for i = A2
            sig_i = find(CellCl(:, i));
            [~, k] = max(CellR(i, sig_i));
            CellCl(:, i) = 0;
            CellCl(sig_i(k), i) = 1;
        end
        garde = sum(CellCl, 2) > MIN_CELL_SIZE;
    case 'v2'
        garde = taille_cl > MIN_CELL_SIZE;
    otherwise
        error('Version ''%s'' inconnue (''v1'' ou ''v2'').', version);
end
M = CellCl(garde, :)' > 0;
end


function finals = activations_jointes(Race, cand, asm, asm_m, graine)
% Activations de chaque assemblée candidate dans chaque SCE (étape 3 de
% RACEORTHO_FDR2 / activations_assemblees), pour tous les jeux de candidats
% à la fois : le tirage au hasard d'un SCE est fait une fois et sert à
% toutes les colonnes. finals(i).cellules, finals(i).n_SCE : assemblées
% gardées du jeu i (au moins min_activations_asm activations significatives).
[NCell, NRace] = size(Race);
finals = struct('cellules', cell(1, numel(cand)), 'n_SCE', []);
Mtout = false(NCell, 0);
for i = 1:numel(cand)
    Mtout = [Mtout, cand(i).M]; %#ok<AGROW>
end
if isempty(Mtout), return; end
[Mu, ~, ic] = unique(double(Mtout)', 'rows');     % colonnes différentes seulement
Mu = Mu' > 0;
U = size(Mu, 2);
obs = double(Mu)' * Race;                          % U x NRace : membres actifs dans chaque SCE
Nrnd_SCE = sum(Race, 1);
Nact = sum(Race, 2);
switch asm.null_cellules
    case 'uniforme'
        taille = sum(Mu, 1)';
        pU = hygecdf(obs - 1, NCell, repmat(taille, 1, NRace), repmat(Nrnd_SCE, U, 1), 'upper');
    case 'pondere'
        activite = Nact';
        NShuf = asm.NShuf;
        Ms = sparse(double(Mu));
        pU = ones(U, NRace);
        parfor j = 1:NRace
            colonne = ones(U, 1);
            if Nrnd_SCE(j) > 0
                fixer_flux_aleatoire(graine + 1, j, NRace);
                cles = rand(NShuf, NCell) .^ (1 ./ activite);    % comme tirages_ponderes
                [~, idx] = maxk(cles, Nrnd_SCE(j), 2);
                Dt = sparse(repmat((1:NShuf)', 1, Nrnd_SCE(j)), idx, 1, NShuf, NCell);
                comptes = full(Dt * Ms);                           % NShuf x U
                colonne = (sum(comptes >= obs(:, j)', 1)' + 1) / (NShuf + 1);
            end
            pU(:, j) = colonne;
        end
    otherwise
        error('null_cellules doit valoir ''pondere'' ou ''uniforme''.');
end
pU(:, Nrnd_SCE == 0) = 1;
fin = 0;
for i = 1:numel(cand)
    nc = size(cand(i).M, 2);
    cols = ic(fin + 1:fin + nc);
    fin = fin + nc;
    if nc == 0, continue; end
    a = asm_m{cand(i).m};
    PCl = bh_colonnes(pU(cols, :), valeur(a, 'FDR_Q_activations', a.FDR_Q));
    PCl(:, Nrnd_SCE == 0) = false;
    n_SCE = sum(PCl, 2)';
    garde = find(n_SCE >= valeur(a, 'min_activations_asm', 1));
    finals(i).cellules = arrayfun(@(g) find(cand(i).M(:, g))', garde, 'UniformOutput', false);
    finals(i).n_SCE = n_SCE(garde);
end
end


function sig = bh_colonnes(P, q)
% Benjamini-Hochberg dans chaque colonne de P (même résultat que
% bh_fdr_mask appliqué colonne par colonne).
[m, n] = size(P);
sp = sort(P, 1);
dessous = sp <= repmat((1:m)' / m * q, 1, n);
kmax = max(dessous .* repmat((1:m)', 1, n), [], 1);
crit = -inf(1, n);
ok = find(kmax > 0);
crit(ok) = sp(sub2ind([m n], kmax(ok), ok));
sig = P <= repmat(crit, m, 1);
end


function v = valeur(s, champ, defaut)
if isfield(s, champ) && ~isempty(s.(champ)), v = s.(champ); else, v = defaut; end
end


function s = fusion_struct(s, modif)
noms = fieldnames(modif);
for k = 1:numel(noms)
    s.(noms{k}) = modif.(noms{k});
end
end
