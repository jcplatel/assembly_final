%% =========================================================================
%  BIAIS CLONAL DES ASSEMBLÉES — analyse d'UN champ d'imagerie
%  -------------------------------------------------------------------------
%  Remplace l'ancien bloc « Test hypergéo + chi2 sur noblack » et le bloc
%  « Calcul color all cells ».
%  À lancer après l'extraction des assemblées, dans le même espace de travail.
%
%  QUESTION POSÉE POUR CHAQUE ASSEMBLÉE
%    « Cette assemblée contient-elle trop de cellules d'UNE couleur, quelle
%      qu'elle soit, par rapport à un tirage au hasard dans le champ ? »
%
%  PRINCIPE : test « min-p » par permutation (Westfall & Young)
%    1. Pour chaque assemblée et chaque couleur : p-value hypergéométrique
%       unilatérale (enrichissement), exactement comme avant.
%    2. On garde la plus petite des 6 : T_obs (la couleur la plus surprenante).
%    3. On mélange les couleurs entre les cellules du champ, les assemblées
%       restant fixes, et on refait exactement 1 et 2 -> T_perm.
%       On répète nPerm fois.
%    4. p_minp = proportion de mélanges où le hasard fait aussi bien ou mieux.
%       Comme le choix de la « meilleure couleur » est refait dans chaque
%       mélange, la correction pour les 6 couleurs est déjà incluse :
%       PAS de Bonferroni à ajouter par-dessus.
%
%    Un SEUL mélange par itération sert à TOUTES les assemblées du champ :
%    elles forment une partition (si une assemblée prend du rouge, les autres
%    en ont moins), et cette contrainte doit être conservée.
%
%  CE QUE LE SCRIPT CALCULE AUSSI
%    - p_min_atteignable : la plus petite p_minp que l'assemblée pourrait
%      obtenir dans le cas le plus favorable. Elle ne dépend que de la taille
%      de l'assemblée et de la composition du champ, jamais des couleurs
%      observées dans l'assemblée -> sert à déclarer une assemblée
%      « non testable » sans tricher.
%    - p_chi2 : test de la composition globale de l'assemblée (toutes les
%      couleurs à la fois : hypergéométrique multivariée), calculé avec les
%      MÊMES mélanges. Remplace l'ancien chi2, dont la p-value (chi2cdf)
%      était fausse avec de petits effectifs attendus.
%      Attention : le chi2 divise par l'effectif attendu, il est donc très
%      sensible aux couleurs RARES (2 cellules cyan à 3 % suffisent à le
%      rendre significatif). C'est pourquoi il reste un test secondaire.
%    - p_champ : test global du champ (« les assemblées de ce champ
%      dépendent-elles de la couleur ? ») = somme des chi2 des assemblées.
%
%  VARIABLES ATTENDUES (celles de ton pipeline)
%    colorcell : couleur de chaque cellule
%                1=R 2=G 3=B 4=RG 5=RB 6=GB 7=All (blanc) 8=none
%    assembly  : cell array, assembly{a} = indices des cellules de l'assemblée a
%    nbassemblies, calcium, outline_gcampx, outline_gcampy, namefull
%    id_animal, id_champ, id_session : identifiants (voir section 0)
%  OPTIONNEL
%    score_sequence : vecteur (une valeur par assemblée), score de
%                     séquentialité défini SANS regarder la couleur
%    Race           : matrice cellules x SCE (prioritaire pour les descripteurs SCE)
%    RasterRace     : raster binaire cellules x trames (utilisé si Race manque)
%    PCl_final      : activations significatives assemblées x SCE (RACEORTHO)
%    id_K           : K de la session (pour l'enregistrer dans la table)
%
%  SORTIE
%    stats_clonales.mat (table « res ») dans namefull, à regrouper ensuite
%    sur tous les champs et animaux avec agregation_stats_clonales.m
%  =========================================================================


%% ---------- 0. PARAMÈTRES (à fixer AVANT de regarder les résultats) -----
% Identifiants de la session : à définir AVANT d'appeler ce script
% (c'est ce que fait lancer_stats_toutes_sessions.m) :
%   id_animal    : identifiant de l'animal
%   id_champ     : identifiant du champ (plan) ; LE MÊME d'un jour à l'autre
%                  si c'est le même champ de vue
%   id_session   : identifiant de cet enregistrement (date)
%   id_condition : optionnel (ex. 'sans_indice' / 'avec_indice')
if ~exist('id_animal', 'var') || ~exist('id_champ', 'var') || ~exist('id_session', 'var')
    error(['Définir id_animal, id_champ et id_session avant de lancer ' ...
        'stats_clonales_champ (sinon toutes les sessions auraient le même nom).']);
end
if ~exist('id_condition', 'var'), id_condition = 'non_precisee'; end

nPerm          = 10000;  % nombre de mélanges. Plus petite p possible = 1/(nPerm+1)
seuilstat      = 0.05;   % seuil de significativité
nmin_couleur   = 5;      % couleur exclue du test si elle a moins de cellules que ça
                         % dans le champ (décidé sur le champ, jamais sur l'assemblée)
nmin_assemblee = 3;      % taille minimale (en cellules colorées) pour tester une assemblée
reference      = 'champ';% 'champ'      : tirage parmi toutes les cellules colorées du champ
                         %                (= ton analyse précédente)
                         % 'assemblees' : tirage parmi les seules cellules colorées
                         %                appartenant à une assemblée
graine         = 1;      % graine aléatoire : mêmes résultats à chaque exécution
faire_figures  = true;

nom_couleur  = {'R','G','B','RG','RB','GB'};
couleurs_rgb = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1];
nC = 6;
% Nombre d'assemblées : compté directement dans « assembly », pour ne pas
% dépendre d'un nbassemblies resté d'une autre session ou d'un autre k.
nA = numel(assembly);
if exist('nbassemblies', 'var') && nbassemblies ~= nA
    warning(['nbassemblies (%d) ne correspond pas au nombre d''assemblées ' ...
        'dans assembly (%d). Le script utilise %d. Vérifie que assembly, ' ...
        'nbassemblies et score_sequence viennent bien de la même session et du même k.'], ...
        nbassemblies, nA, nA);
end
colorcell = colorcell(:);

% Contrôle : assembly et colorcell doivent décrire la même liste de cellules
idx_max = max(cellfun(@(x) max([x(:); 0]), assembly));
if idx_max > numel(colorcell)
    error(['Les assemblées utilisent la cellule n°%d, mais colorcell n''a que %d éléments : ' ...
        'assembly et colorcell ne décrivent pas la même liste de cellules ' ...
        '(segmentation différente, ou données de deux analyses différentes).'], ...
        idx_max, numel(colorcell));
end
if faire_figures && idx_max > min(numel(outline_gcampx), numel(outline_gcampy))
    error(['Les assemblées utilisent la cellule n°%d, mais les contours n''en ont que %d : ' ...
        'assembly et outline_gcampx/y ne décrivent pas la même liste de cellules.'], ...
        idx_max, min(numel(outline_gcampx), numel(outline_gcampy)));
end


%% ---------- 1. Composition du champ (toutes les cellules) ----------------
nbR  = sum(colorcell==1);  nbG  = sum(colorcell==2);  nbB  = sum(colorcell==3);
nbRG = sum(colorcell==4);  nbRB = sum(colorcell==5);  nbGB = sum(colorcell==6);
nball  = sum(colorcell==7);   % blanc (3 canaux) : exclu du test, comme avant
nbnone = sum(colorcell==8);   % sans couleur     : exclu du test
color_all = [nbR nbG nbB nbRG nbRB nbGB nball nbnone];

if faire_figures && isfolder(namefull)
    fig = figure('visible','off');
    ax = axes(fig);
    dessiner_camembert(ax, color_all, [], [couleurs_rgb; 1 1 1; 0 0 0]);
    titre_obj = sgtitle(fig, 'color All cells', 'FontSize', 12, 'FontWeight', 'bold');
    titre_obj.HorizontalAlignment = 'left';
    exportgraphics(fig, strcat(namefull, 'color_all_cells.png'), ...
        'Resolution', 150, 'ContentType', 'image');
    close(fig)
end


%% ---------- 2. Appartenance des cellules aux assemblées ------------------
% On liste toutes les paires (assemblée, cellule). Cela marche même si une
% cellule appartenait à deux assemblées, mais on le signale.
paires_asm = [];
paires_cell = [];
for a = 1:nA
    cells_a = unique(assembly{a}(:));
    paires_asm  = [paires_asm;  repmat(a, numel(cells_a), 1)]; %#ok<AGROW>
    paires_cell = [paires_cell; cells_a];                      %#ok<AGROW>
end
nb_appart = accumarray(paires_cell, 1, [numel(colorcell) 1]);
if any(nb_appart > 1)
    warning(['%d cellule(s) appartiennent à plusieurs assemblées : ' ...
        'le test reste valide, mais vérifie l''orthogonalisation.'], sum(nb_appart > 1));
end


%% ---------- 3. Population de référence pour le tirage --------------------
est_coloree = colorcell >= 1 & colorcell <= 6;   % 7 (blanc) et 8 (aucune) exclues
switch reference
    case 'champ'
        idx_ref = find(est_coloree);
    case 'assemblees'
        idx_ref = find(est_coloree & nb_appart > 0);
    otherwise
        error('reference doit valoir ''champ'' ou ''assemblees''.');
end
Nref    = numel(idx_ref);                       % nb de cellules dans l'urne
col_ref = colorcell(idx_ref);                   % leurs couleurs (1..6)
Kcol    = accumarray(col_ref, 1, [nC 1])';      % nb de cellules de chaque couleur
couleur_valide = Kcol >= nmin_couleur;
if Nref == 0
    error('Aucune cellule colorée dans la population de référence.');
end

% Position de chaque cellule d'assemblée dans l'urne (0 = pas dans l'urne :
% cellule blanche ou sans couleur, ignorée comme avant)
pos_dans_ref = zeros(numel(colorcell), 1);
pos_dans_ref(idx_ref) = 1:Nref;
garde = pos_dans_ref(paires_cell) > 0;
pA    = paires_asm(garde);                      % assemblée de chaque paire
pPos  = pos_dans_ref(paires_cell(garde));       % position de la cellule dans l'urne
if isempty(pA)
    error('Aucune cellule colorée dans les assemblées.');
end

n_col = accumarray(pA, 1, [nA 1]);                        % cellules colorées par assemblée
n_tot = cellfun(@(x) numel(unique(x)), assembly(:));      % taille totale de l'assemblée

% Comptes observés : Obs(a, c) = nb de cellules de couleur c dans l'assemblée a
Obs = accumarray([pA, col_ref(pPos)], 1, [nA nC]);


%% ---------- 4. Table des p-values hypergéométriques (calculée une fois) ---
% Ptab(a, c, k+1) = P(au moins k cellules de couleur c dans l'assemblée a)
% Elle ne dépend que de la taille de l'assemblée et de la composition de
% l'urne, qui ne changent pas quand on mélange : on la calcule une seule
% fois, puis chaque mélange se réduit à des recomptages.
nmax = max([n_col; 1]);
Ptab = nan(nA, nC, nmax + 1);
for a = 1:nA
    for c = 1:nC
        if couleur_valide(c)
            kv = 0:n_col(a);
            Ptab(a, c, 1:n_col(a)+1) = hygecdf(kv - 1, Nref, Kcol(c), n_col(a), 'upper');
            Ptab(a, c, 1) = 1;               % P(au moins 0 cellule) = 1
        else
            Ptab(a, c, 1:n_col(a)+1) = 1;    % couleur exclue : ne gagne jamais
        end
    end
end
% Indice linéaire de Ptab pour (a, c, k) = base(a,c) + k*nA*nC
base = reshape(1:nA*nC, nA, nC);

% Statistique observée : plus petite p-value sur les couleurs
p_brut = Ptab(base + Obs*nA*nC);                 % p hypergéo de chaque couleur (non corrigée)
T_obs  = min(p_brut, [], 2);

% Cas le plus favorable : toutes les cellules possibles de la même couleur
m_max = min(repmat(n_col, 1, nC), repmat(Kcol, nA, 1));
T_min = min(Ptab(base + m_max*nA*nC), [], 2);

% Chi2 de composition (effectifs attendus sous tirage au hasard)
Eatt = n_col * (Kcol / Nref);
Eatt(:, ~couleur_valide) = NaN;
chi2_obs = sum((Obs - Eatt).^2 ./ Eatt, 2, 'omitnan');


%% ---------- 5. Permutations : on mélange les couleurs dans l'urne --------
rng(graine);
T_perm    = zeros(nA, nPerm);
chi2_perm = zeros(nA, nPerm);
for b = 1:nPerm
    col_melange = col_ref(randperm(Nref));                     % même mélange pour toutes
    cnt = accumarray([pA, col_melange(pPos)], 1, [nA nC]);     % les assemblées du champ
    T_perm(:, b)    = min(Ptab(base + cnt*nA*nC), [], 2);
    chi2_perm(:, b) = sum((cnt - Eatt).^2 ./ Eatt, 2, 'omitnan');
end


%% ---------- 6. p-values ---------------------------------------------------
tol = 1e-10;   % tolérance numérique pour les égalités
% p_minp : fraction des mélanges où la meilleure couleur est au moins aussi
% surprenante que dans l'assemblée réelle (+1 : on compte l'observation)
p_minp            = (1 + sum(T_perm    <= T_obs * (1 + tol), 2)) / (nPerm + 1);
p_min_atteignable = (1 + sum(T_perm    <= T_min * (1 + tol), 2)) / (nPerm + 1);
p_chi2            = (1 + sum(chi2_perm >= chi2_obs * (1 - tol), 2)) / (nPerm + 1);

% Assemblée testable : assez grande, ET capable d'atteindre le seuil
% dans le cas le plus favorable. Critère indépendant des couleurs observées.
testable = n_col >= nmin_assemblee & p_min_atteignable < seuilstat;

% Test global du champ
S_obs   = sum(chi2_obs);
S_perm  = sum(chi2_perm, 1);
p_champ = (1 + sum(S_perm >= S_obs * (1 - tol))) / (nPerm + 1);
z_champ = (S_obs - mean(S_perm)) / std(S_perm);

% Couleur la plus enrichie (pour l'annotation seulement)
p_brut(:, ~couleur_valide) = NaN;
[~, idx_dom] = min(p_brut, [], 2);
couleur_dominante = nom_couleur(idx_dom)';
couleur_dominante(n_col == 0) = {'-'};
% Nombre et part des cellules de la couleur dominante, part de cellules colorées
n_couleur_dominante    = Obs(sub2ind(size(Obs), (1:nA)', idx_dom(:)));
part_couleur_dominante = n_couleur_dominante ./ n_col;
part_couleur_dominante(n_col == 0) = NaN;
part_colorees          = n_col ./ n_tot;

% Score de séquentialité, s'il existe
if exist('score_sequence', 'var') && numel(score_sequence) == nA
    sc = score_sequence(:);
else
    sc = nan(nA, 1);
end


%% ---------- 6b. Descripteurs des assemblées (exploration) ---------------
% Ne servent PAS au test clonal : ils décrivent les assemblées pour pouvoir
% comparer ensuite les biaisées et les autres.
n_cell_champ = numel(colorcell);
pct_reseau   = 100 * n_tot / n_cell_champ;      % taille en % des cellules du champ
if exist('id_K', 'var'), K_val = id_K; else, K_val = NaN; end

% Étalement dans le champ : distance moyenne entre les cellules de
% l'assemblée (pixels), et son percentile parmi 1 000 groupes aléatoires de
% même taille tirés dans le champ (petit percentile = assemblée compacte).
% ratio_etalement = D(assemblée) / moyenne des D tirés au hasard : taille
% d'effet (1 = dispersée comme au hasard, < 1 = regroupée), peu dépendante
% du nombre de cellules, contrairement au percentile (qui est une p-value).
% Le tirage se fait parmi toutes les cellules du champ (ROI suite2p).
etalement_px    = nan(nA, 1);
compacite_pct   = nan(nA, 1);
ratio_etalement = nan(nA, 1);
if exist('outline_gcampx', 'var') && exist('outline_gcampy', 'var')
    nc = min([numel(outline_gcampx), numel(outline_gcampy), n_cell_champ]);
    centres = nan(nc, 2);
    for n = 1:nc
        centres(n, :) = [mean(cell2mat(outline_gcampx(n))), mean(cell2mat(outline_gcampy(n)))];
    end
    ok_c = find(all(~isnan(centres), 2));
    rng(graine + 1);
    for a = 1:nA
        cells_a = unique(assembly{a}(:));
        cells_a = cells_a(ismember(cells_a, ok_c));
        m_a = numel(cells_a);
        if m_a < 3, continue; end
        etalement_px(a) = mean(pdist(centres(cells_a, :)));
        d_alea = zeros(1000, 1);
        for b = 1:1000
            d_alea(b) = mean(pdist(centres(ok_c(randperm(numel(ok_c), m_a)), :)));
        end
        compacite_pct(a)   = 100 * mean(d_alea <= etalement_px(a));
        ratio_etalement(a) = etalement_px(a) / mean(d_alea);
    end
end

% SCE : nombre total de SCE, nombre d'activations significatives
% (RACEORTHO, si PCl_final est fourni), nombre de SCE où l'assemblée est la
% plus active (en fraction de ses cellules), et fraction moyenne de ses
% cellules actives dans ces SCE.
% Matrice cellules x SCE utilisée :
%   - Race (cellules x SCE) si elle est fournie ;
%   - sinon reconstruite à partir de RasterRace (cellules x TRAMES, 1 si la
%     cellule participe à un SCE à cette trame) : chaque suite de trames
%     consécutives avec de l'activité = un SCE. Deux SCE accolés sans trame
%     vide entre eux sont alors comptés comme un seul.
n_SCE_total     = NaN;
n_SCE_signif    = nan(nA, 1);
n_SCE_dominante = nan(nA, 1);
frac_active_moy = nan(nA, 1);
source_SCE      = 'aucune';
if exist('PCl_final', 'var') && size(PCl_final, 1) == nA
    n_SCE_signif = sum(PCl_final > 0, 2);
end
idx_max_asm = max(cellfun(@(x) max([x(:); 0]), assembly));
S_sce = [];
if exist('Race', 'var') && ~isempty(Race) && size(Race, 1) >= idx_max_asm
    S_sce = Race > 0;
    source_SCE = 'Race';
elseif exist('RasterRace', 'var') && ~isempty(RasterRace) && size(RasterRace, 1) >= idx_max_asm
    Rr = RasterRace(:, :, 1) > 0;
    act = any(Rr, 1);
    debut = find(diff([false act]) == 1);
    fin   = find(diff([act false]) == -1);
    S_sce = false(size(Rr, 1), numel(debut));
    for e = 1:numel(debut)
        S_sce(:, e) = any(Rr(:, debut(e):fin(e)), 2);
    end
    source_SCE = 'RasterRace (trames regroupées)';
elseif exist('Race', 'var') || exist('RasterRace', 'var')
    warning('Race / RasterRace ne correspondent pas aux indices des assemblées : descripteurs SCE non calculés.');
end
if ~isempty(S_sce)
    n_SCE_total = size(S_sce, 2);
    frac = zeros(nA, n_SCE_total);
    for a = 1:nA
        cells_a = unique(assembly{a}(:));
        if ~isempty(cells_a)
            frac(a, :) = mean(S_sce(cells_a, :), 1);
        end
    end
    [vmax, lab] = max(frac, [], 1);
    lab(vmax == 0) = 0;
    for a = 1:nA
        n_SCE_dominante(a) = sum(lab == a);
        if any(lab == a)
            frac_active_moy(a) = mean(frac(a, lab == a));
        end
    end
end


%% ---------- 7. Table de résultats et sauvegarde --------------------------
res = table();
res.animal            = repmat({id_animal},    nA, 1);
res.champ             = repmat({id_champ},     nA, 1);
res.session           = repmat({id_session},   nA, 1);
res.condition         = repmat({id_condition}, nA, 1);
res.assemblee         = (1:nA)';
res.n_total           = n_tot;
res.n_colorees        = n_col;
res.testable          = testable;
res.p_min_atteignable = p_min_atteignable;
res.p_minp            = p_minp;      % TEST PRINCIPAL (corrigé pour les 6 couleurs)
res.p_chi2            = p_chi2;      % composition globale (secondaire)
res.couleur_dominante = couleur_dominante;
res.n_couleur_dominante    = n_couleur_dominante;     % cellules de la couleur dominante
res.part_couleur_dominante = part_couleur_dominante;  % ... parmi les cellules colorées
res.part_colorees          = part_colorees;           % cellules colorées / toutes les cellules
res.compte            = Obs;         % nA x 6 : R G B RG RB GB
res.p_brut            = p_brut;      % nA x 6 : p hypergéo NON corrigées (affichage)
res.score_sequence    = sc;
res.p_champ           = repmat(p_champ, nA, 1);
res.z_champ           = repmat(z_champ, nA, 1);
res.reference         = repmat({reference}, nA, 1);
res.nPerm             = repmat(nPerm, nA, 1);
% Descripteurs (exploration, section 6b)
res.K                 = repmat(K_val, nA, 1);
res.n_cellules_champ  = repmat(n_cell_champ, nA, 1);
res.pct_reseau        = pct_reseau;
res.etalement_px      = etalement_px;
res.compacite_pct     = compacite_pct;
res.ratio_etalement   = ratio_etalement;
res.n_SCE_total       = repmat(n_SCE_total, nA, 1);
res.n_SCE_signif      = n_SCE_signif;
res.n_SCE_dominante   = n_SCE_dominante;
res.frac_active_moy   = frac_active_moy;
res.source_SCE        = repmat({source_SCE}, nA, 1);

if isfolder(namefull)
    save(strcat(namefull, 'stats_clonales.mat'), 'res', 'Kcol', 'Nref');
else
    warning(['namefull (%s) n''est pas un dossier existant : ' ...
        'stats_clonales.mat N''A PAS été enregistré.'], namefull);
end

% Résumé dans la console
fprintf('\n=== %s / %s / %s : %d assemblées, %d cellules colorées dans l''urne ===\n', ...
    id_animal, id_champ, id_session, nA, Nref);
if any(~couleur_valide)
    fprintf('Couleurs exclues (moins de %d cellules) : %s\n', nmin_couleur, ...
        strjoin(nom_couleur(~couleur_valide), ', '));
end
fprintf('Assemblées testables : %d / %d\n', sum(testable), nA);
fprintf('Testables avec p_minp < %.2f : %d  (attendu par hasard : %.1f)\n', ...
    seuilstat, sum(p_minp(testable) < seuilstat), seuilstat * sum(testable));
fprintf('Test global du champ : p = %.4f  (z = %.2f)\n\n', p_champ, z_champ);
disp(res(:, {'assemblee','n_colorees','testable','p_minp','p_chi2','couleur_dominante'}))
% NB : ne pas conclure ici sur la proportion d'assemblées biaisées. La
% conclusion se fait sur tous les champs, au niveau des animaux
% (agregation_stats_clonales.m).


%% ---------- 8. Figures par assemblée -------------------------------------
if faire_figures && isfolder(namefull)
    for a = 1:nA
        if testable(a), suffixe = ''; else, suffixe = ' | non testable'; end

        % --- Carte de l'assemblée sur l'image calcique ---
        fig = figure('visible','off');
        colormap(gray); imagesc(calcium); hold on
        for n = unique(assembly{a}(:))'
            cc = colorcell(n);
            if cc >= 1 && cc <= 6, col = couleurs_rgb(cc, :); else, col = [0 0 0]; end
            x = cell2mat(outline_gcampx(n)); y = cell2mat(outline_gcampy(n));
            kk = convhull(x, y);
            patch(x(kk), y(kk), col, 'FaceAlpha', 1, 'EdgeColor', 'w', 'LineWidth', 1)
        end
        hold off; axis off; axis image;
        title(sprintf('Assemblée %d : %d cellules (%d colorées) | p_{min-p} = %.4f%s', ...
            a, n_tot(a), n_col(a), p_minp(a), suffixe), 'FontSize', 12);
        exportgraphics(fig, strcat(namefull, ['assembly', num2str(a), '.png']), ...
            'Resolution', 150, 'Padding', 30, 'ContentType', 'image');
        close(fig)

        % --- Camembert des couleurs de l'assemblée ---
        if n_col(a) > 0
            labels = cell(1, nC);
            for c = 1:nC
                if couleur_valide(c)
                    labels{c} = sprintf('%.1f%% (p brut=%.3f)', ...
                        100 * Obs(a, c) / n_col(a), p_brut(a, c));
                else
                    labels{c} = sprintf('%.1f%% (exclue)', 100 * Obs(a, c) / n_col(a));
                end
            end
            fig = figure('visible','off');
            ax = axes(fig);
            dessiner_camembert(ax, Obs(a, :), labels, couleurs_rgb);
            titre_obj = sgtitle(fig, sprintf(['Assemblée %d : %d cellules colorées | ' ...
                'p_{min-p} = %.4f | p_{chi2} = %.4f%s'], ...
                a, n_col(a), p_minp(a), p_chi2(a), suffixe), ...
                'FontSize', 11, 'FontWeight', 'bold');
            titre_obj.HorizontalAlignment = 'left';
            exportgraphics(fig, strcat(namefull, ['piechart_noblack', num2str(a), '.png']), ...
                'Resolution', 150, 'ContentType', 'image');
            close(fig)
        end
    end
end


%% ---------- Fonction locale ----------------------------------------------
function dessiner_camembert(ax, X, labels, rgb)
% Camembert dont chaque part garde SA couleur, même quand certaines parts
% valent 0 (elles sont retirées avant le tracé).
%   labels = []  -> pourcentages par défaut de pie
nz = find(X > 0);
if isempty(nz), return; end
if isempty(labels)
    h = pie(ax, X(nz), ones(1, numel(nz)));
else
    h = pie(ax, X(nz), ones(1, numel(nz)), labels(nz));
end
for k = 1:numel(nz)
    h(2*k-1).FaceColor = rgb(nz(k), :);          % part (patch)
    h(2*k).FontSize    = 8;                      % étiquette (texte)
    h(2*k).Position    = h(2*k).Position * 1.1;  % éloigne un peu le texte
end
end
