function L = resume_screening(S, ordre)
%RESUME_SCREENING  Texte récapitulatif du screening (resume_screening.txt) :
%   ce que donne le hasard après le même tri, la meilleure combinaison, les
%   10 premières, l'effet moyen de chaque réglage, les réglages conseillés
%   et les lignes de config pour les relancer sur tout le jeu de données.
%   S : sortie de screening_agregation ; ordre : combinaisons de la
%   meilleure à la moins bonne (les non éligibles à la fin).

C = S.C; G = S.G; F = S.F;
el = C.eligible > 0;
nEl = sum(el);
B = numel(S.max_nul);
L = {};
L{end+1} = sprintf('=== SCREENING DES RÉGLAGES (pipeline %s, %s) ===', S.version, S.date);
L{end+1} = sprintf('Animaux : %s ; %d sessions : %s', strjoin(S.animaux, ', '), numel(S.sessions), ...
    strjoin(S.sessions, ', '));
L{end+1} = sprintf('Prétraitements finis (%d) : %s', numel(S.pre_complets), strjoin(S.pre_complets, ', '));
if ~isempty(S.pre_incomplets)
    L{end+1} = sprintf('Prétraitements PAS finis, non comptés (%d) : %s', numel(S.pre_incomplets), ...
        strjoin(S.pre_incomplets, ', '));
end
L{end+1} = sprintf(['Combinaisons essayées : %d, dont %d éligibles (au moins %d assemblées testables ' ...
    'au total et %d par animal).'], numel(el), nEl, S.opts.min_testables, S.opts.min_par_animal);
L{end+1} = sprintf(['Mesure : EXCÈS = %% d''assemblées biaisées (p_minp < %g, test champ_globale) - %% attendu ' ...
    'par hasard, calculé par animal puis moyenné sur les animaux (en points).'], S.seuil);
% 23 nov. : calcul de p_minp (mid-p / Tarone, comme le pipeline, ou ancien)
calc = 'p_minp classique, sans Tarone (ancien calcul)';
if isfield(S, 'midp') && S.midp && isfield(S, 'tarone') && S.tarone
    calc = 'p_minp mid-p avec Tarone (comme le pipeline)';
elseif isfield(S, 'midp') && (S.midp || S.tarone)
    calc = sprintf('p_minp %s, %s', ternaire_r(S.midp, 'mid-p', 'classique'), ternaire_r(S.tarone, 'avec Tarone', 'sans Tarone'));
end
L{end+1} = sprintf('Calcul : %s ; le taux attendu suit la même règle.', calc);
L{end+1} = '';
if nEl == 0
    L{end+1} = 'Aucune combinaison éligible : rien à classer.';
    return
end

%% 1. Le hasard après le même tri
L{end+1} = '1. CE QUE DONNE LE HASARD APRÈS LE MÊME TRI';
L{end+1} = sprintf(['   Couleurs mélangées dans chaque champ (%d fois), tout le screening refait à chaque fois : ' ...
    'sans aucun effet réel,'], B);
L{end+1} = sprintf(['   le MEILLEUR réglage atteint quand même %+.1f points d''excès en médiane, et jusqu''à ' ...
    '%+.1f points dans 95 %% des cas.'], 100 * median(S.max_nul), 100 * prctile(S.max_nul, 95));
L{end+1} = '   Un meilleur réglage observé en dessous de cette valeur ne prouve rien à lui seul.';
L{end+1} = '';

%% 2. Meilleure combinaison
% Beaucoup de combinaisons donnent exactement les mêmes assemblées (ex. un
% critère de l'étape 4 qui ne retire rien) : une seule est montrée par résultat.
ke = ordre(1:nEl);
[~, ~, ic] = unique([round(1e6 * C.exces(ke)), C.n_test(ke), C.n_sig(ke), round(1e6 * C.attendu(ke))], 'rows');
premier = accumarray(ic(:), (1:nEl)', [], @min);      % première combinaison (la mieux classée) de chaque résultat
[premier, o] = sort(premier);
nb = accumarray(ic(:), 1);
rep = ke(premier);
n_meme = nb(o) - 1;
k1 = ordre(1);
L{end+1} = '2. MEILLEURE COMBINAISON';
L{end+1} = ['   ' texte_reglages(S, k1)];
if n_meme(1) > 0
    L{end+1} = sprintf('   (%d autre(s) combinaison(s) donnent exactement le même résultat : voir screening_top.csv)', ...
        n_meme(1));
end
L{end+1} = ['   ' texte_resultat(S, k1)];
L{end+1} = sprintf(['   p de cette combinaison seule : %s ; p corrigée pour les %d combinaisons essayées ' ...
    '(p_selection) : %s'], texte_p(C.p_seul(k1), B), nEl, texte_p(C.p_selection(k1), B));
if C.p_selection(k1) < 0.05
    L{end+1} = ['   -> Le meilleur réglage fait mieux que ce que le tri seul explique. ' ...
        'À confirmer sur les autres animaux (voir 7).'];
else
    L{end+1} = ['   -> On ne peut pas distinguer ce résultat de ce que donne le tri du meilleur réglage ' ...
        'parmi tous ceux essayés.'];
    L{end+1} = '      Seule la confirmation sur les autres animaux dira s''il est réel (voir 7).';
end
L{end+1} = '';

%% 3. Les 10 premiers résultats différents
L{end+1} = '3. LES 10 PREMIERS RÉSULTATS DIFFÉRENTS (les 50 premières combinaisons : screening_top.csv)';
for i = 1:min(10, numel(rep))
    k = rep(i);
    autres = '';
    if n_meme(i) > 0, autres = sprintf('  (+ %d combinaison(s) au même résultat)', n_meme(i)); end
    L{end+1} = sprintf('   %2d. %s%s', i, texte_reglages(S, k), autres); %#ok<AGROW>
    L{end+1} = sprintf('       %s ; p seule %s, p_selection %s', texte_resultat(S, k), ...
        texte_p(C.p_seul(k), B), texte_p(C.p_selection(k), B)); %#ok<AGROW>
end
L{end+1} = '';

%% 4. Effet moyen de chaque réglage
L{end+1} = '4. EFFET MOYEN DE CHAQUE RÉGLAGE';
L{end+1} = ['   Excès moyen de toutes les combinaisons éligibles qui ont ce niveau (les autres réglages ' ...
    'confondus) ; entre crochets : ce que donne le hasard (95 %).'];
L{end+1} = '   Plus fiable que la meilleure combinaison : une moyenne sur des centaines de combinaisons.';
% 4a. prétraitement, un facteur à la fois (la base est rappelée dans chaque facteur)
f1 = F(1);
[~, i_base] = ismember('base', f1.niveaux);
facteurs_pre = unique({G.pre(~strcmp({G.pre.facteur}, 'base')).facteur}, 'stable');
modif_combine = struct();
n_modif = 0;
nom_seul = '';
for jf = 1:numel(facteurs_pre)
    fa = facteurs_pre{jf};
    ip = find(strcmp({G.pre.facteur}, fa));
    ip = ip(f1.n(ip) > 0);
    if isempty(ip), continue; end
    nb = '?';
    if isfield(G, 'niveau_base') && isfield(G.niveau_base, fa), nb = G.niveau_base.(fa); end
    L{end+1} = sprintf('   %s (prétraitement)', fa); %#ok<AGROW>
    cand = ip; moy = f1.moyenne(ip)';
    if i_base > 0 && f1.n(i_base) > 0
        L{end+1} = ['      ' ligne_niveau([nb ' (base)'], f1, i_base, false)]; %#ok<AGROW>
        cand = [i_base, ip]; moy = [f1.moyenne(i_base), moy];
    end
    [~, im] = max(moy);
    for q = ip
        L{end+1} = ['      ' ligne_niveau(G.pre(q).niveau, f1, q, cand(im) == q)]; %#ok<AGROW>
    end
    if cand(im) ~= i_base || i_base == 0        % un niveau fait mieux que la base : gardé pour la combinaison
        mo = G.pre(cand(im)).modif;
        noms = fieldnames(mo);
        for z = 1:numel(noms)
            modif_combine.(noms{z}) = mo.(noms{z});
        end
        n_modif = n_modif + 1;                   % nombre de FACTEURS où un niveau fait mieux que la base
        nom_seul = G.pre(cand(im)).nom;
    end
end
% 4b. réglages des assemblées
for f = 2:numel(F)
    L{end+1} = sprintf('   %s', F(f).nom); %#ok<AGROW>
    moy = F(f).moyenne; moy(F(f).n == 0) = -Inf;
    [~, im] = max(moy);
    for q = 1:numel(F(f).niveaux)
        if F(f).n(q) == 0, continue; end
        L{end+1} = ['      ' ligne_niveau(F(f).niveaux{q}, F(f), q, q == im)]; %#ok<AGROW>
    end
end
L{end+1} = '';

%% 5. Réglages conseillés par les effets moyens
% Parmi les combinaisons classées, celle dont la SOMME des effets moyens de
% ses niveaux est la plus haute : chaque niveau y est bon en moyenne, pas
% seulement dans cette combinaison.
L{end+1} = '5. RÉGLAGES CONSEILLÉS PAR LES EFFETS MOYENS';
L{end+1} = ['   La combinaison classée dont chaque niveau est bon EN MOYENNE (somme des effets moyens du point 4 ' ...
    'la plus haute) :'];
niveaux = {C.ip, C.r, C.m, C.v, C.c, C.imc, C.ims};
prevu = zeros(numel(el), 1);
for f = 1:numel(F)
    mf = F(f).moyenne;
    mf(isnan(mf)) = -Inf;
    prevu = prevu + mf(niveaux{f});
end
prevu(~el) = -Inf;
[~, k_conseil] = max(prevu);
L{end+1} = ['   ' texte_reglages(S, k_conseil)];
L{end+1} = sprintf('      %s ; rang %d sur %d ; p seule %s, p_selection %s', texte_resultat(S, k_conseil), ...
    C.rang(k_conseil), nEl, texte_p(C.p_seul(k_conseil), B), texte_p(C.p_selection(k_conseil), B));
if n_modif == 1
    L{end+1} = sprintf(['   Prétraitement : un seul réglage fait mieux que la base en moyenne (%s), et il a déjà ' ...
        'été essayé : rien à combiner.'], nom_seul);
elseif n_modif > 1
    L{end+1} = sprintf(['   Prétraitement : %d facteurs ont un niveau meilleur que la base en moyenne (flèches du ' ...
        'point 4), mais ils n''ont pas encore été essayés ENSEMBLE.'], n_modif);
    L{end+1} = '   (Un écart plus petit que la bande du hasard entre crochets ne veut rien dire.)';
    L{end+1} = '   Pour les essayer (étape 1 à refaire pour ce réglage), copier :';
    L{end+1} = '      opts = struct();';
    L{end+1} = ['      opts.pre_supplementaires = struct(''nom'', ''combine'', ''facteur'', ''combine'', ' ...
        '''niveau'', ''meilleurs niveaux'', ''modif'', ' texte_struct(modif_combine) ');'];
    L{end+1} = '      screening_reglages(config, opts)';
else
    L{end+1} = '   Prétraitement : aucun niveau ne fait mieux que la base en moyenne.';
end
L{end+1} = '';

%% 6. Lignes de config
L{end+1} = '6. LIGNES DE CONFIG';
L{end+1} = '   a) la meilleure combinaison (point 2) :';
L = [L, cellfun(@(t) ['      ' t], lignes_config(S, k1), 'UniformOutput', false)];
if k_conseil ~= k1
    L{end+1} = '   b) les réglages conseillés par les effets moyens (point 5) :';
    L = [L, cellfun(@(t) ['      ' t], lignes_config(S, k_conseil), 'UniformOutput', false)];
end
L{end+1} = '';

%% 7. La suite
L{end+1} = '7. LA SUITE';
L{end+1} = ['   Ces réglages ont été choisis EN REGARDANT les couleurs de ' strjoin(S.animaux, ', ') '.'];
L{end+1} = '   Le test qui compte est celui des AUTRES animaux et des sessions avec indice, avec ces réglages';
L{end+1} = '   fixés à l''avance. Dans ce test, les animaux du screening sont à mettre à part (ou à montrer';
L{end+1} = '   séparément) : leur résultat est gonflé par le tri.';
L{end+1} = '   Le screening ne teste que le mélange des couleurs dans tout le champ (champ_globale) : les contrôles';
L{end+1} = '   blocs et territoire seront faits par le pipeline normal (étapes 3 et 4) sur le réglage retenu.';
end


%% =========================================================================
function t = texte_reglages(S, k)
C = S.C; G = S.G;
t = sprintf('prétraitement %s | clusters %s %g | %s | %s | %s | >= %g cellules, >= %g SCE', ...
    G.pre(C.ip(k)).nom, G.clusters{C.r(k), 1}, G.clusters{C.r(k), 2}, G.membres(C.m(k)).nom, ...
    G.versions{C.v(k)}, G.choix_K(C.c(k)).nom, G.min_cellules(C.imc(k)), G.min_SCE(C.ims(k)));
end


function t = texte_resultat(S, k)
C = S.C;
par = cell(1, numel(S.animaux));
for a = 1:numel(S.animaux)
    if isnan(C.exces_animal(k, a))
        par{a} = sprintf('%s : aucune', S.animaux{a});
    else
        par{a} = sprintf('%s %+.0f (n = %d)', S.animaux{a}, 100 * C.exces_animal(k, a), C.nt_animal(k, a));
    end
end
t = sprintf(['%d biaisées / %d testables = %.1f %% (attendu %.1f %%) ; excès %+.1f points [%s] ; ' ...
    '%d après FDR'], C.n_sig(k), C.n_test(k), 100 * C.n_sig(k) / max(C.n_test(k), 1), ...
    100 * C.attendu(k) / max(C.n_test(k), 1), 100 * C.exces(k), strjoin(par, ', '), C.n_fdr(k));
end


function t = ligne_niveau(nom, Ff, q, meilleur)
t = sprintf('%-26s %+5.1f points  [hasard %+.1f à %+.1f]  meilleur %+5.1f  (%d combinaisons)', nom, ...
    100 * Ff.moyenne(q), 100 * Ff.nul_bas(q), 100 * Ff.nul_haut(q), 100 * Ff.meilleur(q), Ff.n(q));
if meilleur, t = [t '  <-']; end
end


function t = texte_p(p, B)
if isnan(p)
    t = '-';
elseif p <= 1 / (B + 1) + 1e-12
    t = sprintf('< %.3g', 1 / B);
else
    t = sprintf('%.3g', p);
end
end


function L = lignes_config(S, k)
% Lignes à copier après config = config_brainbow(); pour refaire cette
% combinaison avec le pipeline normal.
C = S.C; G = S.G;
L = {};
p = C.ip(k);
mo = G.pre(p).modif;
noms = fieldnames(mo);
% L'étape 1 change-t-elle par rapport à la base ? (config.pre ou config.asm)
base_asm = S.config_base.asm;
mm = G.membres(C.m(k)).modif;
nm = fieldnames(mm);
change_asm = ~strcmp(G.clusters{C.r(k), 1}, base_asm.seuil_clusters) || ...
    G.clusters{C.r(k), 2} ~= base_asm.percentile_seuil;
for z = 1:numel(nm)
    change_asm = change_asm || ~isfield(base_asm, nm{z}) || ~isequal(base_asm.(nm{z}), mm.(nm{z}));
end
if ~isempty(noms) || change_asm
    L{end+1} = ['config.PathSave = ''...'';     % un NOUVEAU dossier d''étape 1 : config.pre ou config.asm ' ...
        'changent, l''étape 1 est refaite (ne pas écraser l''ancien)'];
    L{end+1} = 'config.analyse = '''';';
else
    L{end+1} = '% étape 1 inchangée : même config.PathSave ; seules les étapes 2 à 4 sont à refaire';
    L{end+1} = 'config.analyse = ''screening'';     % nom au choix';
end
if isempty(noms)
    L{end+1} = '% prétraitement : celui de la base (config.pre inchangé)';
end
for z = 1:numel(noms)
    L{end+1} = sprintf('config.pre.%s = %s;', noms{z}, texte_valeur(mo.(noms{z}))); %#ok<AGROW>
end
L{end+1} = sprintf('config.asm.seuil_clusters = ''%s'';', G.clusters{C.r(k), 1});
L{end+1} = sprintf('config.asm.percentile_seuil = %g;', G.clusters{C.r(k), 2});
mo = G.membres(C.m(k)).modif;
noms = fieldnames(mo);
for z = 1:numel(noms)
    L{end+1} = sprintf('config.asm.%s = %s;', noms{z}, texte_valeur(mo.(noms{z}))); %#ok<AGROW>
end
L{end+1} = sprintf('config.version_assemblees = ''%s'';', G.versions{C.v(k)});
ch = G.choix_K(C.c(k));
switch ch.type
    case 'fixe'
        L{end+1} = sprintf('config.liste_K = %s;', mat2str(G.liste_K));
        L{end+1} = 'config.choixK.regle = ''exces'';';
        L{end+1} = sprintf('config.choixK.K_fixe = %d;', ch.K);
        L{end+1} = sprintf('config.choixK.k_start = %d;', min(G.liste_K));
    case 'exces'
        L{end+1} = sprintf('config.liste_K = %s;', mat2str(G.liste_K));
        L{end+1} = 'config.choixK.regle = ''exces'';';
        L{end+1} = 'config.choixK.K_fixe = [];';
        L{end+1} = sprintf('config.choixK.k_start = %d;', min(G.liste_K));
    case 'consensus'
        L{end+1} = sprintf('config.liste_K = %s;', mat2str(ch.liste_K));
        L{end+1} = 'config.choixK.regle = ''consensus'';';
        L{end+1} = 'config.choixK.K_fixe = [];';
        L{end+1} = sprintf('config.choixK.k_start = %d;', min(ch.liste_K));
        L{end+1} = sprintf('config.choixK.consensus.seuil = %g;', ch.seuil);
        L{end+1} = sprintf('config.choixK.consensus.min_K = %g;', ch.min_K);
        L{end+1} = sprintf('config.choixK.consensus.normalisation = ''%s'';', ch.normalisation);
end
L{end+1} = sprintf('config.agreg.min_cellules_asm = %g;', G.min_cellules(C.imc(k)));
L{end+1} = sprintf('config.agreg.min_SCE_asm = %g;', G.min_SCE(C.ims(k)));
end


function t = texte_valeur(v)
if ischar(v)
    t = ['''' v ''''];
elseif islogical(v)
    if v, t = 'true'; else, t = 'false'; end
elseif isscalar(v)
    t = sprintf('%g', v);
else
    t = mat2str(v);
end
end


function t = texte_struct(s)
noms = fieldnames(s);
parts = cell(1, numel(noms));
for z = 1:numel(noms)
    parts{z} = sprintf('''%s'', %s', noms{z}, texte_valeur(s.(noms{z})));
end
t = ['struct(' strjoin(parts, ', ') ')'];
end


function t = ternaire_r(c, a, b)
if c, t = a; else, t = b; end
end
