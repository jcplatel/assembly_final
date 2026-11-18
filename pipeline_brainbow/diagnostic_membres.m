function T = diagnostic_membres(dossier_K, nRep)
%DIAGNOSTIC_MEMBRES  Combien de faux membres et de fausses assemblées chaque
%   réglage de RACEORTHO trouve-t-il par hasard ?
%   Sur un dossier k<K> déjà calculé, on garde les clusters de SCE (IDX2,
%   NClOK) et on mélange les cellules par curveball : chaque cellule garde
%   son nombre de SCE, chaque SCE sa taille, mais tout lien entre une
%   cellule et un cluster est détruit. Toute cellule déclarée membre dans
%   ces matrices mélangées est donc un faux membre.
%   Aucune couleur n'est lue.
%
%   diagnostic_membres('E:\...\_retrouver_ancien\nouveau\444119_220919_plane0_sans_indice\k9')
%   diagnostic_membres(dossier_K, 20)      % 20 mélanges (défaut : 10)
%
%   Réglages comparés (le reste comme dans assemblees.mat de ce dossier) :
%     global 10 %, pondéré       (nouveau)
%     global 20 %, pondéré
%     par cellule 5 %, pondéré
%     par cellule 5 %, uniforme  (ancien)
%   Pour chaque réglage :
%     donnees_membres : cellules significatives de chaque cluster, vraies données
%                       (avant le filtre de taille)
%     donnees_asm     : tailles des assemblées finales, vraies données
%     hasard_membres  : faux membres par cluster (moyenne des mélanges)
%     fausses_asm     : fausses assemblées par mélange (après TOUS les filtres)
%                       et leur taille moyenne
%   Lecture : par cellule, le nombre de faux membres ne dépend pas des autres
%   cellules : hasard_membres estime directement les faux membres des vraies
%   assemblées. En global, BH s'adapte (plus il y a de vrais membres, plus le
%   seuil est souple) : le mélange sous-estime les faux membres, mais leur
%   PROPORTION parmi les membres déclarés reste garantie <= FDR_Q.
%   Durée : quelques minutes par réglage pondéré (10 mélanges).

if nargin < 2 || isempty(nRep), nRep = 10; end
dossier_K = regexprep(dossier_K, '[\\/]+$', '');
dossier_session = fileparts(dossier_K);
P = load(fullfile(dossier_session, 'pretraitement.mat'), 'Race');
Race = double(P.Race > 0);
A = load(fullfile(dossier_K, 'assemblees.mat'), 'IDX2', 'NClOK', 'config_asm');
NClOK = A.NClOK;
if NClOK < 1
    error('Aucun cluster de SCE significatif dans %s.', dossier_K);
end
base = A.config_asm;
if isempty(gcp('nocreate')), parpool; end
etat_rng = rng;
remettre_rng = onCleanup(@() rng(etat_rng)); %#ok<NASGU>

noms    = {'global 10 % ponderé', 'global 20 % ponderé', 'par cellule 5 % ponderé', 'par cellule 5 % uniforme'};
famille = {'globale', 'globale', 'par_cellule', 'par_cellule'};
q       = [0.10 0.20 0.05 0.05];
nul     = {'pondere', 'pondere', 'pondere', 'uniforme'};

% Les mêmes matrices mélangées pour tous les réglages
melanges = cell(1, nRep);
for r = 1:nRep
    rng(1000 + r, 'twister');
    melanges{r} = double(curveball(Race > 0));
end
fprintf('%s : %d cellules, %d SCE, %d clusters significatifs, %d mélanges curveball\n', ...
    dossier_K, size(Race, 1), size(Race, 2), NClOK, nRep);

nS = numel(noms);
lignes = cell(nS, 6);
for s = 1:nS
    t0 = tic;
    cfg = base;
    cfg.fdr_famille   = famille{s};
    cfg.FDR_Q         = q(s);
    cfg.null_cellules = nul{s};

    Rv = RACEORTHO_FDR2(Race, A.IDX2, NClOK, cfg, 1);
    m_vrais = n_significatifs(Rv.pValues_cells, Race, A.IDX2, NClOK, cfg);

    faux = zeros(nRep, NClOK);
    n_fausses = zeros(nRep, 1);
    tailles_fausses = [];
    for r = 1:nRep
        R0 = RACEORTHO_FDR2(melanges{r}, A.IDX2, NClOK, cfg, 100 + r);
        faux(r, :) = n_significatifs(R0.pValues_cells, melanges{r}, A.IDX2, NClOK, cfg);
        n_fausses(r) = numel(R0.assemblyortho);
        tailles_fausses = [tailles_fausses, cellfun(@numel, R0.assemblyortho)]; %#ok<AGROW>
    end
    t_fausses = '-';
    if ~isempty(tailles_fausses), t_fausses = sprintf('%.1f', mean(tailles_fausses)); end
    lignes(s, :) = {noms{s}, texte(m_vrais, '%d'), texte(cellfun(@numel, Rv.assemblyortho), '%d'), ...
        texte(mean(faux, 1), '%.1f'), sprintf('%.2f (taille %s)', mean(n_fausses), t_fausses), ...
        round(toc(t0) / 60, 1)};
    fprintf('  %s : fait (%.1f min)\n', noms{s}, toc(t0) / 60);
end

T = cell2table(lignes, 'VariableNames', {'reglage', 'donnees_membres', 'donnees_asm', ...
    'hasard_membres', 'fausses_asm', 'duree_min'});
fprintf('\n');
disp(T(:, 1:5))
fprintf(['donnees_membres / hasard_membres : par cluster de SCE (rang 1, 2, ...), avant le filtre de taille.\n' ...
    'fausses_asm : fausses assemblées par mélange (après tous les filtres) et leur taille moyenne.\n' ...
    'Par cellule : hasard_membres ~ faux membres contenus dans les vraies assemblées.\n' ...
    'Global : proportion de faux membres garantie <= FDR_Q (le mélange la sous-estime).\n']);
end


function n = n_significatifs(pv, Race, IDX2, NClOK, cfg)
% Cellules significatives de chaque cluster (même règle que RACEORTHO_FDR2,
% avant le filtre de taille)
n = zeros(1, NClOK);
if isempty(pv), return; end
Nact = sum(Race, 2);
CellP = zeros(size(Race, 1), NClOK);
for i = 1:NClOK
    CellP(:, i) = sum(Race(:, IDX2 == i), 2);
end
fam = Nact > 0;
sig = false(size(pv));
switch cfg.fdr_famille
    case 'globale'
        p = pv(fam, :);
        m = bh_fdr_mask(p(:)', cfg.FDR_Q);
        sig(fam, :) = reshape(m, size(p));
    case 'par_cellule'
        for j = find(fam)'
            sig(j, :) = bh_fdr_mask(pv(j, :), cfg.FDR_Q);
        end
end
sig = sig & (CellP >= cfg.min_participations);
n = sum(sig, 1);
end


function t = texte(v, fmt)
t = strjoin(arrayfun(@(x) sprintf(fmt, x), v, 'UniformOutput', false), ' ');
end
