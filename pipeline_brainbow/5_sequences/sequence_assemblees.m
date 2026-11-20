function A = sequence_assemblees(ordre, NCell, assemblees, nPerm, min_paires, graine)
%SEQUENCE_ASSEMBLEES  Les assemblées (étape 2, SCE au repos) suivent-elles
%   la séquence de course ? (20 nov. 2026, demandé par l'utilisateur)
%   ordre      : cellules de la séquence, dans leur ordre
%   assemblees : cell array, cellules de chaque assemblée (assemblyortho)
%   a. Enrichissement : chaque assemblée contient-elle plus de cellules de
%      la séquence que le hasard ? Test hypergéométrique unilatéral
%      (NCell cellules, n_seq dans la séquence, taille de l'assemblée), FDR
%      de Benjamini-Hochberg 5 % sur les assemblées de la session. Aussi,
%      pour toute la session : part des cellules d'assemblées dans la
%      séquence contre part de toutes les cellules.
%   b. Regroupement : les cellules d'une même assemblée sont-elles proches
%      dans la séquence ? Même test que la couleur
%      (distance_sequence_couleur) : position = rang / (n - 1) ; paires de
%      cellules de la séquence qui sont dans au moins une assemblée ;
%      « même assemblée » = au moins une assemblée en commun ;
%      delta = moyenne(même) - moyenne(différentes), < 0 = plus proches ;
%      p : appartenances aux assemblées mélangées entre ces cellules (nPerm
%      fois), unilatéral.

nA = numel(assemblees);
A = struct('n_assemblees', nA, 'n_sequence', numel(ordre), 'taille', nan(nA, 1), ...
    'n_dans_sequence', nan(nA, 1), 'p_enrichissement', nan(nA, 1), 'enrichie', false(nA, 1), ...
    'n_enrichies', 0, 'frac_seq_assemblees', NaN, 'frac_seq_toutes', NaN, ...
    'n_cellules_paires', 0, 'n_meme', 0, 'n_diff', 0, 'd_meme', NaN, 'd_diff', NaN, ...
    'delta', NaN, 'pct', NaN, 'p', NaN);
if nA == 0, return; end
dans_seq = false(NCell, 1);
dans_seq(ordre) = true;
n_seq = numel(ordre);
A.frac_seq_toutes = n_seq / NCell;

%% a. enrichissement
membre = false(NCell, nA);
for a = 1:nA
    c = assemblees{a}(:);
    c = c(c >= 1 & c <= NCell);
    membre(c, a) = true;
    A.taille(a) = numel(c);
    A.n_dans_sequence(a) = sum(dans_seq(c));
    A.p_enrichissement(a) = hygecdf(A.n_dans_sequence(a) - 1, NCell, n_seq, A.taille(a), 'upper');
end
A.enrichie = bh_fdr_mask(A.p_enrichissement(:)', 0.05)';
A.n_enrichies = sum(A.enrichie);
dans_asm = any(membre, 2);
A.frac_seq_assemblees = sum(dans_seq & dans_asm) / max(1, sum(dans_asm));

%% b. regroupement dans la séquence
if n_seq < 3, return; end
position = nan(NCell, 1);
position(ordre) = (0:n_seq - 1)' / (n_seq - 1);
cel = ordre(dans_asm(ordre));
cel = cel(:);
A.n_cellules_paires = numel(cel);
if numel(cel) < 3, return; end
M = membre(cel, :);
[i, j] = find(triu(true(numel(cel)), 1));
d_seq = abs(position(cel(i)) - position(cel(j)));
meme = any(M(i, :) & M(j, :), 2);
A.n_meme = sum(meme);
A.n_diff = sum(~meme);
if A.n_meme < min_paires || A.n_diff < min_paires, return; end
A.d_meme = mean(d_seq(meme));
A.d_diff = mean(d_seq(~meme));
A.delta = A.d_meme - A.d_diff;
A.pct = 100 * A.delta / A.d_diff;
fixer_flux_aleatoire(graine, 9, 1);
delta_perm = zeros(nPerm, 1);
for k = 1:nPerm
    M2 = M(randperm(numel(cel)), :);
    m2 = any(M2(i, :) & M2(j, :), 2);
    delta_perm(k) = mean(d_seq(m2)) - mean(d_seq(~m2));
end
A.p = (1 + sum(delta_perm <= A.delta + 1e-12)) / (nPerm + 1);
end
