function [p_detection, p_F0, p_composition_membres, compte_membres, compte_non_membres] = ...
    test_biais_detection(colorcell, assembly, th_detection, F0_mediane, nPerm, graine)
% Biais de détection lié à la couleur.
%  p_detection : le bruit de détection (th_detection, mesuré au repos)
%                diffère-t-il entre classes de couleur ? (Kruskal-Wallis)
%                Une protéine qui fuit dans le canal GCaMP augmente F0 et le
%                bruit en ΔF/F : la cellule a moins de transients détectés.
%  p_F0        : même question pour F0 (médiane de Fzero).
%  p_composition_membres : les membres d'assemblées ont-ils la même
%                répartition de couleurs que les non-membres ? (chi2,
%                permutation). ATTENTION : un vrai effet clonal sur quelques
%                couleurs rend aussi ce test significatif (vérifié sur
%                données simulées) ; il ne prouve donc pas à lui seul un biais.

nC = 6;
colorcell = colorcell(:);
est_coloree = colorcell >= 1 & colorcell <= 6;
p_detection = kruskal_couleurs(th_detection(:), colorcell, est_coloree);
p_F0        = kruskal_couleurs(F0_mediane(:), colorcell, est_coloree);

membre = false(numel(colorcell), 1);
membre([assembly{:}]) = true;
idx_col = find(est_coloree);
col_c  = colorcell(idx_col);
memb_c = membre(idx_col);
T_obs = tableau_couleurs(col_c, memb_c, nC);
compte_membres     = T_obs(:, 1)';
compte_non_membres = T_obs(:, 2)';
X_obs = chi2_tableau(T_obs);
rng(graine);
X_nul = zeros(1, nPerm);
for b = 1:nPerm
    X_nul(b) = chi2_tableau(tableau_couleurs(col_c, memb_c(randperm(numel(memb_c))), nC));
end
p_composition_membres = (1 + sum(X_nul >= X_obs * (1 - 1e-10))) / (nPerm + 1);
end


function p = kruskal_couleurs(valeurs, colorcell, est_coloree)
% Kruskal-Wallis entre les classes de couleur d'au moins 3 cellules
ok = est_coloree & ~isnan(valeurs);
classes = unique(colorcell(ok));
nb = arrayfun(@(c) sum(colorcell(ok) == c), classes);
garde = ok & ismember(colorcell, classes(nb >= 3));
if numel(unique(colorcell(garde))) < 2
    p = NaN;
else
    p = kruskalwallis(valeurs(garde), colorcell(garde), 'off');
end
end


function T = tableau_couleurs(col_c, memb_c, nC)
T = [accumarray(col_c(memb_c), 1, [nC 1]), accumarray(col_c(~memb_c), 1, [nC 1])];
end


function X = chi2_tableau(T)
attendu = sum(T, 2) * sum(T, 1) / sum(T(:));
X = sum((T - attendu).^2 ./ attendu, 'all', 'omitnan');
end
