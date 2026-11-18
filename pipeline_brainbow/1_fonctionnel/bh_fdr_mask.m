function sig_mask = bh_fdr_mask(pvals, q)
% Correction FDR de Benjamini-Hochberg (ta fonction de RACEORTHO_FDR).
%   pvals    : vecteur ligne de p-values (une famille de tests)
%   q        : taux de fausses découvertes toléré (ex. 0.05)
%   sig_mask : vecteur logique, true si significatif après correction
m = numel(pvals);
[sorted_p, ~] = sort(pvals);
ranks = 1:m;
bh_thresh = (ranks / m) * q;
below = sorted_p <= bh_thresh;
sig_mask = false(1, m);
if any(below)
    k_max = find(below, 1, 'last');   % plus grand rang tel que p(k) <= (k/m)*q
    crit_p = sorted_p(k_max);
    sig_mask = pvals <= crit_p;
end
end
