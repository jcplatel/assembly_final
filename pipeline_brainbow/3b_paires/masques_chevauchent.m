function chevauche = masques_chevauchent(Coomasqx, Coomasqy)
% Paires de ROI qui partagent au moins un pixel
N = numel(Coomasqx);
Lx = max(cellfun(@(x) max([x(:); 1]), Coomasqx));
Ly = max(cellfun(@(y) max([y(:); 1]), Coomasqy));
lin = cell(N, 1);
qui = cell(N, 1);
for m = 1:N
    x = round(Coomasqx{m}(:));  y = round(Coomasqy{m}(:));
    ok = x >= 1 & y >= 1;
    lin{m} = sub2ind([Ly Lx], y(ok), x(ok));
    qui{m} = repmat(m, numel(lin{m}), 1);
end
Pix = sparse(vertcat(lin{:}), vertcat(qui{:}), 1, Ly * Lx, N);
chevauche = full((Pix' * Pix) > 0);
end
