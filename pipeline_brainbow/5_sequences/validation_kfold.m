function V = validation_kfold(DF, debuts, fins, opts, graine)
%% =========================================================================
%  VALIDATION_KFOLD (19 nov. 2026, demandée par l'utilisateur)
%  Reproductibilité de la séquence de course, validation croisée répétée :
%    - n_courses >= opts.kfold_min_courses (10) : opts.kfold_k paquets (5),
%      opts.kfold_repetitions mélanges (20) -> 100 coupures 80 / 20 ;
%    - sinon : 2 paquets (50 / 50), 100 / 2 = 50 mélanges -> 100 coupures.
%    Chaque course est testée exactement une fois par mélange.
%  Pour chaque coupure : cellules de la séquence et décalage sur les courses
%  d'apprentissage (decalage_cellules), décalage sur les courses de test ;
%  rho de Spearman (apprentissage contre test) sur les cellules de la
%  séquence d'apprentissage (au moins 5).
%  Session :
%    rho   = moyenne des rho des coupures ;
%    p_rho = permutation : dans chaque coupure, décalages de test mélangés
%            entre les cellules de la séquence, moyenne des rho refaite
%            (opts.kfold_nPerm fois) ; unilatéral (rho > 0).
%  Cellules validées : retenues sur l'apprentissage ET dans la séquence du
%  test, dans au moins opts.kfold_seuil_valide (50 %) des coupures. Leur
%  ordre = décalage moyen mesuré sur les courses de test (coupures où elles
%  sont retenues sur l'apprentissage) ; rangées comme decalage_cellules (du
%  plus proche de 0 au plus tardif).
%  =========================================================================
n = numel(debuts);
NCell = size(DF, 1);
if n >= opts.kfold_min_courses
    K = opts.kfold_k; nrep = opts.kfold_repetitions;
else
    K = 2; nrep = round(opts.kfold_k * opts.kfold_repetitions / 2);
end
ncoup = K * nrep;
rho = nan(ncoup, 1);
rangs = cell(ncoup, 2);                      % rangs (apprentissage, test) des cellules de la séquence
n_app = zeros(NCell, 1); n_val = zeros(NCell, 1);
somme_dec_test = zeros(NCell, 1);
n_sel = nan(ncoup, 1); n_deux = nan(ncoup, 1); n_deux_hasard = nan(ncoup, 1);
fixer_flux_aleatoire(graine, 7, 1);
c = 0;
for r = 1:nrep
    melange = randperm(n);
    paquet = mod(0:n - 1, K) + 1;
    for k = 1:K
        c = c + 1;
        test = sort(melange(paquet == k));
        app = sort(melange(paquet ~= k));
        [~, d_app, ~, s_app] = decalage_cellules(DF, debuts, fins, app, opts.max_lag);
        [~, d_test, ~, s_test] = decalage_cellules(DF, debuts, fins, test, opts.max_lag);
        s_app = s_app(:); s_test = s_test(:);
        sel = s_app & isfinite(d_test(:));        % décalage de test mesurable
        n_app = n_app + sel;
        n_val = n_val + (s_app & s_test);
        somme_dec_test(sel) = somme_dec_test(sel) + d_test(sel)';
        n_sel(c) = sum(s_app);
        n_deux(c) = sum(s_app & s_test);
        n_deux_hasard(c) = sum(s_app) * sum(s_test) / NCell;
        if sum(sel) >= 5
            ra = tiedrank(d_app(sel)'); rt = tiedrank(d_test(sel)');
            if std(ra) > 0 && std(rt) > 0
                rho(c) = corr(ra, rt);
                rangs(c, :) = {ra, rt};
            end
        end
    end
end
V = struct();
V.validation = sprintf('%d-fold x %d', K, nrep);
V.K = K; V.n_repetitions = nrep; V.n_coupures = ncoup;
V.rho_coupures = rho;
V.rho = mean(rho, 'omitnan');
V.n_coupures_rho = sum(isfinite(rho));
V.n_sequence = mean(n_sel);
V.n_deux = mean(n_deux);
V.n_deux_hasard = mean(n_deux_hasard);

% permutation : décalages de test mélangés entre les cellules de la séquence
V.p_rho = NaN;
ok = find(isfinite(rho));
if ~isempty(ok)
    nulle = zeros(opts.kfold_nPerm, 1);
    for p = 1:opts.kfold_nPerm
        s = 0;
        for q = ok'
            ra = rangs{q, 1}; rt = rangs{q, 2};
            s = s + corr(ra, rt(randperm(numel(rt))));
        end
        nulle(p) = s / numel(ok);
    end
    V.p_rho = (1 + sum(nulle >= V.rho - 1e-12)) / (opts.kfold_nPerm + 1);
end

% cellules validées et leur ordre (décalage moyen sur les courses de test)
V.frac_validee = n_val / ncoup;
V.validees = V.frac_validee >= opts.kfold_seuil_valide;
V.decalage_test_moyen = nan(NCell, 1);
V.decalage_test_moyen(n_app > 0) = somme_dec_test(n_app > 0) ./ n_app(n_app > 0);
idx = find(V.validees);
[~, o] = sort(V.decalage_test_moyen(idx), 'descend');
V.ordre_valide = idx(o)';
V.n_validees = numel(idx);
end
