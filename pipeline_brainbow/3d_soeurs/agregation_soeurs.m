function T = agregation_soeurs(sessions, config, K)
%AGREGATION_SOEURS  Étape 3d, toutes les sessions : les paires probablement
%   sœurs ont-elles une activité plus semblable ?
%   Par animal : moyenne des pentes beta de ses sessions (test_soeurs_session).
%   UN test par mesure : Wilcoxon unilatéral sur les animaux (beta > 0).
%   Effet en % = beta / similarité brute moyenne des paires à <= 100 µm :
%   différence attendue entre une paire de vraies sœurs et une paire non
%   sœur.
%   Dose-réponse : résidu moyen par catégorie de pi, par animal (moyenne
%   des sessions pondérée par le nombre de paires), toutes couleurs et
%   couleur par couleur.
%   Par couleur (descriptif, 6 tests sans correction) : beta_c moyen par
%   animal, Wilcoxon unilatéral sur les animaux (table Bc).
%   Écrit dans <PathSave>\soeurs_agregation\ : soeurs_sessions.csv,
%   soeurs_agregation.mat, resume_soeurs.txt, soeurs.png.

opts = options_soeurs(config);
dossier = fullfile(config.PathSave, 'soeurs_agregation');
if ~isfolder(dossier), mkdir(dossier); end
morceaux = {};
morceaux_c = {};
doses = {};
for s = 1:height(sessions)
    f = fullfile(sessions.dossier{s}, 'soeurs', 'soeurs.mat');
    if ~isfile(f), continue; end
    vars = who('-file', f);
    if ~ismember('T', vars), continue; end
    R = load(f, 'T', 'dose');
    morceaux{end+1, 1} = R.T; %#ok<AGROW>
    if ismember('Tc', vars)
        Rc = load(f, 'Tc');
        morceaux_c{end+1, 1} = Rc.Tc; %#ok<AGROW>
    end
    for nom = fieldnames(R.dose)'
        dd = R.dose.(nom{1});
        if ~isfield(dd, 'moyenne_c'), dd.moyenne_c = nan(6, numel(dd.n)); dd.n_c = zeros(6, numel(dd.n)); end
        doses{end+1, 1} = table({sessions.id_animal{s}}, nom, {dd.moyenne}, {dd.n}, {dd.moyenne_c}, ...
            {dd.n_c}, 'VariableNames', {'animal', 'mesure', 'moyenne', 'n', 'moyenne_c', 'n_c'}); %#ok<AGROW>
    end
end
T = table();
if isempty(morceaux)
    fprintf('Sœurs : aucune session avec des résultats (lancer main_soeurs).\n');
    return
end
T = vertcat(morceaux{:});
writetable(T, fullfile(dossier, 'soeurs_sessions.csv'));
Dz = vertcat(doses{:});

% par animal, puis entre animaux
[G, mes, anim] = findgroups(T.mesure, T.animal);
A = table(mes, anim, splitapply(@(x) mean(x, 'omitnan'), T.beta, G), ...
    splitapply(@(x) mean(x, 'omitnan'), T.beta_simple, G), splitapply(@(x) mean(x, 'omitnan'), T.effet_pct, G), ...
    splitapply(@numel, T.beta, G), splitapply(@(p) sum(p < 0.05), T.p, G), ...
    'VariableNames', {'mesure', 'animal', 'beta', 'beta_simple', 'effet_pct', 'n_sessions', 'n_sig'});
mesures = unique(T.mesure, 'stable');
B = table();
for k = 1:numel(mesures)
    a = A(strcmp(A.mesure, mesures{k}) & isfinite(A.beta), :);
    p = NaN; p_simple = NaN;
    if height(a) >= 2
        p = signrank(a.beta, 0, 'tail', 'right');
        p_simple = signrank(a.beta_simple, 0, 'tail', 'right');
    end
    B = [B; table(mesures(k), height(a), mean(a.effet_pct), std(a.effet_pct) / sqrt(height(a)), ...
        mean(a.beta), std(a.beta) / sqrt(height(a)), sum(a.beta > 0), p, p_simple, sum(a.n_sig), ...
        sum(a.n_sessions), 'VariableNames', {'mesure', 'n_animaux', 'effet_pct', 'effet_sem', 'beta_moyen', ...
        'beta_sem', 'n_pos', 'p_wilcoxon', 'p_simple', 'n_sessions_sig', 'n_sessions'})]; %#ok<AGROW>
end

% dose-réponse par animal (moyenne pondérée par le nombre de paires)
Ad = table();
for k = 1:numel(mesures)
    for an = unique(Dz.animal)'
        m = strcmp(Dz.animal, an{1}) & strcmp(Dz.mesure, mesures{k});
        if ~any(m), continue; end
        mo = vertcat(Dz.moyenne{m});
        nn = vertcat(Dz.n{m});
        mo(isnan(mo)) = 0;
        moy = sum(mo .* nn, 1) ./ sum(nn, 1);
        mc = cat(3, Dz.moyenne_c{m});               % couleur x catégorie x session
        nc = cat(3, Dz.n_c{m});
        mc(isnan(mc)) = 0;
        moy_c = sum(mc .* nc, 3) ./ sum(nc, 3);
        Ad = [Ad; table(mesures(k), an, {moy}, {sum(nn, 1)}, {moy_c}, {sum(nc, 3)}, 'VariableNames', ...
            {'mesure', 'animal', 'moyenne', 'n', 'moyenne_c', 'n_c'})]; %#ok<AGROW>
    end
end

% par couleur (descriptif)
Bc = table();
if ~isempty(morceaux_c)
    Tc = vertcat(morceaux_c{:});
    writetable(Tc, fullfile(dossier, 'soeurs_sessions_couleurs.csv'));
    [G, mes, anim, coul] = findgroups(Tc.mesure, Tc.animal, Tc.couleur);
    Ac = table(mes, anim, coul, splitapply(@(x) mean(x, 'omitnan'), Tc.beta, G), ...
        splitapply(@(x) mean(x, 'omitnan'), Tc.effet_pct, G), splitapply(@sum, Tc.n_meme_pi, G), ...
        'VariableNames', {'mesure', 'animal', 'couleur', 'beta', 'effet_pct', 'n_meme_pi'});
    for k = 1:numel(mesures)
        for c = {'R', 'G', 'B', 'RG', 'RB', 'GB'}
            a = Ac(strcmp(Ac.mesure, mesures{k}) & strcmp(Ac.couleur, c{1}) & isfinite(Ac.beta), :);
            p = NaN;
            if height(a) >= 2, p = signrank(a.beta, 0, 'tail', 'right'); end
            Bc = [Bc; table(mesures(k), c, height(a), mean(a.effet_pct), ...
                std(a.effet_pct) / sqrt(max(1, height(a))), mean(a.beta), sum(a.beta > 0), p, ...
                sum(a.n_meme_pi), 'VariableNames', {'mesure', 'couleur', 'n_animaux', 'effet_pct', ...
                'effet_sem', 'beta_moyen', 'n_pos', 'p_wilcoxon', 'n_paires_proches'})]; %#ok<AGROW>
        end
    end
end

L = resume(B, K, opts, Bc);
fid = fopen(fullfile(dossier, 'resume_soeurs.txt'), 'w', 'n', 'UTF-8');
if fid >= 0
    fprintf(fid, '%s\n', L{:});
    fclose(fid);
end
fprintf('\n%s\n', L{:});
lignes_resume = L; %#ok<NASGU>
date_soeurs = datestr(now, 'yyyy-mm-dd HH:MM'); %#ok<NASGU>
save(fullfile(dossier, 'soeurs_agregation.mat'), 'T', 'A', 'B', 'Ad', 'Bc', 'opts', 'lignes_resume', 'date_soeurs');
try
    figure_soeurs(dossier);
catch err
    warning('figure_soeurs non exécutée : %s', err.message);
end
end


function L = resume(B, K, opts, Bc)
L = {};
L{end+1} = sprintf(['PROBABILITÉ D''ÊTRE SŒURS (étape 3d, couleurs et positions seulement, %d sessions) : ' ...
    'part de vraies sœurs parmi les paires de même couleur (toutes couleurs ; mise à 0 au-delà de %g µm ' ...
    'dans le test) :'], K.n_sessions, opts.pi_max_um);
pb = K.pi_b;
if isfield(K, 'pi_brut_b'), pb = K.pi_brut_b; end
t = arrayfun(@(i) sprintf('%s µm : %.0f %%', texte_tranche(opts.bornes_um(i), opts.bornes_um(i + 1)), ...
    100 * pb(i)), 1:numel(pb), 'UniformOutput', false);
L{end+1} = ['   ' strjoin(t, ' ; ')];
if isfield(K, 'angles') && ~isempty(K.angles)
    g = K.angles;
    L{end+1} = sprintf(['DIRECTION DES BANDES DE COULEUR (angles_bandes.csv, angles_bandes_<animal>.png ; ' ...
        'chaque couleur mesurée sur l''angle trouvé avec les autres couleurs) : paires de même couleur le ' ...
        'long / en travers des bandes (observées / attendues, %g-%g µm), médianes des %d sessions : %.2f / ' ...
        '%.2f ; %d sessions avec le long > en travers.'], opts.angle_um(1), opts.angle_um(2), height(g), ...
        median(g.meme_couleur_le_long), median(g.meme_couleur_en_travers), ...
        sum(g.meme_couleur_le_long > g.meme_couleur_en_travers));
end
L{end+1} = ['LES PAIRES PROBABLEMENT SŒURS ONT-ELLES UNE ACTIVITÉ PLUS SEMBLABLE ? (résidus ~ pi + ' ...
    'même couleur par couleur ; un test par mesure : Wilcoxon unilatéral sur les animaux) :'];
for k = 1:height(B)
    switch B.mesure{k}
        case 'sttc',        nm = 'STTC au repos';
        case 'coassemblee', nm = 'même assemblée';
        case 'sce',         nm = 'co-participation aux SCE (phi)';
        otherwise,          nm = B.mesure{k};
    end
    if isfinite(B.effet_pct(k))
        effet = sprintf('%s %% de plus qu''une paire non sœur (± %s, sem des animaux)', ...
            nb(B.effet_pct(k), '%+.0f'), nb(B.effet_sem(k), '%.0f'));
    else
        effet = sprintf('un phi plus élevé de %s qu''une paire non sœur (± %s, sem des animaux)', ...
            nb(B.beta_moyen(k), '%+.4f'), nb(B.beta_sem(k), '%.4f'));
    end
    L{end+1} = sprintf(['   « %s : une paire de vraies sœurs a %s, %d/%d animaux > 0, %s ; sans les ' ...
        'indicateurs de couleur : %s ; %d/%d sessions à p < 0,05. »'], nm, effet, B.n_pos(k), ...
        B.n_animaux(k), pv(B.p_wilcoxon(k)), pv(B.p_simple(k)), B.n_sessions_sig(k), B.n_sessions(k)); %#ok<AGROW>
    if ~isempty(Bc)
        b = Bc(strcmp(Bc.mesure, B.mesure{k}), :);
        t = arrayfun(@(r) sprintf('%s %s (%s, %d paires proches)', b.couleur{r}, ...
            effet_court(b.effet_pct(r), b.beta_moyen(r)), pv(b.p_wilcoxon(r)), b.n_paires_proches(r)), ...
            1:height(b), 'UniformOutput', false);
        L{end+1} = ['      par couleur (descriptif, sans correction) : ' strjoin(t, ' ; ')]; %#ok<AGROW>
    end
end
end


function t = texte_tranche(a, b)
if isinf(b), t = sprintf('> %g', a); else, t = sprintf('%g-%g', a, b); end
end


function t = nb(x, fmt)
if isnan(x), t = 'NA'; else, t = sprintf(fmt, x); end
end


function t = pv(p)
if isnan(p)
    t = 'p = NA';
elseif p < 0.001
    t = 'p < 0,001';
else
    t = ['p = ' strrep(sprintf('%.3f', p), '.', ',')];
end
end


function t = effet_court(pct, beta)
% en % si la mesure a une échelle positive, sinon phi brut (co-participation aux SCE)
if isfinite(pct), t = sprintf('%+.0f %%', pct); else, t = sprintf('%+.4f', beta); end
end
