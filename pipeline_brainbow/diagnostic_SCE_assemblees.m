function [Ts, Ta] = diagnostic_SCE_assemblees(config, K, nNull)
%% =========================================================================
%  DIAGNOSTIC_SCE_ASSEMBLEES : les petits SCE, en nombre proche de ce que
%  donnerait le hasard, sont-ils quand même de vraies activations des
%  assemblées ? (7 nov., question de l'utilisateur : « si la participation
%  moyenne est de 30 %, pour avoir des assemblées de 30 il me faut des SCE
%  de 10 »)
%  -------------------------------------------------------------------------
%  NE RECALCULE ET NE MODIFIE RIEN. Pour chaque session : pretraitement.mat
%  (Raster, Race, repos, seuil des SCE) et les assemblées (assemblees.mat :
%  membres et activations retenues PCl_final, au K retenu par l'étape 2 ou
%  au K donné).
%  1. SCE du hasard : le train de chaque cellule est décalé au hasard (bloc
%     circulaire) à l'intérieur du repos ; même détection que l'étape 1
%     (cellules actives entre t et t + synchronous_frames, pics au-dessus du
%     seuil de la session, MinPeakDistancesce). nNull tirages.
%  2. Pour chaque assemblée : membres présents dans chacune de ses
%     activations retenues ; « activation minimale » = 10e percentile de ce
%     nombre. Activations que le hasard pourrait imiter = SCE du hasard qui
%     contiennent au moins autant de membres (moyenne des tirages). C'est un
%     maximum : sous ce modèle toute l'activité est du hasard.
%  3. Petits SCE = au plus 3 cellules au-dessus du seuil (10-12 cellules
%     pour un seuil de 9) : combien sont des activations d'une assemblée ?
%  Résultat sur 444175_221121_plane1, K = 10 : 241 activations, au plus 2
%  imitables par le hasard (1 %) ; les activations à 10-12 cellules
%  contiennent 3,5 à 6 membres (médiane), les SCE du hasard de même taille
%  0 à 1 (95e percentile 1 à 3) ; participation des membres 28-57 %.
%
%  Lancer :   diagnostic_SCE_assemblees                 % K retenu, toutes les sessions
%             diagnostic_SCE_assemblees([], 10)          % K = 10 partout
%             config = config_brainbow(); config.seulement = {'444175'};
%             diagnostic_SCE_assemblees(config, 10)
%  Écrit dans <PathSave>\diagnostic_SCE\ : diagnostic_SCE_sessions.csv,
%  diagnostic_SCE_assemblees.csv, resume_diagnostic_SCE.txt, diagnostic_SCE.png
%  (figure_diagnostic_SCE) ; le résumé est aussi affiché.
%  =========================================================================

if nargin < 1 || isempty(config), config = config_brainbow(); end
if nargin < 2, K = []; end
if nargin < 3 || isempty(nNull), nNull = 20; end
addpath(genpath(fileparts(mfilename('fullpath'))));
etat_rng = rng;
remettre_rng = onCleanup(@() rng(etat_rng)); %#ok<NASGU>
config = dossiers_analyse(config);
sessions = lire_table_sessions(config);
a_traiter = 1:height(sessions);
if isfield(config, 'seulement') && ~isempty(config.seulement)
    seulement = cellstr(config.seulement);
    a_traiter = find(ismember(sessions.nom, seulement) | ismember(sessions.id_animal, seulement))';
end
lignes_s = cell(0, 13);
la = cell(0, 12);
for k = a_traiter
    session = table2struct(sessions(k, :));
    try
        [A, K_s] = charger_assemblees(session, K);
        if isempty(A), continue; end
        P = load(fullfile(session.dossier_source, 'pretraitement.mat'), 'Raster', 'Race', 'WinRest', 'opts', ...
            'diagnostic');
        Race = P.Race > 0;
        [NCell, NRace] = size(Race);
        s = seuil_session(P);
        sync = P.opts.synchronous_frames;
        dmin = P.opts.MinPeakDistancesce;
        PCl = A.PCl_final > 0;
        if size(PCl, 2) ~= NRace && size(PCl, 1) == NRace, PCl = PCl'; end
        asm = A.assemblyortho;
        nA = numel(asm);
        taille = sum(Race, 1);
        petit = taille <= s + 3;

        % 1. SCE du hasard
        R = P.Raster(:, P.WinRest) > 0;
        Nr = size(R, 2);
        nul = cell(1, nNull);
        for n = 1:nNull
            fixer_flux_aleatoire(session.graine, 20 + n, 1);
            Rn = false(NCell, Nr);
            for c = 1:NCell
                Rn(c, :) = circshift(R(c, :), randi(Nr));
            end
            nul{n} = race_du_hasard(Rn, s, sync, dmin);
        end
        n_hasard = mean(cellfun(@(x) size(x, 2), nul));

        % 2. chaque assemblée
        imitables = zeros(1, nA);
        act_tous = false(1, NRace);
        for a = 1:nA
            m = unique(asm{a}(:))';
            m = m(m >= 1 & m <= NCell);
            act = find(PCl(a, :));
            act_tous(act) = true;
            if isempty(act) || isempty(m), continue; end
            pres = sum(Race(m, act), 1);
            seuil_m = prctile(pres, 10);
            pres_nul = cellfun(@(x) sum(x(m, :), 1), nul, 'UniformOutput', false);
            taille_nul = cellfun(@(x) sum(x, 1), nul, 'UniformOutput', false);
            imitables(a) = mean(cellfun(@(p) sum(p >= seuil_m), pres_nul));
            p_petits = cell2mat(cellfun(@(p, t) p(t <= s + 3), pres_nul, taille_nul, 'UniformOutput', false));
            if isempty(p_petits), p_petits = NaN; end
            la(end+1, :) = {session.id_animal, session.nom, K_s, a, numel(m), numel(act), ...
                median(pres) / numel(m), sum(petit(act)), median(pres(petit(act))), seuil_m, ...
                prctile(p_petits, 95), imitables(a)}; %#ok<AGROW>
        end

        % 3. la session
        lignes_s(end+1, :) = {session.id_animal, session.nom, K_s, s, NRace, n_hasard, sum(petit), ...
            sum(petit & act_tous), nA, nnz(PCl), ...
            sum(imitables), sum(~act_tous), sum(petit & ~act_tous)}; %#ok<AGROW>
        fprintf(['[%s] K = %d, seuil %d : %d SCE (%.0f attendus par hasard) ; %d petits, dont %d ' ...
            'activations d''une assemblée ; %d activations, au plus %.1f imitables par le hasard\n'], ...
            session.nom, K_s, s, NRace, n_hasard, sum(petit), sum(petit & act_tous), lignes_s{end, 10}, sum(imitables));
    catch err
        fprintf(2, '[%s] : ÉCHEC : %s\n', session.nom, err.message);
    end
end
Ts = cell2table(lignes_s, 'VariableNames', {'animal', 'session', 'K', 'seuil_SCE', 'n_SCE', 'n_SCE_hasard', ...
    'n_petits', 'n_petits_activations', 'n_assemblees', 'n_activations', 'n_activations_imitables', ...
    'n_SCE_hors_activations', 'n_petits_hors_activations'});
Ta = cell2table(la, 'VariableNames', {'animal', 'session', 'K', 'assemblee', 'n_membres', 'n_activations', ...
    'participation_mediane', 'n_activations_petits_SCE', 'membres_dans_petites_activations', ...
    'activation_minimale_membres', 'membres_petits_SCE_hasard_95e', 'activations_imitables'});
dossier = fullfile(config.PathSave, 'diagnostic_SCE');
if ~isfolder(dossier), mkdir(dossier); end
writetable(Ts, fullfile(dossier, 'diagnostic_SCE_sessions.csv'));
writetable(Ta, fullfile(dossier, 'diagnostic_SCE_assemblees.csv'));
if isempty(Ts), return; end
L = {};
L{end+1} = sprintf('PETITS SCE ET ASSEMBLÉES (diagnostic_SCE_assemblees, %d sessions) :', height(Ts));
L{end+1} = sprintf(['   %d SCE, dont %d petits (au plus 3 cellules au-dessus du seuil) ; %d petits ' ...
    'sont des activations d''une assemblée (%.0f %% des petits).'], sum(Ts.n_SCE), ...
    sum(Ts.n_petits), sum(Ts.n_petits_activations), 100 * sum(Ts.n_petits_activations) / max(1, sum(Ts.n_petits)));
L{end+1} = sprintf(['   %d activations d''assemblées : %.0f %% dans de petits SCE ; au plus %.0f imitables par le ' ...
    'hasard (%.1f %%).'], sum(Ts.n_activations), 100 * sum(Ta.n_activations_petits_SCE) / ...
    max(1, sum(Ts.n_activations)), sum(Ts.n_activations_imitables), 100 * sum(Ts.n_activations_imitables) / ...
    max(1, sum(Ts.n_activations)));
L{end+1} = sprintf(['   Participation des membres à leurs activations : médiane %.0f %% ; membres dans les activations ' ...
    'à petits SCE : médiane %.1f, contre %.1f au 95e percentile des petits SCE du hasard.'], ...
    100 * median(Ta.participation_mediane, 'omitnan'), median(Ta.membres_dans_petites_activations, 'omitnan'), ...
    median(Ta.membres_petits_SCE_hasard_95e, 'omitnan'));
ratio = Ts.n_SCE_hasard ./ Ts.n_SCE;
L{end+1} = sprintf(['   Sessions avec plus de SCE attendus par hasard qu''observés : %d sur %d ; activations ' ...
    'imitables dans ces sessions : %.1f %%, dans les autres : %.1f %%.'], sum(ratio >= 1), height(Ts), ...
    100 * sum(Ts.n_activations_imitables(ratio >= 1)) / max(1, sum(Ts.n_activations(ratio >= 1))), ...
    100 * sum(Ts.n_activations_imitables(ratio < 1)) / max(1, sum(Ts.n_activations(ratio < 1))));
fid = fopen(fullfile(dossier, 'resume_diagnostic_SCE.txt'), 'w', 'n', 'UTF-8');
if fid >= 0
    fprintf(fid, '%s\n', L{:});
    fclose(fid);
end
fprintf('\n%s\n', L{:});
try
    figure_diagnostic_SCE(dossier);
catch err
    warning('figure_diagnostic_SCE non exécutée : %s', err.message);
end
end


function [A, K_s] = charger_assemblees(session, K)
% assemblees.mat au K donné, ou au K retenu par l'étape 2 (0 = consensus)
A = [];
K_s = K;
if isempty(K)
    f = fullfile(session.dossier, 'choix_K.mat');
    if ~isfile(f), return; end
    c = load(f);
    if ~isfield(c, 'statut_choix') || ~strcmp(c.statut_choix, 'ok'), return; end
    K_s = c.best_K;
end
fa = fichier_assemblees(session.dossier, session.dossier_source, K_s);   % v1 ou v2 selon l'analyse
if ~isfile(fa), return; end
A = load(fa, 'assemblyortho', 'PCl_final');
if isempty(A.assemblyortho), A = []; end
end


function s = seuil_session(P)
% seuil des SCE utilisé à l'étape 1
if isfield(P.opts, 'sce_n_cells_threshold_calcule')
    s = P.opts.sce_n_cells_threshold_calcule;
elseif isfield(P, 'diagnostic') && isfield(P.diagnostic, 'seuil_SCE')
    s = P.diagnostic.seuil_SCE;
else
    s = P.opts.sce_n_cells_threshold;
end
end


function Rc = race_du_hasard(R, s, sync, dmin)
% même détection que SCE_detection5 : cellules actives entre t et t + sync,
% pics strictement au-dessus de s, séparés d'au moins dmin trames
Nz = size(R, 2);
Rc = false(size(R, 1), 0);
if Nz <= sync + dmin + 1, return; end
F = movmax(double(R), [0 sync], 2);
MAct = sum(F(:, 1:Nz - sync), 1);
if max(MAct) <= s, return; end
[pks, locs] = findpeaks(MAct, 'MinPeakHeight', s, 'MinPeakDistance', dmin);
locs = locs(pks <= 150);
Rc = F(:, locs) > 0;
end
