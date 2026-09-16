%% ========================================================================
% ANALYSE COMPARATIVE : PROÉMINENCE ET LARGEUR DES TRANSIENTS
% (à partir des sorties déjà calculées par transient_detection_shot_theta)
% ========================================================================
[NCell, Nz] = size(Raster);
sync_win = opts.synchronous_frames;
sampling_rate_local = sampling_rate; % suppose déjà chargée dans le workspace

% =========================================================================
% 1. RECONSTRUCTION DIRECTE DES MATRICES (aucun findpeaks relancé)
% =========================================================================
Mat_Prominence = zeros(NCell, Nz);
Mat_Width      = zeros(NCell, Nz);

for c = 1:NCell
    if ~isempty(Acttmp2{c})
        Mat_Prominence(c, Acttmp2{c}) = PromALL{c};
        Mat_Width(c, Acttmp2{c})      = WidthALL{c};
    end
end

% =========================================================================
% 2. MASQUES TEMPORELS VECTORISÉS
% =========================================================================
sce_asm_idx     = [Cl1, Cl2, Cl3];
sce_non_asm_idx = Cl0;
mask_sce_asm     = false(1, Nz);
mask_sce_non_asm = false(1, Nz);

t_start_asm = TRace(sce_asm_idx);
t_end_asm   = min(Nz, t_start_asm + sync_win);
for s = 1:length(sce_asm_idx)
    mask_sce_asm(t_start_asm(s):t_end_asm(s)) = true;
end

t_start_non = TRace(sce_non_asm_idx);
t_end_non   = min(Nz, t_start_non + sync_win);
for s = 1:length(sce_non_asm_idx)
    mask_sce_non_asm(t_start_non(s):t_end_non(s)) = true;
end

mask_sce_non_asm(mask_sce_asm) = false;

% =========================================================================
% 3. MASQUE "CORE"
% =========================================================================
CellCl_mat = zeros(NCl, NCell);
for k = 1:NCl
    CellCl_mat(k, cell2mat(assembly(k))) = 1;
end

is_core_mask = false(NCell, Nz);
for s = 1:length(sce_asm_idx)
    orig_sce_id = sce_asm_idx(s);
    active_asms = find(PCl_final(:, orig_sce_id) > 0);
    if isempty(active_asms)
        continue;
    end
    cells_core = any(CellCl_mat(active_asms, :), 1);
    is_core_mask(cells_core, t_start_asm(s):t_end_asm(s)) = true;
end

% =========================================================================
% 4. EXTRACTION PAR CATÉGORIE : PROÉMINENCE + LARGEUR
% =========================================================================
prom_isolated = []; prom_sce_non_asm = []; prom_sce_asm_all = []; prom_sce_asm_core = [];
width_isolated = []; width_sce_non_asm = []; width_sce_asm_all = []; width_sce_asm_core = [];

med_isolated_prom_percell  = nan(NCell, 1);
med_core_prom_percell      = nan(NCell, 1);
med_isolated_width_percell = nan(NCell, 1);
med_core_width_percell     = nan(NCell, 1);

for c = 1:NCell
    row_p = Mat_Prominence(c, :);
    row_w = Mat_Width(c, :);
    active_frames = find(row_p > 0);
    if isempty(active_frames)
        continue;
    end
    
    vals_p     = row_p(active_frames);
    vals_w_sec = row_w(active_frames) / sampling_rate_local;
    
    is_asm     = mask_sce_asm(active_frames);
    is_non_asm = mask_sce_non_asm(active_frames);
    is_iso     = ~is_asm & ~is_non_asm;
    is_core    = is_asm & is_core_mask(c, active_frames);
    
    prom_isolated     = [prom_isolated,     vals_p(is_iso)];     %#ok<AGROW>
    prom_sce_non_asm  = [prom_sce_non_asm,  vals_p(is_non_asm)]; %#ok<AGROW>
    prom_sce_asm_all  = [prom_sce_asm_all,  vals_p(is_asm)];     %#ok<AGROW>
    prom_sce_asm_core = [prom_sce_asm_core, vals_p(is_core)];    %#ok<AGROW>
    
    width_isolated     = [width_isolated,     vals_w_sec(is_iso)];     %#ok<AGROW>
    width_sce_non_asm  = [width_sce_non_asm,  vals_w_sec(is_non_asm)]; %#ok<AGROW>
    width_sce_asm_all  = [width_sce_asm_all,  vals_w_sec(is_asm)];     %#ok<AGROW>
    width_sce_asm_core = [width_sce_asm_core, vals_w_sec(is_core)];    %#ok<AGROW>
    
    if any(is_iso)
        med_isolated_prom_percell(c)  = median(vals_p(is_iso));
        med_isolated_width_percell(c) = median(vals_w_sec(is_iso));
    end
    if any(is_core)
        med_core_prom_percell(c)  = median(vals_p(is_core));
        med_core_width_percell(c) = median(vals_w_sec(is_core));
    end
end

% =========================================================================
% 5. STATISTIQUES DESCRIPTIVES
% =========================================================================
fprintf('\n================== BILAN PROÉMINENCE (dF/F0) ==================\n');
fprintf('1. Isolés       : N=%5d | Médiane=%.4f | Moyenne=%.4f\n', length(prom_isolated), median(prom_isolated), mean(prom_isolated));
fprintf('2. SCE Non-Asm  : N=%5d | Médiane=%.4f | Moyenne=%.4f\n', length(prom_sce_non_asm), median(prom_sce_non_asm), mean(prom_sce_non_asm));
fprintf('3. SCE Asm All  : N=%5d | Médiane=%.4f | Moyenne=%.4f\n', length(prom_sce_asm_all), median(prom_sce_asm_all), mean(prom_sce_asm_all));
fprintf('4. SCE Asm Core : N=%5d | Médiane=%.4f | Moyenne=%.4f\n', length(prom_sce_asm_core), median(prom_sce_asm_core), mean(prom_sce_asm_core));

paired_p = ~isnan(med_isolated_prom_percell) & ~isnan(med_core_prom_percell);
if sum(paired_p) >= 5
    [p_wilcoxon_prom, ~] = signrank(med_isolated_prom_percell(paired_p), med_core_prom_percell(paired_p));
    fprintf('Wilcoxon appairé Proéminence (N=%d) : p = %.4f\n', sum(paired_p), p_wilcoxon_prom);
else
    p_wilcoxon_prom = NaN;
end

fprintf('\n================== BILAN LARGEUR (secondes) ==================\n');
fprintf('1. Isolés       : N=%5d | Médiane=%.3f s | Moyenne=%.3f s\n', length(width_isolated), median(width_isolated), mean(width_isolated));
fprintf('2. SCE Non-Asm  : N=%5d | Médiane=%.3f s | Moyenne=%.3f s\n', length(width_sce_non_asm), median(width_sce_non_asm), mean(width_sce_non_asm));
fprintf('3. SCE Asm All  : N=%5d | Médiane=%.3f s | Moyenne=%.3f s\n', length(width_sce_asm_all), median(width_sce_asm_all), mean(width_sce_asm_all));
fprintf('4. SCE Asm Core : N=%5d | Médiane=%.3f s | Moyenne=%.3f s\n', length(width_sce_asm_core), median(width_sce_asm_core), mean(width_sce_asm_core));

paired_w = ~isnan(med_isolated_width_percell) & ~isnan(med_core_width_percell);
if sum(paired_w) >= 5
    [p_wilcoxon_width, ~] = signrank(med_isolated_width_percell(paired_w), med_core_width_percell(paired_w));
    fprintf('Wilcoxon appairé Largeur (N=%d) : p = %.4f\n', sum(paired_w), p_wilcoxon_width);
else
    p_wilcoxon_width = NaN;
end
fprintf('===============================================================\n\n');

%% ========================================================================
% EXPORT EXCEL
% ========================================================================
Categorie = {'1. Isolés (Hors SCE)'; '2. SCE Hors Assemblée (Cl0)'; ...
    '3. SCE Assemblée (Tous)'; '4. SCE Assemblée (Core)'};

T_stats_prom = table(Categorie, ...
    [length(prom_isolated); length(prom_sce_non_asm); length(prom_sce_asm_all); length(prom_sce_asm_core)], ...
    [median(prom_isolated); median(prom_sce_non_asm); median(prom_sce_asm_all); median(prom_sce_asm_core)], ...
    [mean(prom_isolated); mean(prom_sce_non_asm); mean(prom_sce_asm_all); mean(prom_sce_asm_core)], ...
    'VariableNames', {'Categorie', 'N', 'Mediane_prom', 'Moyenne_prom'});

T_stats_width = table(Categorie, ...
    [length(width_isolated); length(width_sce_non_asm); length(width_sce_asm_all); length(width_sce_asm_core)], ...
    [median(width_isolated); median(width_sce_non_asm); median(width_sce_asm_all); median(width_sce_asm_core)], ...
    [mean(width_isolated); mean(width_sce_non_asm); mean(width_sce_asm_all); mean(width_sce_asm_core)], ...
    'VariableNames', {'Categorie', 'N', 'Mediane_sec', 'Moyenne_sec'});

Parametre = {'Test Proéminence'; 'p-value'; 'Test Largeur'; 'p-value'};
Valeur = {'Wilcoxon Isolé vs Core'; num2str(p_wilcoxon_prom, '%.4f'); ...
          'Wilcoxon Isolé vs Core'; num2str(p_wilcoxon_width, '%.4f')};
T_test = table(Parametre, Valeur);

CellID = (1:NCell)';
T_percell = table(CellID, med_isolated_prom_percell, med_core_prom_percell, ...
    med_isolated_width_percell, med_core_width_percell, ...
    'VariableNames', {'CellID', 'Mediane_Prom_Isolee', 'Mediane_Prom_Core', ...
                      'Mediane_Largeur_Isolee_sec', 'Mediane_Largeur_Core_sec'});

if exist('namefull', 'var') && isfolder(namefull)
    excel_path = fullfile(namefull, 'bilan_prominence_largeur_transients.xlsx');
else
    excel_path = 'bilan_prominence_largeur_transients.xlsx';
end

writetable(T_stats_prom,  excel_path, 'Sheet', 'Statistiques_Prominence');
writetable(T_stats_width, excel_path, 'Sheet', 'Statistiques_Largeur');
writetable(T_test,        excel_path, 'Sheet', 'Tests_Wilcoxon');
writetable(T_percell,     excel_path, 'Sheet', 'Detail_Par_Cellule');
fprintf('Fichier Excel enregistré : %s\n', excel_path);

%% ========================================================================
% FIGURE 1 : PROÉMINENCE
% ========================================================================
fig_prom = figure('Color', 'w', 'Position', [100, 100, 1100, 450], 'Visible', 'off');
subplot(1, 2, 1); hold on;
[f_iso, x_iso] = ecdf(prom_isolated);
[f_non, x_non] = ecdf(prom_sce_non_asm);
[f_asm, x_asm] = ecdf(prom_sce_asm_all);
[f_cor, x_cor] = ecdf(prom_sce_asm_core);
plot(x_iso, f_iso, 'Color', [0.6 0.6 0.6], 'LineWidth', 1.8, 'DisplayName', 'Isolés');
plot(x_non, f_non, 'Color', [0.2 0.6 0.8], 'LineWidth', 1.8, 'DisplayName', 'SCE Non-Asm');
plot(x_asm, f_asm, 'Color', [0.9 0.4 0.1], 'LineWidth', 2.0, 'DisplayName', 'SCE Asm (All)');
plot(x_cor, f_cor, 'Color', [0.8 0.1 0.1], 'LineWidth', 2.2, 'DisplayName', 'SCE Asm (Core)');
xlabel('Proéminence \DeltaF/F_0'); ylabel('Probabilité cumulée');
title('Distribution cumulée des proéminences'); legend('Location', 'southeast'); grid on; box off;

subplot(1, 2, 2);
data_box = [prom_isolated(:); prom_sce_non_asm(:); prom_sce_asm_all(:); prom_sce_asm_core(:)];
group_box = [repmat({'1. Isolés'}, length(prom_isolated), 1); repmat({'2. SCE Non-Asm'}, length(prom_sce_non_asm), 1); ...
             repmat({'3. SCE Asm (All)'}, length(prom_sce_asm_all), 1); repmat({'4. SCE Asm (Core)'}, length(prom_sce_asm_core), 1)];
boxplot(data_box, group_box, 'Symbol', 'r.', 'OutlierSize', 3);
ylabel('Proéminence \DeltaF/F_0'); title('Comparaison des proéminences'); grid on; box off;

if exist('namefull', 'var') && isfolder(namefull)
    exportgraphics(fig_prom, fullfile(namefull, 'comparaison_prominence_transients.png'), 'Resolution', 300, 'BackgroundColor', 'white');
end

%% ========================================================================
% FIGURE 2 : LARGEUR
% ========================================================================
fig_width = figure('Color', 'w', 'Position', [100, 100, 1100, 450], 'Visible', 'off');
subplot(1, 2, 1); hold on;
[fw_iso, xw_iso] = ecdf(width_isolated);
[fw_non, xw_non] = ecdf(width_sce_non_asm);
[fw_asm, xw_asm] = ecdf(width_sce_asm_all);
[fw_cor, xw_cor] = ecdf(width_sce_asm_core);
plot(xw_iso, fw_iso, 'Color', [0.6 0.6 0.6], 'LineWidth', 1.8, 'DisplayName', 'Isolés');
plot(xw_non, fw_non, 'Color', [0.2 0.6 0.8], 'LineWidth', 1.8, 'DisplayName', 'SCE Non-Asm');
plot(xw_asm, fw_asm, 'Color', [0.9 0.4 0.1], 'LineWidth', 2.0, 'DisplayName', 'SCE Asm (All)');
plot(xw_cor, fw_cor, 'Color', [0.8 0.1 0.1], 'LineWidth', 2.2, 'DisplayName', 'SCE Asm (Core)');
xlabel('Largeur du pic (secondes)'); ylabel('Probabilité cumulée');
title('Distribution cumulée des largeurs'); legend('Location', 'southeast'); grid on; box off;

subplot(1, 2, 2);
data_box_w = [width_isolated(:); width_sce_non_asm(:); width_sce_asm_all(:); width_sce_asm_core(:)];
group_box_w = [repmat({'1. Isolés'}, length(width_isolated), 1); repmat({'2. SCE Non-Asm'}, length(width_sce_non_asm), 1); ...
               repmat({'3. SCE Asm (All)'}, length(width_sce_asm_all), 1); repmat({'4. SCE Asm (Core)'}, length(width_sce_asm_core), 1)];
boxplot(data_box_w, group_box_w, 'Symbol', 'r.', 'OutlierSize', 3);
ylabel('Largeur du pic (secondes)'); title('Comparaison des largeurs'); grid on; box off;

if exist('namefull', 'var') && isfolder(namefull)
    exportgraphics(fig_width, fullfile(namefull, 'comparaison_largeur_transients.png'), 'Resolution', 300, 'BackgroundColor', 'white');
end