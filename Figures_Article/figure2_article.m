% =========================================================================
% SCRIPT : FIGURES ARTICLE - CARACTÉRISATION DES ASSEMBLÉES + LMM
% Lecture depuis l'Excel des meilleurs K et extraction des .mat
% Exportation standardisée : 600 DPI, SVG, Police Arial
% Modèle LMM avec facteurs emboîtés (Session dans Animal)
% =========================================================================

clear; clc; close all;

%% 1. PARAMÈTRES ET LECTURE EXCEL
dossier_mere = "E:\Data\Aurelie\analysis\Final\assembly\nocues\";
fichier_excel = fullfile(dossier_mere, 'analysis_brainbow_best_final.xlsx');

if ~isfile(fichier_excel)
    error('Le fichier %s est introuvable.', fichier_excel);
end

opts = detectImportOptions(fichier_excel);
opts.VariableNamingRule = 'preserve';
Data_Excel = readtable(fichier_excel, opts);

if ~ismember('Chemin_Dossier', Data_Excel.Properties.VariableNames)
    error('La colonne "Chemin_Dossier" est absente du fichier Excel.');
end

% Vérification et création de la variable Animal_ID si manquante
if ~ismember('Animal_ID', Data_Excel.Properties.VariableNames)
    warning('Colonne Animal_ID absente de l''Excel. Création d''IDs fictifs pour tester le code.');
    Data_Excel.Animal_ID = categorical([1;1;2;2;3;4;4;5;5;6;6;7;7;7;7;8;8]);
else
    Data_Excel.Animal_ID = categorical(Data_Excel.Animal_ID); 
end
Data_Excel.Session_ID = categorical(1:height(Data_Excel))';

nb_sessions = height(Data_Excel);
fprintf('Lecture de %d sessions depuis l''Excel.\n', nb_sessions);
% %%%
% % Création d'une table simplifiée avec uniquement la métadonnée
% Table_Identifiers = table(Data_Excel.Session_ID, Data_Excel.Animal_ID, Data_Excel.Chemin_Dossier, ...
%     'VariableNames', {'Session_ID', 'Animal_ID', 'Chemin_Dossier'});
% 
% % Option 1 : Sauvegarde en .mat (Idéal pour recharger dans un autre script MATLAB)
% fichier_mat_ids = fullfile(dossier_mere, 'Session_Identifiers.mat');
% save(fichier_mat_ids, 'Table_Identifiers');
% fprintf('Identifiants sauvegardés en format MATLAB dans : %s\n', fichier_mat_ids);
% 
% %%%
%% 2. INITIALISATION DES VARIABLES
nb_assemblies_per_session = zeros(nb_sessions, 1);
pourcent_biaisees_per_session = zeros(nb_sessions, 1);

max_assemblies = nb_sessions * 20; 
cells_per_assembly = NaN(max_assemblies, 1);
spatial_dispersion = NaN(max_assemblies, 1);
correlation_assemblies = NaN(max_assemblies, 1);
sce_freq_assembly = NaN(max_assemblies, 1);
is_biased = NaN(max_assemblies, 1); 

% Variables de traçabilité pour LMM et SuperPlots
session_id_ass = NaN(max_assemblies, 1);
animal_id_ass = NaN(max_assemblies, 1);

index_ass = 1; 

%% 3. EXTRACTION DEPUIS LES .MAT
h_wait = waitbar(0, 'Chargement des données biologiques...');

for s = 1:nb_sessions
    chemin_dossier = Data_Excel.Chemin_Dossier{s};
    matFilePath = fullfile(chemin_dossier, 'results.mat');
    brainbowmatFilePath = fullfile(chemin_dossier, 'brainbow.mat');
    
    if ~isfile(matFilePath) || ~isfile(brainbowmatFilePath)
        continue;
    end
    
    data_res = load(matFilePath, 'assemblyortho');
    data_bb = load(brainbowmatFilePath, 'successrate', 'color_assembly', 'distance_assembly', 'freqsceok', 'correlation_assembly');
    
    if isempty(data_res.assemblyortho)
        nb_ass = 0;
    else
        nb_ass = length(data_res.assemblyortho);
    end
    nb_assemblies_per_session(s) = nb_ass;
    
    if isfield(data_bb, 'successrate') && ~isnan(data_bb.successrate)
        pourcent_biaisees_per_session(s) = data_bb.successrate;
    end
    
    curr_animal_id = double(Data_Excel.Animal_ID(s)); 
    
    for a = 1:nb_ass
        cells_per_assembly(index_ass) = length(data_res.assemblyortho{a});
        session_id_ass(index_ass) = double(Data_Excel.Session_ID(s));
        animal_id_ass(index_ass) = curr_animal_id; 
        
        if isfield(data_bb, 'color_assembly') && a <= length(data_bb.color_assembly)
            if isempty(data_bb.color_assembly{a}) || strcmp(data_bb.color_assembly{a}, 'Not significant')
                is_biased(index_ass) = 0;
            else
                is_biased(index_ass) = 1;
            end
        end
        
        if isfield(data_bb, 'distance_assembly') && a <= length(data_bb.distance_assembly)
            spatial_dispersion(index_ass) = data_bb.distance_assembly(a);
        end
        
        if isfield(data_bb, 'correlation_assembly') && a <= length(data_bb.correlation_assembly)
            correlation_assemblies(index_ass) = data_bb.correlation_assembly(a);
        end
        
        if isfield(data_bb, 'freqsceok') && a <= length(data_bb.freqsceok)
            sce_freq_assembly(index_ass) = data_bb.freqsceok(a);
        end
        
        index_ass = index_ass + 1;
    end
    waitbar(s / nb_sessions, h_wait);
end
close(h_wait);

% Nettoyage des NaN
valid_idx = ~isnan(cells_per_assembly);
cells_per_assembly = cells_per_assembly(valid_idx);
is_biased = is_biased(valid_idx);
spatial_dispersion = spatial_dispersion(valid_idx);
correlation_assemblies = correlation_assemblies(valid_idx);
sce_freq_assembly = sce_freq_assembly(valid_idx);
session_id_ass = categorical(session_id_ass(valid_idx));
animal_id_ass = categorical(animal_id_ass(valid_idx));
is_biased_cat = categorical(is_biased, [0 1], {'Mixed', 'Biased'});

%% 4. PREPARATION STATISTIQUES LMM
T_stats = table(cells_per_assembly, spatial_dispersion, correlation_assemblies, sce_freq_assembly, ...
    is_biased_cat, session_id_ass, animal_id_ass, ...
    'VariableNames', {'Cells', 'Dispersion', 'Correlation', 'SCE', 'BiasType', 'Session', 'Animal'});

%% 5. PARAMETRES GRAPHIQUES (Mise en page Article)
dossier_export = fullfile(dossier_mere, 'Figures_Article');
if ~exist(dossier_export, 'dir'), mkdir(dossier_export); end

format_panel = @(ax) set(ax, ...
    'FontName', 'Arial', 'FontSize', 11, 'LineWidth', 1.2, ...
    'Box', 'off', 'TickDir', 'out', 'XColor', 'k', 'YColor', 'k');

col_unbiased = [0.6 0.6 0.6]; % Gris
col_biased = [0.85 0.32 0.1]; % Rouge Brique

n_animals = length(categories(T_stats.Animal));
cmap_animals = lines(n_animals); % Palette de couleurs pour les animaux dans les SuperPlots

%% ========================================================================
% FIG A : Distribution générale (Nb assemblées et Nb Cellules)
% ========================================================================
figA = figure('Name', 'Stats Generales', 'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 14 10]);

ax1 = subplot(1,2,1);
boxchart(ones(size(nb_assemblies_per_session)), nb_assemblies_per_session, 'BoxFaceColor', 'k', 'MarkerStyle', 'none');
hold on;
swarmchart(ones(size(nb_assemblies_per_session)), nb_assemblies_per_session, 25, [0.4 0.4 0.4], 'filled', 'MarkerEdgeAlpha', 0.5);
ylabel('Assemblies per session'); xticks(1); xticklabels({''});
format_panel(ax1);

ax2 = subplot(1,2,2);
boxchart(ones(size(cells_per_assembly)), cells_per_assembly, 'BoxFaceColor', 'k', 'MarkerStyle', 'none');
hold on;
swarmchart(ones(size(cells_per_assembly)), cells_per_assembly, 25, [0.4 0.4 0.4], 'filled', 'MarkerEdgeAlpha', 0.5);
ylabel('Cells per assembly'); xticks(1); xticklabels({''});
format_panel(ax2);

exportgraphics(figA, fullfile(dossier_export, 'FigA_GeneralStats.svg'), 'ContentType', 'vector');
exportgraphics(figA, fullfile(dossier_export, 'FigA_GeneralStats.png'), 'Resolution', 600);

%% ========================================================================
% FIG B : Pourcentage d'assemblées avec un biais de couleur
% ========================================================================
figB = figure('Name', 'Proportion Biaisees', 'Color', 'w', 'Units', 'centimeters', 'Position', [17 5 7 10]);
ax_bias = axes(figB);

boxchart(ones(size(pourcent_biaisees_per_session)), pourcent_biaisees_per_session, ...
    'BoxFaceColor', col_biased, 'BoxFaceAlpha', 0.4, 'MarkerStyle', 'none', 'LineWidth', 1.5);
hold on;
swarmchart(ones(size(pourcent_biaisees_per_session)), pourcent_biaisees_per_session, ...
    30, col_biased, 'filled', 'MarkerEdgeAlpha', 0.6, 'XJitterWidth', 0.4);

ylabel('% assemblies with color bias', 'FontWeight', 'bold');
ylim([0 100]); xticks(1); xticklabels({'All sessions'});
format_panel(ax_bias);

exportgraphics(figB, fullfile(dossier_export, 'FigB_ProportionBiased.svg'), 'ContentType', 'vector');
exportgraphics(figB, fullfile(dossier_export, 'FigB_ProportionBiased.png'), 'Resolution', 600);

%% ========================================================================
% FIG C : SUPERPLOTS Biaisées vs Non-Biaisées (4 métriques avec LMM)
% ========================================================================
figC = figure('Name', 'SuperPlots Comparaison', 'Color', 'w', 'Units', 'centimeters', 'Position', [2 5 28 10]);

metrics = {'Cells', 'Dispersion', 'Correlation', 'SCE'};
ylabels = {'Cells per assembly', 'Spatial dispersion (\mum)', 'Mean pairwise correlation', 'SCE frequency (Hz)'};

for m = 1:4
    metric_name = metrics{m};
    ax = subplot(1, 4, m);
    hold on;
    
    % 1. STATISTIQUES LMM
    formule_lme = sprintf('%s ~ BiasType + (1|Animal) + (1|Animal:Session)', metric_name);
    lme_mod = fitlme(T_stats, formule_lme, 'FitMethod', 'reml');
    anova_res = anova(lme_mod, 'dfmethod', 'satterthwaite');
    p_val = anova_res.pValue(2); % P-value associée à 'BiasType'
    
    % 2. TRACER L'ARRIERE-PLAN (Données individuelles)
    for a = 1:n_animals
        curr_anim = categories(T_stats.Animal); curr_anim = curr_anim{a};
        
        idx_mixed = (T_stats.Animal == curr_anim) & (T_stats.BiasType == 'Mixed');
        if any(idx_mixed)
            swarmchart(ones(sum(idx_mixed),1), T_stats.(metric_name)(idx_mixed), ...
                10, cmap_animals(a,:), 'filled', 'XJitterWidth', 0.3, 'MarkerFaceAlpha', 0.2);
        end
        
        idx_biased = (T_stats.Animal == curr_anim) & (T_stats.BiasType == 'Biased');
        if any(idx_biased)
            swarmchart(2*ones(sum(idx_biased),1), T_stats.(metric_name)(idx_biased), ...
                10, cmap_animals(a,:), 'filled', 'XJitterWidth', 0.3, 'MarkerFaceAlpha', 0.2);
        end
    end
    
    % 3. TRACER LE PREMIER PLAN (Moyennes par animal et lignes connectrices)
    for a = 1:n_animals
        curr_anim = categories(T_stats.Animal); curr_anim = curr_anim{a};
        
        % omitnan évite un plantage si un animal n'a pas du tout d'assemblées d'une des catégories
        mean_mixed = mean(T_stats.(metric_name)((T_stats.Animal == curr_anim) & (T_stats.BiasType == 'Mixed')), 'omitnan');
        mean_biased = mean(T_stats.(metric_name)((T_stats.Animal == curr_anim) & (T_stats.BiasType == 'Biased')), 'omitnan');
        
        plot(1, mean_mixed, 'o', 'MarkerSize', 9, 'MarkerFaceColor', cmap_animals(a,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
        plot(2, mean_biased, 'o', 'MarkerSize', 9, 'MarkerFaceColor', cmap_animals(a,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
        
        if ~isnan(mean_mixed) && ~isnan(mean_biased)
            plot([1 2], [mean_mixed mean_biased], '-', 'Color', [cmap_animals(a,:) 0.6], 'LineWidth', 1.5);
        end
    end
    
    % 4. FORMATAGE DU GRAPHIQUE
    ylabel(ylabels{m});
    xticks([1 2]); xticklabels({'Mixed', 'Biased'});
    xlim([0.5 2.5]);
    
    % Format de l'affichage de la p-value
    if p_val < 0.001
        title('LMM p < 0.001', 'FontWeight', 'normal');
    else
        title(sprintf('LMM p = %.3f', p_val), 'FontWeight', 'normal');
    end
    format_panel(ax);
end

% Export
exportgraphics(figC, fullfile(dossier_export, 'FigC_ComparisonBiased.svg'), 'ContentType', 'vector');
exportgraphics(figC, fullfile(dossier_export, 'FigC_ComparisonBiased.png'), 'Resolution', 600);