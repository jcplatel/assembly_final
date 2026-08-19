% =========================================================================
% SCRIPT : FIGURES ARTICLE - CORRÉLATIONS (SAME vs DIFF)
% Lecture du résumé des corrélations (summary_animals.csv)
% Filtrage selon les sessions validées (Excel)
% Modèle LMM et SuperPlot (Same Color vs Different Color)
% =========================================================================

clear; clc; close all;

%% 1. PARAMÈTRES ET CHEMINS
dossier_corr = "E:\Data\Aurelie\analysis\Final\correlation\";
fichier_summary = fullfile(dossier_corr, 'summary_animals.csv');
fichier_excel = "E:\Data\Aurelie\analysis\Final\assembly\nocues\analysis_brainbow_best_final.xlsx";

if ~isfile(fichier_summary)
    error('Le fichier de corrélation %s est introuvable.', fichier_summary);
end
if ~isfile(fichier_excel)
    error('Le fichier Excel des sessions validées %s est introuvable.', fichier_excel);
end

% Charger les résumés de corrélation
Opts_csv = detectImportOptions(fichier_summary);
T_Corr = readtable(fichier_summary, Opts_csv);

% Charger l'Excel pour filtrer et récupérer Animal_ID
Opts_xls = detectImportOptions(fichier_excel);
Opts_xls.VariableNamingRule = 'preserve';
Data_Excel = readtable(fichier_excel, Opts_xls);

if ~ismember('Animal_ID', Data_Excel.Properties.VariableNames)
    warning('Animal_ID manquant dans l''Excel. Création d''IDs fictifs.');
    Data_Excel.Animal_ID = categorical(randi([1 5], height(Data_Excel), 1));
else
    Data_Excel.Animal_ID = categorical(Data_Excel.Animal_ID); 
end
Data_Excel.Session_ID = categorical(1:height(Data_Excel))';

%% 2. FILTRAGE ET FUSION DES DONNÉES
% Votre fichier correlation contient "identifier".
% Il faut s'assurer que ce "identifier" matche ou peut être relié aux chemins de l'Excel.
% Ici, on va faire une jointure simple si possible.
% (NB: adaptez la méthode de jointure selon le nom exact de vos identifiants)

% Supposons que l'on extrait un "identifier" court depuis le Chemin_Dossier de l'Excel
Data_Excel.identifier_court = string(cellfun(@(x) extractBetween_custom(x), Data_Excel.Chemin_Dossier, 'UniformOutput', false));

% Jointure
T_Merged = innerjoin(T_Corr, Data_Excel, 'LeftKeys', 'identifier', 'RightKeys', 'identifier_court');
nb_sessions = height(T_Merged);
fprintf('Sessions conservées après filtrage : %d\n', nb_sessions);

%% 3. REFORMATAGE LONG POUR LMM (Linear Mixed Model)
% Le LMM a besoin des données au format "Long" (une ligne = une observation)
% On va empiler Same_Pearson et Diff_Pearson

SessionArray = repmat(T_Merged.Session_ID, 2, 1);
AnimalArray  = repmat(T_Merged.Animal_ID, 2, 1);

% Vecteur des valeurs (Pearson)
CorrValues = [T_Merged.mean_same_pearson ; T_Merged.mean_diff_pearson];

% Vecteur des catégories
ConditionArray = categorical([repmat({'Same Color'}, nb_sessions, 1) ; repmat({'Different Color'}, nb_sessions, 1)]);

T_LMM = table(CorrValues, ConditionArray, SessionArray, AnimalArray, ...
    'VariableNames', {'Correlation', 'PairType', 'Session', 'Animal'});

% Nettoyage NaN potentiels
T_LMM = rmmissing(T_LMM);

%% 4. STATISTIQUES LMM
fprintf('\n--- Modèle Linéaire Mixte (LMM) ---\n');
% Effet fixe: PairType. Effets aléatoires: Intercept par Animal, et Session dans Animal
lme = fitlme(T_LMM, 'Correlation ~ PairType + (1|Animal) + (1|Animal:Session)', 'FitMethod', 'reml');
anova_res = anova(lme, 'dfmethod', 'satterthwaite');
p_val = anova_res.pValue(2);
disp(anova_res);

%% 5. SUPERPLOT (Same vs Diff)
format_panel = @(ax) set(ax, 'FontName', 'Arial', 'FontSize', 11, 'LineWidth', 1.2, ...
    'Box', 'off', 'TickDir', 'out', 'XColor', 'k', 'YColor', 'k');

n_animals = length(categories(T_LMM.Animal));
cmap_animals = lines(n_animals);

figure('Name', 'SuperPlot Corrélation', 'Color', 'w', 'Units', 'centimeters', 'Position', [5 5 12 10]);
hold on;

x_cats = {'Same Color', 'Different Color'};

for a = 1:n_animals
    curr_anim = categories(T_LMM.Animal); curr_anim = curr_anim{a};
    
    % Extraction des données pour cet animal
    idx_same = (T_LMM.Animal == curr_anim) & (T_LMM.PairType == 'Same Color');
    idx_diff = (T_LMM.Animal == curr_anim) & (T_LMM.PairType == 'Different Color');
    
    val_same = T_LMM.Correlation(idx_same);
    val_diff = T_LMM.Correlation(idx_diff);
    
    % Swarmcharts (Arrière-plan)
    if any(idx_same)
        swarmchart(ones(sum(idx_same),1), val_same, 15, cmap_animals(a,:), 'filled', 'XJitterWidth', 0.2, 'MarkerFaceAlpha', 0.4);
    end
    if any(idx_diff)
        swarmchart(2*ones(sum(idx_diff),1), val_diff, 15, cmap_animals(a,:), 'filled', 'XJitterWidth', 0.2, 'MarkerFaceAlpha', 0.4);
    end
    
    % Moyennes (Premier plan)
    mean_s = mean(val_same);
    mean_d = mean(val_diff);
    
    plot(1, mean_s, 'o', 'MarkerSize', 8, 'MarkerFaceColor', cmap_animals(a,:), 'MarkerEdgeColor', 'k');
    plot(2, mean_d, 'o', 'MarkerSize', 8, 'MarkerFaceColor', cmap_animals(a,:), 'MarkerEdgeColor', 'k');
    
    % Lignes de connexion
    if ~isnan(mean_s) && ~isnan(mean_d)
        plot([1 2], [mean_s mean_d], '-', 'Color', [cmap_animals(a,:) 0.6], 'LineWidth', 1.5);
    end
end

xticks([1 2]);
xticklabels(x_cats);
xlim([0.5 2.5]);
ylabel('Pearson Correlation');
title(sprintf('LMM p-value = %.4f', p_val), 'FontWeight', 'normal');
format_panel(gca);

dossier_export = fullfile(dossier_corr, 'Figures_Article');
if ~exist(dossier_export, 'dir'), mkdir(dossier_export); end
exportgraphics(gcf, fullfile(dossier_export, 'SuperPlot_Correlation_Pearson.svg'), 'ContentType', 'vector');
exportgraphics(gcf, fullfile(dossier_export, 'SuperPlot_Correlation_Pearson.png'), 'Resolution', 600);

% --- Fonction utilitaire pour extraire le nom du dossier ---
function id = extractBetween_custom(path_str)
    % A adapter selon le format exact de vos dossiers
    % Ex: si le dossier est '.../P2M_444118_220914/'
    parts = strsplit(path_str, filesep);
    id = parts{end-1}; % ou parts{end} selon s'il y a un / à la fin
end