% =========================================================================
% SCRIPT : FIGURES ARTICLE - CARACTÉRISATION DES ASSEMBLÉES
% Lecture depuis l'Excel des meilleurs K et extraction des .mat
% Exportation standardisée : 600 DPI, SVG, Police Arial
% =========================================================================

clear; clc; close all;

%% 1. PARAMÈTRES ET LECTURE EXCEL
dossier_mere = "E:\Data\Aurelie\analysis\Final\assembly\nocues\";
fichier_excel = fullfile(dossier_mere, 'analysis_brainbow_best_final.xlsx');

if ~isfile(fichier_excel)
    error('Le fichier %s est introuvable.', fichier_excel);
end

% Lecture de la table (conserve les noms de colonnes)
opts = detectImportOptions(fichier_excel);
opts.VariableNamingRule = 'preserve';
Data_Excel = readtable(fichier_excel, opts);

% On s'assure que la colonne Chemin_Dossier existe
if ~ismember('Chemin_Dossier', Data_Excel.Properties.VariableNames)
    error('La colonne "Chemin_Dossier" est absente du fichier Excel.');
end

% Filtrer pour garder uniquement les lignes que vous avez validées/choisies 
% (Par exemple, si vous avez ajouté une colonne "Garder" avec des 1 et 0, 
% décommentez la ligne ci-dessous)
% Data_Excel = Data_Excel(Data_Excel.Garder == 1, :);

nb_sessions = height(Data_Excel);
fprintf('Lecture de %d sessions depuis l''Excel.\n', nb_sessions);

%% 2. INITIALISATION DES VARIABLES POUR LES FIGURES
% Stats par session
nb_assemblies_per_session = zeros(nb_sessions, 1);
pourcent_biaisees_per_session = zeros(nb_sessions, 1);

% Stats par assemblée (on pré-alloue grand, on coupera à la fin)
max_assemblies = nb_sessions * 20; 
cells_per_assembly = NaN(max_assemblies, 1);
spatial_dispersion = NaN(max_assemblies, 1);
correlation_assemblies = NaN(max_assemblies, 1);
is_biased = NaN(max_assemblies, 1); % 1 = Biaisée (colorée), 0 = Mixte

index_ass = 1; % Compteur global d'assemblées

%% 3. BOUCLE D'EXTRACTION DES DONNÉES DEPUIS LES .MAT
h_wait = waitbar(0, 'Chargement des données biologiques...');

for s = 1:nb_sessions
    chemin_dossier = Data_Excel.Chemin_Dossier{s};
    
    matFilePath = fullfile(chemin_dossier, 'results.mat');
    brainbowmatFilePath = fullfile(chemin_dossier, 'brainbow.mat');
    
    if ~isfile(matFilePath) || ~isfile(brainbowmatFilePath)
        fprintf('⚠️ Fichiers manquants pour la session %d : %s\n', s, chemin_dossier);
        continue;
    end
    
    % Chargement des variables dont on a besoin
    % - Depuis results.mat : les cellules par assemblée, etc.
    % - Depuis brainbow.mat : distance, couleurs, succès...
    data_res = load(matFilePath, 'assemblyortho');
    data_bb = load(brainbowmatFilePath, 'successrate', 'color_assembly', 'distance_assembly');
    
    % 1. Nombre d'assemblées dans cette session
    if isempty(data_res.assemblyortho)
        nb_ass = 0;
    else
        nb_ass = length(data_res.assemblyortho);
    end
    nb_assemblies_per_session(s) = nb_ass;
    
    % 2. Pourcentage de colorées (biaisées) pour la session entière
    if isfield(data_bb, 'successrate') && ~isnan(data_bb.successrate)
        pourcent_biaisees_per_session(s) = data_bb.successrate;
    end
    
    % 3. Extraction des données par assemblée (Boucle sur les assemblées du K)
    for a = 1:nb_ass
        % Taille (Nb cellules)
        cells_per_assembly(index_ass) = length(data_res.assemblyortho{a});
        
        % Est-elle biaisée ? (on regarde si elle a une couleur assignée)
        if isfield(data_bb, 'color_assembly') && a <= length(data_bb.color_assembly)
            if isempty(data_bb.color_assembly{a}) || strcmp(data_bb.color_assembly{a}, 'Not significant')
                is_biased(index_ass) = 0;
            else
                is_biased(index_ass) = 1;
            end
        end
        
        % Dispersion spatiale de l'assemblée
        if isfield(data_bb, 'distance_assembly') && a <= length(data_bb.distance_assembly)
            spatial_dispersion(index_ass) = data_bb.distance_assembly(a);
        end
        
        % Corrélation (À AJOUTER SI VOUS L'AVEZ CALCULÉE DANS BRAINBOW.MAT)
        % Si vous avez une variable 'correlation_assembly' par exemple :
        % if isfield(data_bb, 'correlation_assembly')
        %    correlation_assemblies(index_ass) = data_bb.correlation_assembly(a);
        % end
        
        index_ass = index_ass + 1;
    end
    
    waitbar(s / nb_sessions, h_wait);
end
close(h_wait);

% Nettoyage des NaN (on coupe les tableaux pré-alloués)
valid_idx = ~isnan(cells_per_assembly);
cells_per_assembly = cells_per_assembly(valid_idx);
is_biased = is_biased(valid_idx);
spatial_dispersion = spatial_dispersion(valid_idx);
correlation_assemblies = correlation_assemblies(valid_idx);

fprintf('Extraction terminée : %d assemblées traitées.\n', sum(valid_idx));

%% 4. FONCTION DE STANDARDISATION (Mise en page Article)
dossier_export = fullfile(dossier_mere, 'Figures_Article');
if ~exist(dossier_export, 'dir'), mkdir(dossier_export); end

format_panel = @(ax) set(ax, ...
    'FontName', 'Arial', 'FontSize', 11, 'LineWidth', 1.2, ...
    'Box', 'off', 'TickDir', 'out', 'XColor', 'k', 'YColor', 'k');

col_unbiased = [0.6 0.6 0.6]; % Gris
col_biased = [0.85 0.32 0.1]; % Rouge Brique

%% ========================================================================
% FIG A : Distribution générale (Nb assemblées et Nb Cellules)
% ========================================================================
figA = figure('Name', 'Stats Generales', 'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 14 10]);

ax1 = subplot(1,2,1);
boxchart(ones(size(nb_assemblies_per_session)), nb_assemblies_per_session, 'BoxFaceColor', 'k', 'MarkerStyle', 'none');
hold on;
swarmchart(ones(size(nb_assemblies_per_session)), nb_assemblies_per_session, 25, [0.4 0.4 0.4], 'filled', 'MarkerEdgeAlpha', 0.5);
ylabel('Assemblies per session');
xticks(1); xticklabels({''});
format_panel(ax1);

ax2 = subplot(1,2,2);
boxchart(ones(size(cells_per_assembly)), cells_per_assembly, 'BoxFaceColor', 'k', 'MarkerStyle', 'none');
hold on;
swarmchart(ones(size(cells_per_assembly)), cells_per_assembly, 25, [0.4 0.4 0.4], 'filled', 'MarkerEdgeAlpha', 0.5);
ylabel('Cells per assembly');
xticks(1); xticklabels({''});
format_panel(ax2);

exportgraphics(figA, fullfile(dossier_export, 'FigA_GeneralStats.svg'), 'ContentType', 'vector');
exportgraphics(figA, fullfile(dossier_export, 'FigA_GeneralStats.png'), 'Resolution', 600);

%% ========================================================================
% FIG B : Proportion d'assemblées biaisées
% ========================================================================
figB = figure('Name', 'Proportion Biaisees', 'Color', 'w', 'Units', 'centimeters', 'Position', [17 2 8 10]);

ax3 = axes(figB);
boxchart(ones(size(pourcent_biaisees_per_session)), pourcent_biaisees_per_session, 'BoxFaceColor', col_biased, 'MarkerStyle', 'none');
hold on;
swarmchart(ones(size(pourcent_biaisees_per_session)), pourcent_biaisees_per_session, 25, col_biased*0.8, 'filled', 'MarkerEdgeAlpha', 0.5);
ylabel('% of color-biased assemblies');
ylim([0 max(100, max(pourcent_biaisees_per_session)+10)]);
xticks(1); xticklabels({''});
format_panel(ax3);

exportgraphics(figB, fullfile(dossier_export, 'FigB_ProportionBiased.svg'), 'ContentType', 'vector');
exportgraphics(figB, fullfile(dossier_export, 'FigB_ProportionBiased.png'), 'Resolution', 600);

%% ========================================================================
% FIG C : Biaisées vs Non-Biaisées (Taille et Dispersion)
% ========================================================================
% (Ici on fait 2 sous-graphes car la corrélation n'est pas encore chargée)
figC = figure('Name', 'Comparaison Biais vs Non Biais', 'Color', 'w', 'Units', 'centimeters', 'Position', [2 14 16 10]);

idx_bias = (is_biased == 1);
idx_unbias = (is_biased == 0);

groupe_tailles = [ones(sum(idx_unbias),1) ; 2*ones(sum(idx_bias),1)];
donnees_tailles = [cells_per_assembly(idx_unbias) ; cells_per_assembly(idx_bias)];

groupe_disp = [ones(sum(idx_unbias),1) ; 2*ones(sum(idx_bias),1)];
donnees_disp = [spatial_dispersion(idx_unbias) ; spatial_dispersion(idx_bias)];

% C1. Taille (Nb Cellules)
ax4 = subplot(1,2,1);
b1 = boxchart(groupe_tailles, donnees_tailles, 'GroupByColor', groupe_tailles, 'MarkerStyle', 'none');
b1(1).BoxFaceColor = col_unbiased; 
b1(2).BoxFaceColor = col_biased;   
hold on;
swarmchart(groupe_tailles(groupe_tailles==1), donnees_tailles(groupe_tailles==1), 15, col_unbiased*0.8, 'filled', 'XJitterWidth', 0.4);
swarmchart(groupe_tailles(groupe_tailles==2), donnees_tailles(groupe_tailles==2), 15, col_biased*0.8, 'filled', 'XJitterWidth', 0.4);
ylabel('Cells per assembly');
xticks([1 2]); xticklabels({'Mixed', 'Biased'});
format_panel(ax4);

% C2. Dispersion Spatiale
ax5 = subplot(1,2,2);
b2 = boxchart(groupe_disp, donnees_disp, 'GroupByColor', groupe_disp, 'MarkerStyle', 'none');
b2(1).BoxFaceColor = col_unbiased; 
b2(2).BoxFaceColor = col_biased;   
hold on;
swarmchart(groupe_disp(groupe_disp==1), donnees_disp(groupe_disp==1), 15, col_unbiased*0.8, 'filled', 'XJitterWidth', 0.4);
swarmchart(groupe_disp(groupe_disp==2), donnees_disp(groupe_disp==2), 15, col_biased*0.8, 'filled', 'XJitterWidth', 0.4);
ylabel('Spatial dispersion (\mum)');
xticks([1 2]); xticklabels({'Mixed', 'Biased'});
format_panel(ax5);

exportgraphics(figC, fullfile(dossier_export, 'FigC_ComparisonBiased.svg'), 'ContentType', 'vector');
exportgraphics(figC, fullfile(dossier_export, 'FigC_ComparisonBiased.png'), 'Resolution', 600);

fprintf('✅ Les 3 figures de l''article sont exportées dans : %s\n', dossier_export);