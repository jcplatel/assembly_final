%check colorcell
clear
load("E:\Data\Aurelie\data\nocues\nocues.mat");%PC
color1 = load("E:\Data\Aurelie\data\chroms\119\220919\registration\colorcell_biapy_2026_05_27.mat");
color2 = load( "E:\Data\Aurelie\data\chroms\119\220923\registration\color_matched_biapyJC_26_06_10.mat");
color3 = load("E:\Data\Aurelie\data\chroms\582\230331plane0\colorcell_exvivo_JC_20260622.mat");


%% 

% files_to_process = [1:6,8:76] ;%1:numel (matfile) ;%22;%  [1:3,23:30];%1:numel (matfile) 
files_to_process =22% 1:numel (matfile);% 7    38    39    40    41    42    43    44    45

for file_num = files_to_process
    try 
    file_num
    filename=string(matfile{file_num});
    [path,name,ext] = fileparts(filename);
    
    if strlength(name)> 23; identifier = extractBetween(name,5,24);else identifier='suite2p';end
    id_propre = strrep(identifier, '_', '-'); 
    path_colorcell=fullfile(path,"colorcell_Maxwell.mat");
    load (path_colorcell)
    fig = proportion_couleur(colorcell,'test',1)
    % p_value_chi2 = proportion(colorcell);
    % p_value_chi2_all (:,file_num) = p_value_chi2;
    nam_all (:,file_num) = identifier;
    alertes = verification_proportions_colorcell(colorcell);
    alertes_all(file_num) = alertes;
    catch 
        error(file_num) = 1; 
        disp (strcat('error',num2str(file_num)))
    end

end
Table_Nouvelle_Exp = table(nam_all',alertes_all');

function p_value = proportion(colorcell)

    % --- 1. DONNÉES OBSERVÉES ---
    % On s'intéresse aux 6 couleurs principales
    counts_obs = zeros(1, 6); 
    
    counts_obs(1) = sum(colorcell == 1); % Red
    counts_obs(2) = sum(colorcell == 2); % Green
    counts_obs(3) = sum(colorcell == 3); % Blue
    counts_obs(4) = sum(colorcell == 4); % Yellow
    counts_obs(5) = sum(colorcell == 5); % Magenta
    counts_obs(6) = sum(colorcell == 6); % Cyan
    
    % Nombre total de cellules colorées (on ignore le Blanc=7 et le Noir=8)
    N_total = sum(counts_obs);
    percent_obs = (counts_obs / N_total) * 100;
    
    % --- 2. DONNÉES THÉORIQUES (CONTRÔLE) ---
    props_theo_raw = [28.4, 7.7, 3.8, 29.4, 20.8, 10.8];
    %[34, 5.6, 4.3, 27.3, 19.3, 9.4];
    %[34+8, 0-17, 0-11, 15-43, 10-27, 2-20];
    props_theo_norm = props_theo_raw / sum(props_theo_raw);
    percent_theo = props_theo_norm * 100;
    counts_expected = N_total * props_theo_norm;
    
    % --- 3. AFFICHAGE DES POURCENTAGES ---
    noms_couleurs = {'Red', 'Green', 'Blue', 'Yellow', 'Magenta', 'Cyan'};
    
    fprintf('====== COMPARAISON DES PROPORTIONS ======\n');
    fprintf('%-10s | %-12s | %-12s\n', 'Couleur', 'Observé (%)', 'Contrôle (%)');
    fprintf('-----------------------------------------\n');
    for i = 1:6
        fprintf('%-10s | %5.1f %% | %5.1f %%\n', ...
            noms_couleurs{i}, percent_obs(i), percent_theo(i));
    end
    fprintf('-----------------------------------------\n');
    
    % --- 4. TEST STATISTIQUE (CHI-SQUARE) ---
    chi2_stat = sum( ((counts_obs - counts_expected).^2) ./ counts_expected );
    df = length(counts_obs) - 1;
    p_value = 1 - chi2cdf(chi2_stat, df);
    
    fprintf('\n====== RÉSULTAT DU TEST STATISTIQUE ======\n');
    fprintf('Chi-2 Statistique = %.2f\n', chi2_stat);
    fprintf('Degrés de liberté = %d\n', df);
    fprintf('p-value = %.2e\n', p_value);
    
    if p_value < 0.05
        fprintf('\nCONCLUSION : Il y a une différence STATISTIQUEMENT SIGNIFICATIVE\n');
        fprintf('entre vos proportions et les proportions de contrôle.\n');
    else
        fprintf('\nCONCLUSION : Aucune différence significative (p >= 0.05).\n');
        fprintf('Vos proportions sont conformes au contrôle.\n');
    end
    
    % --- 5. GRAPHIQUE BAR PLOT (AVEC LES VRAIES COULEURS) ---
    
    % Définition des couleurs RGB pour [Red, Green, Blue, Yellow, Magenta, Cyan]
    true_colors = [
        1.0, 0.0, 0.0; % Red
        0.0, 0.8, 0.0; % Green (légèrement assombri pour être visible sur fond blanc)
        0.0, 0.0, 1.0; % Blue
        0.9, 0.9, 0.0; % Yellow (légèrement assombri pour être lisible)
        1.0, 0.0, 1.0; % Magenta
        0.0, 1.0, 1.0  % Cyan
    ];
return
    figure('Name', 'Comparaison des Proportions', 'Color', 'w');
    hold on;
    
    % Largeur des barres et espacement
    bar_width = 0.35;
    x = 1:6; % Positions sur l'axe X
    
    % Création des barres "Observé" (avec les vraies couleurs pleines)
    for i = 1:6
        b1 = bar(x(i) - bar_width/2, percent_obs(i), bar_width);
        b1.FaceColor = true_colors(i, :);
        b1.EdgeColor = 'k'; % Bord noir pour délimiter
    end
    
    % Création des barres "Contrôle" (avec les vraies couleurs, mais transparentes/hachurées)
    % Pour MATLAB, on utilise l'attribut FaceAlpha pour rendre le contrôle semi-transparent
    for i = 1:6
        b2 = bar(x(i) + bar_width/2, percent_theo(i), bar_width);
        b2.FaceColor = true_colors(i, :);
        b2.EdgeColor = 'k';
        b2.FaceAlpha = 0.3; % Rend la couleur pâle (30% d'opacité) pour différencier du groupe "Observé"
    end
    
    % Paramétrage de l'axe X
    set(gca, 'XTick', x, 'XTickLabel', noms_couleurs);
    ylabel('Proportion (%)');
    title('Proportions Observées (Plein) vs Théoriques (Transparent)');
    
    % Création d'une légende "fantôme" pour expliquer la différence opacité/plein
    h_obs = bar(nan, nan, 'FaceColor', 'k', 'FaceAlpha', 1);
    h_theo = bar(nan, nan, 'FaceColor', 'k', 'FaceAlpha', 0.3);
    legend([h_obs, h_theo], {'Observé', 'Théorique'}, 'Location', 'northeast');
    
    grid on;
    hold off;

end

function alertes = verification_proportions_colorcell(colorcell)
    % Convertir colorcell en format numérique s'il est en categorical/string
    if iscategorical(colorcell) || isstring(colorcell) || iscell(colorcell)
        colorcell_num = str2double(string(colorcell));
    else
        colorcell_num = double(colorcell);
    end

    % 1. Calcul des proportions à partir de colorcell
    % Vos classes : 1=rouge, 2=vert, 3=bleu, 4=jaune, 5=magenta, 6=cyan, 7=blanc, 8=noir
    % On ne garde que les cellules des 6 couleurs principales (1 à 6)
    masque_valides = (colorcell_num >= 1 & colorcell_num <= 6);
    N_total_valides = sum(masque_valides);
    
    if N_total_valides == 0
        error('Aucune cellule valide (couleurs 1 à 6) trouvée dans colorcell.');
    end

    % Extraction des pourcentages (Attention à l'ordre pour matcher le graphe)
    prop_rouge   = sum(colorcell_num == 1) / N_total_valides * 100;
    prop_vert    = sum(colorcell_num == 2) / N_total_valides * 100;
    prop_bleu    = sum(colorcell_num == 3) / N_total_valides * 100;
    prop_magenta = sum(colorcell_num == 5) / N_total_valides * 100; % 5 = magenta
    prop_jaune   = sum(colorcell_num == 4) / N_total_valides * 100; % 4 = jaune
    prop_cyan    = sum(colorcell_num == 6) / N_total_valides * 100; % 6 = cyan
    
    proportions_champ_test = [prop_rouge, prop_vert, prop_bleu, prop_magenta, prop_jaune, prop_cyan];

    % 2. Définition des distributions de référence (issues du graphique)
    couleurs = {'Rouge', 'Vert', 'Bleu', 'Magenta', 'Jaune', 'Cyan'};
    moyennes_ref = [34.0, 7.0, 5.0, 28.0, 18.0, 11.0];
    sigma_ref = [10.0, 7.0, 5.0, 13.0, 8.0, 9]; % Écart-types estimés
    %[34+8, 0-17, 0-11, 15-43, 10-27, 2-20];

    % 3. Analyse et Affichage
    fprintf('\n=== ANALYSE DES PROPORTIONS DE COLORCELL ===\n');
    fprintf('Nombre de cellules analysées : %d (hors blanc/noir)\n\n', N_total_valides);
    fprintf('%-10s | %-12s | %-15s | %-10s | %-15s\n', ...
            'Couleur', 'Votre Champ', 'Réf (Moy ± \sigma)', 'Z-score', 'Statut');
    fprintf(repmat('-', 1, 72));
    fprintf('\n');

    alertes = 0; % Compteur d'anomalies

    for i = 1:6
        val_test = proportions_champ_test(i);
        val_ref = moyennes_ref(i);
        std_ref = sigma_ref(i);
        
        % Calcul du Z-score et p-value
        z_score = (val_test - val_ref) / std_ref;
        p_value = 2 * (1 - normcdf(abs(z_score)));
        
        % Statut (seuil p < 0.05, soit Z > 1.96)
        if abs(z_score) <= 1.96
            statut = 'Dans la norme';
        else
            statut = '!! ANOMALIE !!';
            alertes = alertes + 1;
        end
        
        % Affichage de la ligne
        fprintf('%-10s | %6.1f %%    | %5.1f ± %-4.1f %%  | %6.2f     | %-15s (p=%.3f)\n', ...
                couleurs{i}, val_test, val_ref, std_ref, z_score, statut, p_value);
    end
    
    fprintf(repmat('-', 1, 72));
    fprintf('\n');

    % 4. Conclusion
    if alertes == 0
        fprintf('Conclusion : Les proportions extraites de colorcell sont plausibles.\n\n');
    else
        fprintf('Conclusion : %d proportion(s) s''éloigne(nt) de la norme de référence.\n\n', alertes);
    end
end