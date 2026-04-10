%% Consensus Clustering

assemblyraw = [];
NRace= size(Race,2);

[IDX2,sCl,M,S,NClini] = kmeansoptnew(Race, kmean_iter, 'var', NClini,[]);
NCl = NClini;

[~,x2] = sort(IDX2); % cluster de SCE
MSort = M(x2,x2);

%% Remove cluster non-statistically significant
% (Le reste de ton code reste identique car il travaille sur IDX2 et M)
% ...
sClrnd = zeros(1,kmeans_surrogate);
% cRace = parallel.pool.Constant(Race); % Toujours utiliser les vraies données pour le shuffle
cRace = parallel.pool.Constant(Race); % Toujours utiliser les vraies données pour le shuffle

All_Shuffle_Scores = [];

parfor i = 1:kmeans_surrogate  
    % Attention : ici on utilise toujours la fonction shuffle standard
    % car faire du consensus clustering DANS le shuffle serait trop long
    scores_du_run = kmeansoptrndnew(cRace.Value, kmeans_rnd_iter, NCl);
    All_Shuffle_Scores = [All_Shuffle_Scores, scores_du_run]; 
end


threshold = prctile(All_Shuffle_Scores, 95); 
NClOK = sum(sCl > threshold);
% NClOK =sum(sCl>prctile(sClrnd,95)); %use max, use 99%  ?

sClOK = sCl(1:NClOK)';

% NClOK = NCl;
% sClOK = sCl(1:NClOK)';

RaceOK = Race(:,IDX2<=NClOK);%ici on revient sur la matrice originale 
NRaceOK = size(RaceOK,2);
% disp(['nSCEOK: ' num2str(NRaceOK)])    


if NClOK>1
     RACE_Orthojcnew2
     NClOK=NCl;
     %%call rasterplot here
else
    NCl=NClOK;
    assemblyortho= cell(0);
    assemblystat= cell(0);
    assemblyraw = [];
end

%% recalcul silhouette cluster finaux
if NCl>1
% 1. Identifier les indices qui appartiennent aux clusters valides (de 1 à NCl)
    valid_indices = (IDX2 <= NCl) & (IDX2 > 0); 
    
    % 2. Filtrer IDX2 et récupérer les indices originaux
    IDX2_ok = IDX2(valid_indices); % Garde seulement les labels valides
    
    % 3. Récupérer les colonnes correspondantes dans M
    % Note: on n'a pas besoin de 'x2' (le tri) pour extraire la sous-matrice,
    % on utilise directement le masque booléen.
    MSort_ok = M(valid_indices, valid_indices);
    
    % Si vous avez besoin que ce soit trié par numéro de cluster pour silh() :
    [IDX2_sorted_ok, sort_idx] = sort(IDX2_ok);
    MSort_sorted_ok = MSort_ok(sort_idx, sort_idx);
    
    % Calcul de la silhouette
    s = silh(MSort_sorted_ok, IDX2_sorted_ok);
    
    sClOK = zeros(1, NCl);
    for i = 1:NCl
        % Vérifier si le cluster n'est pas vide avant de calculer la médiane
        if any(IDX2_sorted_ok == i)
            sClOK(i) = median(s(IDX2_sorted_ok == i)); 
        else
            sClOK(i) = NaN; % Ou 0, selon votre préférence
        end
    end
    mean_sClOK = mean(sClOK, 'omitnan'); % Moyenne en ignorant les NaN
else
    mean_sClOK = 'NA';
end
%

%% ========================================================================
%% NOUVEAU BLOC : Calcul de la Spécificité Spatio-Temporelle (Purity & Fidelity)
%% ========================================================================
if NClOK > 0 && exist('assemblyortho', 'var') && ~isempty(assemblyortho)
    Event_Purity      = NaN(1, NClOK); % Vertical
    Cellular_Fidelity = NaN(1, NClOK); % Horizontal
    
    % Masque global de tous les SCE appartenant à un cluster validé
    valid_sce_mask = (IDX2 <= NClOK) & (IDX2 > 0);
    
    for c = 1:NClOK
        cells_c = assemblyortho{c};
        idx_SCE_c = (IDX2 == c); % Masque des SCEs affectés à ce cluster 'c'
        
        if sum(idx_SCE_c) > 0 && ~isempty(cells_c)
            
            % --- CALCUL VERTICAL (Event Purity) ---
            % Total des tirs de TOUTES les cellules pendant les SCE de 'c'
            tirs_totaux_colonne = sum(Race(:, idx_SCE_c), 'all');
            
            % Tirs UNIQUEMENT des cellules de l'assemblée 'c' pendant ces SCE
            tirs_in_box = sum(Race(cells_c, idx_SCE_c), 'all');
            
            if tirs_totaux_colonne > 0
                Event_Purity(c) = tirs_in_box / tirs_totaux_colonne;
            end
            
            % --- CALCUL HORIZONTAL (Cellular Fidelity) ---
            % Tirs des cellules de l'assemblée 'c' dans TOUS les SCE valides (tous clusters)
            tirs_totaux_ligne = sum(Race(cells_c, valid_sce_mask), 'all');
            
            if tirs_totaux_ligne > 0
                Cellular_Fidelity(c) = tirs_in_box / tirs_totaux_ligne;
            end
            
        end
    end
    
    % On calcule la moyenne pour ce run
    mean_Event_Purity      = mean(Event_Purity, 'omitnan');
    mean_Cellular_Fidelity = mean(Cellular_Fidelity, 'omitnan');
else
    mean_Event_Purity      = NaN;
    mean_Cellular_Fidelity = NaN;
end
%% ========================================================================

if savefig==1 &&  NClOK>1
    figure('visible','off');
    % figure
    subplot(1,2,1)
    imagesc(MSort)
    colormap jet
    axis image
    xlabel('sorted SCE #')     %was RACE
    ylabel('sorted SCE #')     %was RACE

    subplot(1,2,2)
    % imagesc(Race(x1,x2),[-1 1.2])
    imagesc(Race(x1,RList),[-1 1.2])
    axis image
    xlabel('sorted SCE #')     %was RACE
    ylabel('sorted Cell #')
    exportgraphics(gcf,strcat(namefull ,'clusters.png'),'Resolution',300)
    close gcf
end
% rastercolor


if NCl ==0
    NClOK=0;
     assemblyortho= cell(0);
     assemblystat= cell(0);
end