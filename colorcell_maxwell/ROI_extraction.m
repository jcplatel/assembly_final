function [R,G,B,Ca,Coomasqx,Coomasqy,Coox,Cooy] = ROI_extraction  (path,aligned_red,aligned_green,aligned_blue,calcium_norm) 

load(fullfile(path, 'Fall.mat'),'stat','iscell');
% L'élément structurant pour l'érosion du masque de la cellule
se = strel('disk', 1);
[imgH,imgW ]= size (aligned_blue);

NCell_total = sum(iscell(:,1) > 0);
B  = zeros(NCell_total, 1);
G = zeros(NCell_total, 1);
R   = zeros(NCell_total, 1);
Ca = zeros(NCell_total, 1);

outline_gcampx = cell(1, NCell_total);
outline_gcampy = cell(1, NCell_total);
neuropil = cell(1, NCell_total);
Coomasqx = cell(1, NCell_total);
Coomasqy = cell(1, NCell_total);

m = 0;

for n = 1 : length(stat)
    if iscell(n,1) == 1
        m = m + 1;
        % Extraction des données brutes (Suite2P est en Python 0-based -> MATLAB 1-based)
        outline_gcampx{m} = double(stat{1, n}.xext) + 1;
        outline_gcampy{m} = double(stat{1, n}.yext) + 1;
        Coox(m)={double(stat{1, n}.xpix+1)};%+1  to adpat suite2p output   ONLY 1 ??????
        Cooy(m)={double(stat{1, n}.ypix+1)};

        neuropil{m} = stat{1, n}.neuropil_mask;
        
        xpix = double(stat{1, n}.xpix) + 1;
        ypix = double(stat{1, n}.ypix) + 1;
        
        % Sécurité : s'assurer qu'aucun pixel ne sort de l'image
        xpix = max(1, min(xpix, imgW));
        ypix = max(1, min(ypix, imgH));
        
        % --- ÉROSION CORRIGÉE ---
        % A. Convertir les coordonnées (x,y) en un index linéaire unique
        lin_idx = sub2ind([imgH, imgW], ypix, xpix);
        
        % B. Créer un masque vide, l'allumer aux coordonnées de la cellule, puis l'éroder
        mask = false(imgH, imgW);
        mask(lin_idx) = true;
        mask_eroded = imerode(mask, se);
        % mask_dilated = imdilate(mask_eroded, se);
        
        
        % C. Récupérer les nouveaux index après érosion
        eroded_idx = find(mask_eroded);
        
        % D. Sécurité : si la cellule était trop petite et a été entièrement effacée
        if isempty(eroded_idx)
            eroded_idx = lin_idx; % On garde l'originale
        end
        
        % Sauvegarde des nouvelles coordonnées X,Y érodées
        [new_y, new_x] = ind2sub([imgH, imgW], eroded_idx);
        Coomasqx{m} = new_x;
        Coomasqy{m} = new_y;
        
        % % --- NOUVEAUTÉ : "CORE PIXELS" ET CALCUL DE LA MÉDIANE ---
        % % On récupère les pixels sur l'image LISSÉE
        % pixels_B = blue_smooth(eroded_idx);
        % pixels_G = green_smooth(eroded_idx);
        % pixels_R = red_smooth(eroded_idx);
        % 
        % % Filtrage des 25% de pixels les plus sombres (bruit de fond restant)
        % if length(eroded_idx) >= 4 % S'il y a assez de pixels
        %     intensite_totale = pixels_B + pixels_G + pixels_R;
        %     seuil_bas = prctile(intensite_totale, 25);
        %     idx_core = (intensite_totale >= seuil_bas);
        % 
        %     % Calcul de la MEDIANE sur le "cœur" lumineux de la cellule
        %     meanBlue(m)  = median(pixels_B(idx_core));
        %     meanGreen(m) = median(pixels_G(idx_core));
        %     meanRed(m)   = median(pixels_R(idx_core));
        % else
        %     % Cellule minuscule : on prend juste la médiane globale
        %     meanBlue(m)  = median(pixels_B);
        %     meanGreen(m) = median(pixels_G);
        %     meanRed(m)   = median(pixels_R);
        % end
        B(m)  = median(aligned_blue(eroded_idx));
        G(m) = median(aligned_green(eroded_idx));
        R(m)   = median(aligned_red(eroded_idx));
        Ca(m)   = median(calcium_norm(eroded_idx));
        % B(m)  = mean(aligned_blue(eroded_idx));
        % G(m) = mean(aligned_green(eroded_idx));
        % R(m)   = mean(aligned_red(eroded_idx));
        % Ca(m)   = mean(calcium_norm(eroded_idx));
    end
end

NCell = m;
% fprintf('Extraction terminée pour %d cellules.\n', NCell);

% R = meanRed(1:NCell); % S'assurer de ne prendre que les cellules valides
% G = meanGreen(1:NCell);
% B = meanBlue(1:NCell);
% Ca = meanCalcium (1:NCell);

% if nargin==3
%     YFP_corrige = G - (k * Ca);
% end
% 
% cutoff = 0.05 * prctile(YFP_corrige, 99);
% YFP_corrige(YFP_corrige < cutoff) = 0;

end
