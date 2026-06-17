function [R, G, B, id_low_intensity,I] = removal_low_intensities2(R, G, B, seuil_manuel)
    % Calcule l'intensité globale (Norme euclidienne)
    I = sqrt(R.^2 + G.^2 + B.^2);
    
    if nargin < 4 || isempty(seuil_manuel)
        % MODE OTSU AUTOMATIQUE (silencieux)
        % graythresh attend des valeurs entre 0 et 1.
        % On normalise I entre 0 et 1 par rapport à son maximum.
        max_I = max(I(:));
        if max_I > 0
            I_norm = I / max_I;
            thr_otsu_norm = graythresh(I_norm);
            seuil_thr = thr_otsu_norm * max_I; % On ramène à l'échelle d'origine
        else
            seuil_thr = 0; % Cas extrême où tout est à 0
        end
    else
        % MODE MANUEL : On utilise le seuil donné en argument
        seuil_thr = seuil_manuel;
    end
    
    % Filtrage basé sur l'intensité globale
    id_low_intensity = (I < seuil_thr);
    
    nb_low_intensity = sum(id_low_intensity);
    fprintf('Seuil utilisé : %.1f -> Nombre de cellules rejetées : %d / %d\n', ...
            seuil_thr, nb_low_intensity, length(I));
    
    % On met les intensités de ces cellules à NaN pour qu'elles soient ignorées
    R(id_low_intensity) = NaN;
    G(id_low_intensity) = NaN;
    B(id_low_intensity) = NaN;
end