function Raster = transient_detection_std (Tr1b,MinPeakDistancesce,MinPeakDistance,threshold_peak,WinRest,WinActive);

[NCell, Nz] = size(Tr1b);
Raster = zeros(NCell, Nz);
Acttmp2 = cell(1, NCell);
ampli = cell(1, NCell);
th = zeros(1, NCell);

for i = 1:NCell
%%%%%%%%%%%%%%%%%%%%
    % trace=Tr1b(i,:);
    trace=Tr1b(i,WinRest);
 % 1. Première estimation grossière du bruit (MAD classique)
    sigma_rough = 1.4826 * mad(trace, 1);
    med_val = median(trace);
    
    % 2. Identifier les zones "calmes" (tout ce qui est < Médiane + 2*sigma)
    %    On est large volontairement pour être sûr d'exclure tout début de pic
    mask_calme = trace < (med_val + 2 * sigma_rough);
    
    % 3. Calcul du "VRAI" bruit uniquement sur les zones calmes
    if sum(mask_calme) > 100 % Sécurité : s'il reste assez de points
        trace_calme = trace(mask_calme);
        % On recalcule le bruit sur cette trace épurée
        % On utilise l'écart-type ici car on a enlevé les outliers
        sigma_reelle = std(trace_calme); 
        baseline_reelle = median(trace_calme);
    else
        % Fallback si tout le signal est considéré comme actif (rare)
        sigma_reelle = sigma_rough;
        baseline_reelle = med_val;
    end
    
    % 4. Définition du seuil final (plus strict maintenant qu'on a le vrai bruit)
    seuil_final = baseline_reelle + 3 * sigma_reelle; % 3.5 ou 4 selon sensibilité voulue

    % th(i) = threshold_peak * mad_trace * 1.4826;
    seuil(i) = seuil_final;
    % th(i) = threshold_peak * sigma_bruit
    
    % Détecter les pics sur la trace complète
    [amplitude, locs] = findpeaks(Tr1b(i, :), 'MinPeakProminence', seuil(i), 'MinPeakDistance', MinPeakDistance, 'MinPeakWidth',5);
    
    % Filtrer : enlever les pics qui tombent sur des bad frames
    % % % valid_mask = ~bad_frames(locs)';
    % % % locs = locs(valid_mask);
    % % % amplitude = amplitude(valid_mask);

    % Filtrer les pics en mouvement
    valeurs_identiques = intersect(locs, WinActive);
    [locs_sans_ide, idx] = setdiff(locs(:), valeurs_identiques);
    ampli_sans_ide = amplitude(idx);
    Acttmp2{i} = locs_sans_ide;

end

% Créer Raster avec les transients valides
for i = 1:NCell
    if ~isempty(Acttmp2{i})
        Raster(i, Acttmp2{i}) = 1;
    end
end