function [aligned_red, aligned_green, aligned_blue, calcium_norm, aligned_red890] = image_alignment_multimodal(images_norm)
    % =========================================================================
    % 1. PRÉPARATION DE L'OPTIMISEUR MULTI-MODAL
    % =========================================================================
    [optimizer, metric] = imregconfig('multimodal');
    % On réduit le pas initial pour éviter les divergences
    optimizer.InitialRadius = optimizer.InitialRadius / 5;
    
    % Vue de référence absolue (le calcium fixe)
    ref_view_calcium = imref2d(size(images_norm.calcium));

    % =========================================================================
    % 2. ALIGNEMENT DU BLOC 1040 (Rouge et Vert) SUR LE CALCIUM
    % =========================================================================
    % On calcule la transformation en utilisant le red1040 (souvent plus propre) 
    % par rapport au calcium (fixe).
    tform_1040_to_calcium = imregtform(images_norm.red1040, images_norm.calcium, ...
                                       'translation', optimizer, metric);
    
    % On applique cette transformation au Rouge 1040
    aligned_red = imwarp(images_norm.red1040, tform_1040_to_calcium, ...
                         'OutputView', ref_view_calcium);
                     
    % On applique EXACTEMENT la même transformation au Vert 
    % (car l'image verte est trop bruitée pour guider son propre alignement)
    aligned_green = imwarp(images_norm.green, tform_1040_to_calcium, ...
                           'OutputView', ref_view_calcium);

    % =========================================================================
    % 3. ALIGNEMENT DU BLOC 890 (Rouge 890 et Bleu) SUR LE ROUGE 1040 ALIGNÉ
    % =========================================================================
    % Ici, le red1040_aligné devient la nouvelle image de référence (fixe).
    % Le red890 est l'image qui bouge pour rattraper le red1040.
    tform_890_to_red1040 = imregtform(images_norm.red890, aligned_red, ...
                                      'translation', optimizer, metric);
                                  
    % On applique la transformation au Rouge 890
    aligned_red890 = imwarp(images_norm.red890, tform_890_to_red1040, ...
                            'OutputView', ref_view_calcium);
                        
    % On applique la même transformation au Bleu (acquis en même temps que le R890)
    aligned_blue = imwarp(images_norm.blue, tform_890_to_red1040, ...
                          'OutputView', ref_view_calcium);

    % =========================================================================
    % 4. FINALISATION ET NETTOYAGE
    % =========================================================================
    calcium_norm = images_norm.calcium; % Le calcium n'a pas bougé
    
    % imwarp peut parfois créer de très légères valeurs négatives lors des 
    % interpolations aux bords. On ramène tout ce qui est sous 0 à 0.
    aligned_red(aligned_red < 0) = 0;
    aligned_green(aligned_green < 0) = 0;
    aligned_blue(aligned_blue < 0) = 0;
    aligned_red890(aligned_red890 < 0) = 0;

end