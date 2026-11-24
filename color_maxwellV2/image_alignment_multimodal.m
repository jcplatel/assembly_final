function [aligned_red, aligned_green, aligned_blue, calcium_norm, aligned_red890, info] = image_alignment_multimodal(images_norm)
    % Ajout du 11 nov. (ne change RIEN aux 5 premières sorties) : 6e sortie
    % facultative 'info', pour repérer un recalage raté.
    %   info.decalage_1040_px : translation [x y] du bloc 1040 vers le calcium
    %   info.decalage_890_px  : translation [x y] du bloc 890 vers le rouge 1040
    %   info.corr_calcium_avant / _apres : corrélation entre l'image calcium
    %       et le rouge 1040, avant et après recalage (doit monter)
    %   info.corr_rouges_avant / _apres : corrélation entre rouge 890 et
    %       rouge 1040 recalé, avant et après (même fluorophore : doit être
    %       haute après ; c'est le meilleur indicateur)
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

    % 17 nov. (demandé par l'utilisateur) : vérification du VERT. Le vert contient la
    % fuite du GCaMP : il doit tomber pile sur le GCaMP. Sur 446178 14/06 plane0, le
    % vert recalé restait décalé de 3 x 2 px alors que le rouge était juste (vérifié
    % par le contraste des cellules sur les masques) : dans cette session, vert et
    % rouge 1040 ne sont pas superposés dans les images brutes. Si le vert recalé
    % reste décalé d'au moins 2 px avec une correspondance nette, on corrige le VERT
    % SEUL (le rouge, et donc le bloc 890 recalé sur lui, ne changent pas).
    [res_vert, force_vert] = decalage_phase(aligned_green, images_norm.calcium);
    correction_vert = norm(res_vert) >= 2 && force_vert > 0.05;
    tform_vert = tform_1040_to_calcium;
    if correction_vert
        t = translation_de(tform_1040_to_calcium) - res_vert;   % le vert est décalé de res_vert
        tform_vert = affine2d([1 0 0; 0 1 0; t(1) t(2) 1]);
        aligned_green = imwarp(images_norm.green, tform_vert, 'OutputView', ref_view_calcium);
    end
    res_vert_apres = decalage_phase(aligned_green, images_norm.calcium);

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

    if nargout >= 6
        cc = @(a, b) subsref(corrcoef(double(a(:)), double(b(:))), struct('type', '()', 'subs', {{1, 2}}));
        info = struct();
        info.decalage_1040_px   = translation_de(tform_1040_to_calcium);
        info.decalage_890_px    = translation_de(tform_890_to_red1040);
        info.corr_calcium_avant = cc(images_norm.calcium, images_norm.red1040);
        info.corr_calcium_apres = cc(images_norm.calcium, aligned_red);
        info.corr_rouges_avant  = cc(images_norm.red890, aligned_red);
        info.corr_rouges_apres  = cc(aligned_red890, aligned_red);
        info.vert_vs_gcamp_avant_px = res_vert;          % décalage du vert recalé, avant vérification
        info.vert_vs_gcamp_force    = force_vert;        % netteté de la correspondance (pic)
        info.correction_par_le_vert = correction_vert;   % true si le vert a été recalé à part
        info.decalage_vert_px       = translation_de(tform_vert);   % décalage appliqué au vert
        info.vert_vs_gcamp_apres_px = res_vert_apres;
    end

end

function [s, force] = decalage_phase(a, b)
    % décalage entier [x y] de l'image a par rapport à b (corrélation de phase) et
    % netteté du pic (0 = aucune correspondance)
    a = double(a) - mean(double(a(:))); b = double(b) - mean(double(b(:)));
    X = fft2(a) .* conj(fft2(b)); X = X ./ max(abs(X), eps);
    C = real(ifft2(X));
    [force, idx] = max(C(:));
    [dy, dx] = ind2sub(size(C), idx);
    n = size(C);
    s = [mod(dx - 1 + floor(n(2) / 2), n(2)) - floor(n(2) / 2), mod(dy - 1 + floor(n(1) / 2), n(1)) - floor(n(1) / 2)];
end

function t = translation_de(tform)
    % translation [x y] quelle que soit la version de MATLAB (affine2d ou transltform2d)
    if isprop(tform, 'Translation')
        t = tform.Translation;
    else
        t = tform.T(3, 1:2);
    end
end
