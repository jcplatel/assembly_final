function [R,G,B,Ca,Coomasqx,Coomasqy,Coox,Cooy,outline_gcampx,outline_gcampy,extra] = ROI_extraction2(path, aligned_red, aligned_green, aligned_blue, calcium_norm)
% Ajout du 11 nov. (ne change RIEN aux 10 premières sorties) : 11e sortie
% facultative 'extra', pour vérifier les couleurs des cellules faibles.
%   extra.R_fond, G_fond, B_fond, Ca_fond : médiane de chaque canal dans un
%       anneau autour de la cellule (de 3 à 8 pixels du masque), sans les
%       pixels d'aucune ROI de suite2p (cellules ou non)
%   extra.R_local, G_local, B_local : intensité de la cellule moins son
%       fond local (coupée à 0)
%   extra.n_pixels : pixels du masque érodé ; extra.n_pixels_fond : de l'anneau
%   extra.roi_suite2p : numéro de ROI de suite2p (à partir de 0) de chaque cellule

% 18 nov. : masques_corriges.mat s'il existe (charger_masques), sinon Fall.mat.
% L'anneau de fond (sortie extra) exclut toujours les masques d'origine.
[stat, iscell] = charger_masques(path);
stat_origine = load(fullfile(path, 'Fall.mat'), 'stat'); stat_origine = stat_origine.stat;

% Élément structurant principal pour l'érosion
se_erode  = strel('disk', 1);

% Éléments structurants pour le nettoyage/lissage
se_open   = strel('disk', 1);   % enlève petits prolongements
se_close  = strel('disk', 1);   % lisse un peu les contours

[imgH, imgW] = size(aligned_blue);

% Toutes les ROI de suite2p (pour les retirer des anneaux de fond local)
faire_extra = nargout >= 11;
if faire_extra
    toutes_roi = false(imgH, imgW);
    for n = 1:length(stat_origine)
        xa = max(1, min(double(stat_origine{1,n}.xpix) + 1, imgW));
        ya = max(1, min(double(stat_origine{1,n}.ypix) + 1, imgH));
        toutes_roi(sub2ind([imgH, imgW], ya, xa)) = true;
    end
    se_in  = strel('disk', 3);
    se_out = strel('disk', 8);
end

NCell_total = sum(iscell(:,1) > 0);
B  = zeros(NCell_total, 1);
G  = zeros(NCell_total, 1);
R  = zeros(NCell_total, 1);
Ca = zeros(NCell_total, 1);

outline_gcampx = cell(1, NCell_total);
outline_gcampy = cell(1, NCell_total);
neuropil       = cell(1, NCell_total);
Coomasqx       = cell(1, NCell_total);
Coomasqy       = cell(1, NCell_total);
Coox           = cell(1, NCell_total);
Cooy           = cell(1, NCell_total);
extra = struct('R_fond', nan(NCell_total, 1), 'G_fond', nan(NCell_total, 1), 'B_fond', nan(NCell_total, 1), ...
    'Ca_fond', nan(NCell_total, 1), 'n_pixels', zeros(NCell_total, 1), 'n_pixels_fond', zeros(NCell_total, 1), ...
    'roi_suite2p', zeros(NCell_total, 1));

m = 0;

for n = 1:length(stat)
    if iscell(n,1) == 1
        m = m + 1;

        outline_gcampx{m} = double(stat{1,n}.xext) + 1;
        outline_gcampy{m} = double(stat{1,n}.yext) + 1;
        Coox{m} = double(stat{1,n}.xpix) + 1;
        Cooy{m} = double(stat{1,n}.ypix) + 1;
        neuropil{m} = stat{1,n}.neuropil_mask;

        xpix = double(stat{1,n}.xpix) + 1;
        ypix = double(stat{1,n}.ypix) + 1;

        xpix = max(1, min(xpix, imgW));
        ypix = max(1, min(ypix, imgH));

        lin_idx = sub2ind([imgH, imgW], ypix, xpix);

        mask = false(imgH, imgW);
        mask(lin_idx) = true;

        % 1. Érosion initiale
        mask_clean = imerode(mask, se_erode);

        % 2. Si plusieurs objets apparaissent, garder uniquement le plus gros
        if any(mask_clean(:))
            mask_clean = bwareafilt(mask_clean, 1);
        end

        % 3. Lisser la forme et enlever les petites pointes/prolongements
        if any(mask_clean(:))
            mask_clean = imopen(mask_clean, se_open);
            mask_clean = imclose(mask_clean, se_close);
        end

        % 4. Sécurité : si trop agressif, revenir au masque original
        if ~any(mask_clean(:))
            mask_clean = mask;
        end

        final_idx = find(mask_clean);

        [new_y, new_x] = ind2sub([imgH, imgW], final_idx);
        Coomasqx{m} = new_x;
        Coomasqy{m} = new_y;

        B(m)  = median(aligned_blue(final_idx));
        G(m)  = median(aligned_green(final_idx));
        R(m)  = median(aligned_red(final_idx));
        Ca(m) = median(calcium_norm(final_idx));

        % Fond local (ajout du 11 nov., sortie 'extra' seulement)
        if faire_extra
            extra.n_pixels(m) = numel(final_idx);
            extra.roi_suite2p(m) = n - 1;
            anneau = imdilate(mask, se_out) & ~imdilate(mask, se_in) & ~toutes_roi;
            extra.n_pixels_fond(m) = nnz(anneau);
            if extra.n_pixels_fond(m) >= 20
                extra.R_fond(m)  = median(aligned_red(anneau));
                extra.G_fond(m)  = median(aligned_green(anneau));
                extra.B_fond(m)  = median(aligned_blue(anneau));
                extra.Ca_fond(m) = median(calcium_norm(anneau));
            end
        end
    end
end
if faire_extra
    extra.R_local = max(0, R - extra.R_fond);
    extra.G_local = max(0, G - extra.G_fond);
    extra.B_local = max(0, B - extra.B_fond);
end
end