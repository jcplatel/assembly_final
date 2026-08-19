function [R,G,B,Ca,Coomasqx,Coomasqy,Coox,Cooy,outline_gcampx,outline_gcampy] = ROI_extraction2(path, aligned_red, aligned_green, aligned_blue, calcium_norm)

load(fullfile(path, 'Fall.mat'),'stat','iscell');

% Élément structurant principal pour l'érosion
se_erode  = strel('disk', 1);

% Éléments structurants pour le nettoyage/lissage
se_open   = strel('disk', 1);   % enlève petits prolongements
se_close  = strel('disk', 1);   % lisse un peu les contours

[imgH, imgW] = size(aligned_blue);

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
    end
end
end