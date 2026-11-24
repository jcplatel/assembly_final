

% path = "E:\Data\Aurelie\data\cues\all\411582\230331_plane0\";
path = "E:\Data\Aurelie\data\nocues\444175\221125_plane0\";
% out_tif_labels = 'E:\Data\Aurelie\data\chroms\582\230331plane0\mask_suite2p2eroded.tif';
out_tif_labels = 'E:\Data\Aurelie\data\chroms\175\221125plane0\mask_suite2p2eroded.tif';

    load(fullfile(path, 'Fall.mat'), 'stat', 'iscell');

    se = strel('disk', 1);
    [imgH, imgW] = size(aligned_blue);

    Labels = zeros(imgH, imgW, 'uint32');
    MaskBinary = zeros(imgH, imgW, 'uint8');

    centroid_x = [];
    centroid_y = [];
    cell_ids = [];

    m = 0;

    for n = 1:length(stat)
        if iscell(n,1) == 1
            m = m + 1;

            xpix = double(stat{1,n}.xpix) + 1;
            ypix = double(stat{1,n}.ypix) + 1;
                
            ypix = imgH - ypix + 1;

            xpix = max(1, min(xpix, imgW));
            ypix = max(1, min(ypix, imgH));

            lin_idx = sub2ind([imgH, imgW], ypix, xpix);

            mask = false(imgH, imgW);
            mask(lin_idx) = true;

            mask_eroded = imerode(mask, se);
            % mask_eroded = mask;
            eroded_idx = find(mask_eroded);

            if isempty(eroded_idx)
                eroded_idx = lin_idx;
                mask_eroded = false(imgH, imgW);
                mask_eroded(eroded_idx) = true;
            end

            % Image de labels napari
            Labels(eroded_idx) = uint32(m);

            % Masque binaire pour image annotée
            MaskBinary(eroded_idx) = 255;

            % Centroïde pour positionner le texte
            props = regionprops(mask_eroded, 'Centroid');
            centroid_x(m,1) = props(1).Centroid(1);
            centroid_y(m,1) = props(1).Centroid(2);
            cell_ids(m,1) = m;
        end
    end

    fprintf('Nombre de cellules labellisées : %d\n', m);

    %% Sauvegarde du vrai fichier labels pour napari
    

    if isfile(out_tif_labels)
        delete(out_tif_labels);
    end

    t = Tiff(out_tif_labels, 'w');
    tagstruct.ImageLength = imgH;
    tagstruct.ImageWidth = imgW;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 32;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Compression = Tiff.Compression.None;

    t.setTag(tagstruct);
    t.write(Labels);
    t.close();

    fprintf('Image de labels sauvée : %s\n', out_tif_labels);

    %% Création de l'image annotée avec numéros en noir
    % insertText retourne une image RGB
    positions = [centroid_x+5, centroid_y];
    text_str = arrayfun(@num2str, cell_ids, 'UniformOutput', false);

    AnnotRGB = insertText( ...
        MaskBinary, ...
        positions, ...
        text_str, ...
        'FontSize', 8, ...
        'TextColor', 'white', ...
        'BoxOpacity', 0, ...
        'AnchorPoint', 'LeftCenter');
out_tif_annot = 'E:\Data\Aurelie\data\chroms\582\230331plane0\mask_suite2p_annotated2.tif';
    if isfile(out_tif_annot)
        delete(out_tif_annot);
    end

    imwrite(AnnotRGB, out_tif_annot);
    fprintf('Image annotée sauvée : %s\n', out_tif_annot);





