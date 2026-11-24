function slides_miniatures(sortie, par_slide)
% Miniatures des images compositecalcium (R,G,B) des sessions incluses
% (inclure = 1 dans sessions_brainbow.xlsx), réparties sur plusieurs slides.
%   slides_miniatures('E:\...\composites_miniatures')        16 images par slide
%   slides_miniatures('E:\...\composites_miniatures', 18)
if nargin < 2, par_slide = 16; end
addpath(genpath('C:\Users\jcplatel\Desktop\matlab\assembly_final\pipeline_brainbow'));
config = config_brainbow();
config.fichier_sessions = 'E:\Data\Aurelie\analysis\Final\assembly\nocues\sessions_brainbow.xlsx';
S = lire_table_sessions(config);
n = height(S);
if ~isfolder(sortie), mkdir(sortie); end
png = cell(n, 1); etiquette = cell(n, 1);
for s = 1:n
    path = avec_separateur(S.dossier_suite2p{s});
    etiquette{s} = strrep(regexprep(S.nom{s}, '_sans_indice$', ''), '_', ' ');
    fc = trouver_composite(path);
    if isempty(fc)
        fprintf(2, '[%s] aucune image composite\n', S.nom{s});
        continue
    end
    [~, nf, ef] = fileparts(fc);
    if ~strcmpi([nf ef], 'compositecalcium.tif'), etiquette{s} = [etiquette{s} ' (' nf ef ')']; end
    info = imfinfo(fc); ch = {};
    for pg = 1:numel(info)
        A = double(imread(fc, pg));
        for c = 1:size(A, 3), ch{end+1} = A(:, :, c); end %#ok<AGROW>
    end
    I = zeros([size(ch{1}) 3]);
    for c = 1:min(3, numel(ch))       % rouge, vert, bleu (le 4e canal = calcium, non affiché)
        v = ch{c}; p = prctile(v(:), [1 99.7]);
        I(:, :, c) = min(max((v - p(1)) / max(p(2) - p(1), eps), 0), 1);
    end
    png{s} = fullfile(sortie, sprintf('%02d_%s.png', s, S.nom{s}));
    imwrite(I, png{s});
    fprintf('[%s] %s\n', S.nom{s}, [nf ef]);
end

% Slides
import mlreportgen.ppt.*
fp = fullfile(sortie, 'composites_miniatures.pptx');
ppt = Presentation(fp); open(ppt);
W = 13.333; H = 7.5; haut = 0.5; lab = 0.22;
nc = ceil(sqrt(par_slide * W / (H - haut))); nr = ceil(par_slide / nc);
cw = (W - 0.3) / nc; rh = (H - haut - 0.1) / nr;
t = min(cw - 0.08, rh - lab - 0.05);              % côté de l'image (carrée)
nslides = ceil(n / par_slide);
for k = 1:nslides
    sl = add(ppt, 'Blank');
    idx = (k-1)*par_slide + 1 : min(k*par_slide, n);
    tb = TextBox(); tb.X = '0.2in'; tb.Y = '0.05in'; tb.Width = '12.9in'; tb.Height = '0.4in';
    p = Paragraph(sprintf('Composite images (compositecalcium, R/G/B) — sessions %d–%d / %d', idx(1), idx(end), n));
    p.FontSize = '16pt'; p.Bold = true; add(tb, p); add(sl, tb);
    for j = 1:numel(idx)
        s = idx(j);
        r = floor((j-1) / nc); c = mod(j-1, nc);
        x0 = 0.15 + c*cw; y0 = haut + r*rh;
        L = TextBox(); L.X = sprintf('%.3fin', x0); L.Y = sprintf('%.3fin', y0);
        L.Width = sprintf('%.3fin', cw); L.Height = sprintf('%.3fin', lab);
        q = Paragraph(etiquette{s}); q.FontSize = '9pt'; q.HAlign = 'center'; add(L, q); add(sl, L);
        if isempty(png{s}), continue; end
        pic = Picture(png{s});
        pic.X = sprintf('%.3fin', x0 + (cw - t)/2); pic.Y = sprintf('%.3fin', y0 + lab);
        pic.Width = sprintf('%.3fin', t); pic.Height = sprintf('%.3fin', t);
        add(sl, pic);
    end
end
close(ppt);

% Retire les champs date / pied de page / numéro du modèle
try
    app = actxserver('PowerPoint.Application');
    pr = app.Presentations.Open(fp, 0, 0, 0);
    for k = 1:pr.Slides.Count
        sh = pr.Slides.Item(k).Shapes;
        for i = sh.Count:-1:1
            if sh.Item(i).Type == 14, sh.Item(i).Delete; end
        end
    end
    pr.Save; pr.Close; app.Quit;
catch err
    warning('Champs du modèle non retirés : %s', err.message);
end
fprintf('%d slides : %s\n', nslides, fp);
end
