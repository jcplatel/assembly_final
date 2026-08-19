% Script pour générer un schéma expérimental modifiable dans PowerPoint
% Nécessite Windows et Microsoft PowerPoint installé.

disp('Ouverture de PowerPoint en cours...');

% 1. Initialisation de l'application PowerPoint via COM
ppt = actxserver('PowerPoint.Application');
ppt.Visible = 1; % Rendre la fenêtre visible

% Création d'une nouvelle présentation
pres = ppt.Presentations.Add();

% Ajout d'une diapositive vide (le layout 12 correspond généralement à une diapo vide)
slide = pres.Slides.Add(1, 12); 

% Constantes pour les couleurs (Format RGB pour COM : R + G*256 + B*65536)
colorBlack = 0;
colorGray = 200 + 200*256 + 200*65536;

% ---------------------------------------------------------
% 2. CRÉATION DE LA FRISE CHRONOLOGIQUE
% ---------------------------------------------------------

% Ligne principale (X_début, Y_début, X_fin, Y_fin)
timeline = slide.Shapes.AddLine(150, 150, 800, 150);
timeline.Line.Weight = 2.5;
timeline.Line.ForeColor.RGB = colorBlack;

% --- Point P0 ---
% Marqueur (Cercle)
slide.Shapes.AddShape(9, 150-6, 150-6, 12, 12).Fill.ForeColor.RGB = colorBlack; % 9 = msoShapeOval
% Titre du point
t_P0 = slide.Shapes.AddTextbox(1, 130, 110, 100, 30); % 1 = msoTextOrientationHorizontal
t_P0.TextFrame.TextRange.Text = 'P0';
t_P0.TextFrame.TextRange.Font.Bold = 1;
% Description
d_P0 = slide.Shapes.AddTextbox(1, 100, 165, 150, 50);
d_P0.TextFrame.TextRange.Text = 'AAV1.hSyn.GCaMP8m\nInjection ICV';

% --- Point P45 ---
slide.Shapes.AddShape(9, 350-6, 150-6, 12, 12).Fill.ForeColor.RGB = colorBlack;
t_P45 = slide.Shapes.AddTextbox(1, 330, 110, 100, 30);
t_P45.TextFrame.TextRange.Text = 'P45';
t_P45.TextFrame.TextRange.Font.Bold = 1;
d_P45 = slide.Shapes.AddTextbox(1, 300, 165, 150, 50);
d_P45.TextFrame.TextRange.Text = 'Chirurgie crânienne\n(Fenêtre optique)';

% --- Point P53-57 ---
slide.Shapes.AddShape(9, 550-6, 150-6, 12, 12).Fill.ForeColor.RGB = colorBlack;
t_P53 = slide.Shapes.AddTextbox(1, 530, 110, 100, 30);
t_P53.TextFrame.TextRange.Text = 'P53-57';
t_P53.TextFrame.TextRange.Font.Bold = 1;
d_P53 = slide.Shapes.AddTextbox(1, 520, 165, 150, 50);
d_P53.TextFrame.TextRange.Text = 'Entraînement\nsur tapis roulant';

% --- Point P60-P64 ---
slide.Shapes.AddShape(9, 750-6, 150-6, 12, 12).Fill.ForeColor.RGB = colorBlack;
t_P60 = slide.Shapes.AddTextbox(1, 730, 110, 100, 30);
t_P60.TextFrame.TextRange.Text = 'P60-P64';
t_P60.TextFrame.TextRange.Font.Bold = 1;
d_P60 = slide.Shapes.AddTextbox(1, 710, 165, 150, 50);
d_P60.TextFrame.TextRange.Text = 'Imagerie Calcium 2P\nlongitudinale';

% ---------------------------------------------------------
% 3. AJOUT DES ENCADRÉS DE TEXTE (Métadonnées)
% ---------------------------------------------------------

% Encadré pour les Génotypes (En haut à gauche)
genoBox = slide.Shapes.AddShape(1, 20, 20, 150, 70); % 1 = msoShapeRectangle
genoBox.Fill.Visible = 0; % Fond transparent
genoBox.Line.ForeColor.RGB = colorBlack;
genoText = slide.Shapes.AddTextbox(1, 25, 25, 140, 60);
genoText.TextFrame.TextRange.Text = 'Genotypes:\n- Tcyt5 +/-\n- Tcyt5 +/- ; Emx1Cre +/-';

% ---------------------------------------------------------
% 4. AJOUT DES ZONES RÉSERVÉES POUR VOS IMAGES
% ---------------------------------------------------------

% Zone pour le setup comportemental (tapis roulant)
img_setup = slide.Shapes.AddShape(1, 100, 250, 450, 200);
img_setup.Fill.ForeColor.RGB = colorGray;
img_setup.Line.Visible = 0;
txt_setup = slide.Shapes.AddTextbox(1, 150, 320, 350, 50);
txt_setup.TextFrame.TextRange.Text = '[Glissez-déposez les illustrations du tapis et de la souris ici]';
txt_setup.TextFrame.TextRange.Font.Color.RGB = colorBlack;

% Zone pour la microscopie et l'hippocampe
img_micro = slide.Shapes.AddShape(1, 600, 250, 320, 250);
img_micro.Fill.ForeColor.RGB = colorGray;
img_micro.Line.Visible = 0;
txt_micro = slide.Shapes.AddTextbox(1, 620, 340, 280, 50);
txt_micro.TextFrame.TextRange.Text = '[Glissez-déposez l''image de l''hippocampe et l''objectif ici]';
txt_micro.TextFrame.TextRange.Font.Color.RGB = colorBlack;

disp('Diapositive générée avec succès ! Vous pouvez maintenant modifier le texte et insérer vos images directement dans PowerPoint.');