function [MdlLin,  MdlSVM, MdlRF,X,Ycat, X_raw,Ycat_raw, cell_IDs] = prediction_color(colorcell,roi_red,roi_green,roi_blue,roi_ca,roi_green_corrected)

if nargin==4
    X = [roi_red, roi_green,roi_blue];
elseif nargin==5
    X = [roi_red, roi_green,roi_blue,roi_ca];
elseif nargin==6
    X = [roi_red, roi_green,roi_blue,roi_green_corrected,roi_ca];
end

I = sqrt(roi_red.^2 + roi_green.^2 + roi_blue.^2);

% rouge_brut_z = zscore(roi_red_raw);
% vert_brut_z  = zscore(roi_green_raw);
% bleu_brut_z  = zscore(roi_blue_raw);
% X = [I_blue, I_red, I_green, I_calcium, I_green_dec];
ratio_RG = roi_red   ./ (roi_green + 0.01);
ratio_RB  = roi_red   ./ (roi_blue + 0.01);
ratio_GB  = roi_green ./ (roi_blue + 0.01);
X = [X,I,ratio_RG,ratio_RB,ratio_GB];

X_raw=X;
% Supprimer la catégorie 'noir' devenue vide (important pour confusionchart)

% --- Partition sur les données filtrées ---

rng(1);
K = 10;

Y_num = str2double(string(colorcell));

% On recrée la variable catégorielle proprement avec les 8 couleurs
% classNames = {'rouge','vert','bleu','jaune','magenta','cyan','blanc','noir'};
classNames = {'red','green','blue','yellow','magenta','cyan','white','black'};
Ycat = categorical(Y_num, 1:8, classNames);
Ycat_raw=Ycat;
cell_IDs_original = (1:length(Ycat))';

maskNoNoir = Ycat ~= 'black';   % vecteur logique : true pour toutes les classes sauf noir
% mask_noir   = (Ycat == 'noir');
% mask_actifs = ~mask_noir;   % toutes les classes sauf noir
X_filt   = X(maskNoNoir, :);   % features filtrées
X = X_filt;

Ycat = Ycat(maskNoNoir);  % labels filtrés
cell_IDs = cell_IDs_original(maskNoNoir);
% cell_IDs = cell_IDs_original;
mask_valide = ~any(isnan(X), 2) & ~any(isinf(X), 2);
% 
% On ne garde que les cellules valides
X = X(mask_valide, :);
% X = zscore(X);
Ycat = Ycat(mask_valide);
cell_IDs = cell_IDs(mask_valide);

% Création de la partition stratifiée
cvp = cvpartition(Ycat, 'KFold', K);


%% ================= 1. MODÈLE LINÉAIRE (LDA) =================
CVMdlLin = fitcdiscr(X, Ycat, 'DiscrimType', 'linear', 'CVPartition', cvp);

accLin = 1 - kfoldLoss(CVMdlLin);
YPredLin = kfoldPredict(CVMdlLin);
% Accuracy sans les noirs


subplot(2, 2, 1);
confusionchart(Ycat, YPredLin, 'Title', sprintf('1. Linéaire (LDA)\nAccuracy: %.1f%%', accLin*100));

% %% ================= 2. MODÈLE QUADRATIQUE (QDA) =================
% % Attention: La QDA peut planter si une classe a trop peu de points par rapport
% % au nombre de variables. MATLAB mettra un 'pseudoLinear' si besoin.
% CVMdlQuad = fitcdiscr(X, Ycat, 'DiscrimType', 'quadratic', 'CVPartition', cvp);
% accQuad = 1 - kfoldLoss(CVMdlQuad);
% YPredQuad = kfoldPredict(CVMdlQuad);
% % Accuracy sans les noirs

% 
% subplot(2, 2, 2);
% confusionchart(Ycat, YPredQuad, 'Title', sprintf('2. Quadratique (QDA)\nAccuracy: %.1f%%', accQuad*100));

%% ================= 3. MODÈLE SVM (Support Vector Machine) =================
% On utilise ECOC (Error-Correcting Output Codes) pour gérer le multiclass
% avec un noyau linéaire pour rester robuste.
tSVM = templateSVM('KernelFunction', 'linear', 'Standardize', true);
CVMdlSVM = fitcecoc(X, Ycat, 'Learners', tSVM, 'CVPartition', cvp);

accSVM = 1 - kfoldLoss(CVMdlSVM);
YPredSVM = kfoldPredict(CVMdlSVM);

subplot(2, 2, 2);
confusionchart(Ycat, YPredSVM, 'Title', sprintf('3. SVM (Linéaire)\nAccuracy: %.1f%%', accSVM*100));

%% ================= 4. MODÈLE RANDOM FOREST (Arbres) =================
% On utilise 'Bag' (Bootstrap Aggregation) avec 100 arbres de décision.
tTrees = templateTree('Reproducible', true);
CVMdlRF = fitcensemble(X, Ycat, 'Method', 'Bag', 'NumLearningCycles', 100, ...
    'Learners', tTrees, 'CVPartition', cvp);

accRF = 1 - kfoldLoss(CVMdlRF);
YPredRF = kfoldPredict(CVMdlRF);

subplot(2, 2, 3);
confusionchart(Ycat, YPredRF, 'Title', sprintf('4. Random Forest\nAccuracy: %.1f%%', accRF*100));

%% --- RÉCAPITULATIF CONSOLE ---
fprintf('\n====== RÉSULTATS DE LA CROSS-VALIDATION ======\n');
fprintf('1. Linéaire (LDA)       : %.2f %%\n', accLin * 100);
% fprintf('2. Quadratique (QDA)    : %.2f %%\n', accQuad * 100);
fprintf('3. SVM Multi-classes    : %.2f %%\n', accSVM * 100);
fprintf('4. Random Forest        : %.2f %%\n', accRF * 100);

% Entraînement sur 100% des données
%% ================= 1. MODÈLE LINÉAIRE (LDA) =================

MdlLin = fitcdiscr(X, Ycat, 'DiscrimType', 'linear');

%% ================= 2. MODÈLE QUADRATIQUE (QDA) =================

% MdlQuad = fitcdiscr(X, Ycat, 'DiscrimType', 'quadratic');

%% ================= 3. MODÈLE SVM (Support Vector Machine) =================
% Entraînement sur 100% des données avec standardisation
tSVM = templateSVM('KernelFunction', 'linear', 'Standardize', true);
MdlSVM = fitcecoc(X, Ycat, 'Learners', tSVM);

%% ================= 4. MODÈLE RANDOM FOREST (Arbres) =================
% Entraînement sur 100% des données avec 100 arbres
tTrees = templateTree('Reproducible', true);
MdlRF = fitcensemble(X, Ycat, 'Method', 'Bag', 'NumLearningCycles', 100, 'Learners', tTrees);

%% --- SAUVEGARDE DES MODÈLES FINAUX ---
% Enregistre les 4 modèles entraînés dans un fichier nommé "Modeles_Brainbow.mat"
% save('E:\Data\Aurelie\data\chroms\119\220919\Modeles_Brainbow.mat', 'MdlLin', 'MdlQuad', 'MdlSVM', 'MdlRF');

end