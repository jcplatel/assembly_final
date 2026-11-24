function [accuracy_active, accuracy_ref, C] = confusion_matrix(colorcell_ref, colorcell_test, nom_ref, nom_test)
% CONFUSION_MATRIX  Compare deux classifications de couleurs (colorcell).
%
%   [accuracy_active, accuracy_ref, C] = confusion_matrix(colorcell_ref, colorcell_test, nom_ref, nom_test)
%
%   colorcell_ref  : classification de référence (ex. in vivo)  -> LIGNES, axe Y
%   colorcell_test : classification comparée (ex. Biapy)         -> COLONNES, axe X
%   nom_ref, nom_test : noms affichés sur les axes (facultatifs)
%
%   C(i,j) = nombre de neurones classés i par la référence et j par la méthode comparée.
%
%   Sorties :
%   accuracy_active : accord sur les neurones colorés selon LES DEUX méthodes
%                     (symétrique : ne change pas si on inverse les arguments)
%   accuracy_ref    : accord sur les neurones colorés selon la RÉFÉRENCE ;
%                     un neurone noir pour la méthode comparée compte comme erreur
%                     (dépend du sens)
%   C               : matrice de confusion 8 x 8
%
%   Exemple : confusion_matrix(colorcell_invivo, colorcell_biapy, 'in vivo', 'Biapy')

if nargin < 3 || isempty(nom_ref),  nom_ref  = 'Référence'; end
if nargin < 4 || isempty(nom_test), nom_test = 'Méthode comparée'; end

colorcell_ref  = colorcell_ref(:);
colorcell_test = colorcell_test(:);
if numel(colorcell_ref) ~= numel(colorcell_test)
    error('Les deux colorcell n''ont pas le même nombre de cellules (%d contre %d).', ...
        numel(colorcell_ref), numel(colorcell_test));
end

% --- Labels ---
labels_str  = {'Red', 'Green', 'Blue', 'Yellow', 'Magenta', 'Cyan', 'White', 'Black'};
all_classes = 1:8;
idx_black   = 8;

% --- Matrice de confusion : lignes = référence, colonnes = méthode comparée ---
C = confusionmat(colorcell_ref, colorcell_test, 'Order', all_classes);

% --- Accord global (tous les neurones, noirs compris) ---
accuracy = trace(C) / sum(C(:)) * 100;
fprintf('Accord global : %.1f%% (%d neurones)\n', accuracy, sum(C(:)));

% --- Accord sur les neurones colorés selon les deux méthodes (symétrique) ---
idxs_to_keep = setdiff(all_classes, idx_black);
C_active = C(idxs_to_keep, idxs_to_keep);
accuracy_active = trace(C_active) / sum(C_active(:)) * 100;
fprintf('Accord sur les neurones colorés par les deux méthodes : %.1f%% (%d neurones)\n', ...
    accuracy_active, sum(C_active(:)));

% --- Accord sur les neurones colorés selon la référence (dépend du sens) ---
C_ref = C(idxs_to_keep, :);
accuracy_ref = trace(C_ref(:, idxs_to_keep)) / sum(C_ref(:)) * 100;
fprintf('Accord sur les neurones colorés selon %s : %.1f%% (%d neurones, dont %d noirs pour %s)\n', ...
    nom_ref, accuracy_ref, sum(C_ref(:)), sum(C_ref(:, idx_black)), nom_test);

% --- Figure ---
figure;
cm = confusionchart(C, labels_str);
sortClasses(cm, labels_str);
cm.Title = sprintf('%s (Y) contre %s (X), accord sans noires : %.1f%%', ...
    nom_ref, nom_test, accuracy_active);
cm.YLabel = nom_ref;    % lignes = référence
cm.XLabel = nom_test;   % colonnes = méthode comparée
cm.RowSummary    = 'row-normalized';     % par ligne : part de chaque classe de la référence retrouvée (rappel)
cm.ColumnSummary = 'column-normalized';  % par colonne : part des neurones classés j confirmée par la référence (précision)
end
