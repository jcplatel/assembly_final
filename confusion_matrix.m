function accuracy_active = confusion_matrix (colorcell_method1,colorcell_method2)
%calculation of confusion matrix between 2 colorcell

% ---Définition des labels ---
labels_str = {'Red', 'Green', 'Blue', 'Yellow', 'Magenta', 'Cyan', 'White', 'Black'};
n_classes = length(labels_str);
all_classes = 1:8; 

% --- 3. Calcul de la Matrice de Confusion ---
% C(i,j) = nombre de neurones classés 'i' par méthode 1 et 'j' par méthode 2
C = confusionmat(colorcell_method1, colorcell_method2,'Order', all_classes);
% --- 5. Calcul du Taux d'Accord Global (Accuracy) ---
accuracy = trace(C) / sum(C(:)) * 100;
fprintf('Accord global entre les deux méthodes : %.2f%%\n', accuracy);

idx_black = 8; 

% On crée un masque de tout ce qui n'est PAS Black
% On veut exclure la ligne Black et la colonne Black
C_no_black = C;
C_no_black(idx_black, :) = []; % Supprime la ligne Black
C_no_black(:, idx_black) = []; % Supprime la colonne Black (attention aux indices décalés !)
% MIEUX : Utiliser les indices pour extraire la sous-matrice
idxs_to_keep = setdiff(1:8, idx_black);
C_active = C(idxs_to_keep, idxs_to_keep);

% Calcul de la nouvelle accuracy
accuracy_active = trace(C_active) / sum(C_active(:)) * 100;

fprintf('Accuracy sur les neurones actifs (sans Black) : %.2f%%\n', accuracy_active);
figure;
cm = confusionchart(C, labels_str);
sortClasses(cm, labels_str);

cm.Title = strcat('Comparaison des Méthodes de Classification, accuracy: ',num2str(ceil(accuracy_active)),'%');

cm.XLabel = 'Biapy';%Méthode 1

% cm.YLabel = 'Biapy corrected';%Méthode 2
cm.YLabel = 'in vivo';
cm.RowSummary = 'row-normalized';
cm.ColumnSummary = 'column-normalized';
end