% 1. Lecture du fichier CSV ligne par ligne sous forme de tableau de chaînes (strings)
clear

% path="E:\Data\Aurelie\data\chroms\119\220923\registration\";
path="E:\Data\Aurelie\data\chroms\582\230331plane0\";
    % name="color_matched_biapy.csv";
name="colorcell_exvivo_romain2026-06-17.csv";

% strColors = readlines("E:\Data\Aurelie\data\chroms\119\220919\registration\colorcell_validated_for_prediction2.csv");
% strColors = readlines("E:\Data\Aurelie\data\chroms\119\220923\registration\color_matched_biapy.csv");

strColors = readlines(strcat(path,name));

%% 


% 2. Nettoyage : suppression des espaces inutiles, des lignes vides et du '0' initial
strColors = strtrim(strColors);
strColors(strColors == "" | strColors == "0") = [];

% 3. Initialisation de la matrice colorcell avec des zéros
colorcell = zeros(length(strColors), 1);

% 4. Assignation des valeurs selon le code couleur demandé
colorcell(strColors == "red")     = 1;
colorcell(strColors == "green")   = 2;
colorcell(strColors == "blue")    = 3;
colorcell(strColors == "yellow")  = 4;
colorcell(strColors == "magenta") = 5;
colorcell(strColors == "cyan")    = 6;
colorcell(strColors == "white")   = 7;
colorcell(strColors == "black")   = 8;

% 5. Affichage des 10 premières valeurs pour vérifier le résultat
% disp('Aperçu des 10 premiers éléments de colorcell :');
% disp(colorcell(1:10));
save (strcat(path,'colorcell_exvivo_romain_20260617'),"colorcell")


