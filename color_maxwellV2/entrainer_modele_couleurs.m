function modele = entrainer_modele_couleurs(dossiers_GT, fichier_modele)
%% =========================================================================
%  ENTRAINER_MODELE_COULEURS (17 nov., demandée par l'utilisateur)
%  Apprend le LDA de prediction_color (code de l'utilisateur) sur les
%  cellules colorées (classes 1 à 6) de colorcell_GT.mat de plusieurs
%  sessions, avec le prétraitement de variables_couleurs.
%  Étapes 1 à 3 comme interface_couleurs : brainbow_image_read,
%  brainbow_image_processing, image_alignment_multimodal, ROI_extraction2.
%
%  Lancer :  entrainer_modele_couleurs({'E:\...\220919_plane0', ...}, ...
%                'E:\...\modele_couleurs_LDA.mat')
%  Validation (17 nov.) : appris sur les autres animaux, ce LDA donne la
%  bonne couleur à 82-87 % des cellules colorées du GT (444119, 444175) ;
%  57 % sur 411582 230331 (bruitée, exclue de l'apprentissage).
%  =========================================================================

X = []; y = []; par_session = struct('dossier', {}, 'n_cellules', {}, 'n_appris', {});
for i = 1:numel(dossiers_GT)
    d = [regexprep(char(dossiers_GT{i}), '[\\/]+$', '') filesep];
    [R, G, B, Ca] = valeurs_roi(d);
    [Xi, valide] = variables_couleurs(R, G, B, Ca);
    g = load([d 'colorcell_GT.mat'], 'colorcell'); g = double(g.colorcell(:));
    if numel(g) ~= numel(R)
        error('%s : colorcell_GT a %d valeurs, ROI_extraction2 %d cellules.', d, numel(g), numel(R));
    end
    k = valide & ismember(g, 1:6);
    X = [X; Xi(k, :)]; y = [y; g(k)]; %#ok<AGROW>
    par_session(end+1) = struct('dossier', d, 'n_cellules', numel(g), 'n_appris', sum(k)); %#ok<AGROW>
    fprintf('%s : %d cellules, %d apprises\n', d, numel(g), sum(k));
end

modele = struct();
modele.Mdl = fitcdiscr(X, y, 'DiscrimType', 'linear');
modele.type = 'LDA (fitcdiscr linear), variables de prediction_color';
modele.sessions = par_session;
modele.effectifs_classes = histcounts(y, 0.5:1:6.5);
modele.seuil_certitude = 0.5;      % en dessous : noire (comme colorcell_maxwell4)
modele.date = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<DATST,TNOW1>
modele.fonctions = 'variables_couleurs, entrainer_modele_couleurs (color_maxwellV2)';
if nargin >= 2 && ~isempty(fichier_modele)
    save(fichier_modele, '-struct', 'modele');
    fprintf('Modèle enregistré : %s\n', fichier_modele);
end
end


function [R, G, B, Ca] = valeurs_roi(d)
[red, red890, green, blue, calcium] = brainbow_image_read(d);
In = brainbow_image_processing(d, red, red890, green, blue, calcium);
[ar, ag, ab] = image_alignment_multimodal(In);
[R, G, B, Ca] = ROI_extraction2(d, ar, ag, ab, In.calcium);
end
