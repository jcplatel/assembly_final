function redessiner_figures_etape4(dossier_variante, config)
%REDESSINER_FIGURES_ETAPE4  Redessine les figures de l'étape 4 d'une variante
%   (fenêtres MATLAB : zoom, curseur), sans refaire aucun calcul : lit
%   resultats_biais_clonal.mat de dossier_variante. Utilisé par le bouton
%   « Results figures » de l'interface.
%
%   redessiner_figures_etape4('E:\...\agregation\champ_globale', config_brainbow())
%
%   Mêmes seuils que l'étape 4 (config.agreg, config.clonal.seuilstat).
%   Les PNG du dossier sont réécrits.

if nargin < 2 || isempty(config), config = config_brainbow(); end
dossier_racine = regexprep(char(dossier_variante), '[\\/]+$', '');
[~, nom] = fileparts(dossier_racine);
% variables lues par le script figures_biais_clonal (comme dans l'étape 4)
fdr_seuil_figure = config.agreg.fdr_seuil_figure; %#ok<NASGU>
seuil_non_biaise = config.agreg.seuil_non_biaise; %#ok<NASGU>
alpha = config.clonal.seuilstat;
v = config.clonal.variantes;
if strcmp(regexprep(nom, '_PROVISOIRE$', ''), [v{1, 1} '_' v{1, 2}])
    nom_variante = [nom ' (TEST PRINCIPAL)'];
else
    nom_variante = [nom ' (contrôle)'];
end
figures_biais_clonal
figure_assemblees_par_session(dossier_racine, config.agreg.fdr_seuil_figure, alpha, nom_variante);
figure_biaisees_par_animal(dossier_racine, config.agreg.fdr_seuil_figure, alpha, nom_variante);
figures_article(dossier_racine, config.agreg.fdr_seuil_figure, alpha, nom_variante);
% contrôles du champ (communs à toutes les variantes)
dossier_ctrl = fullfile(fileparts(dossier_racine), 'controles_champ');
if isfile(fullfile(dossier_ctrl, 'controles_champ.mat'))
    figures_controles_champ(dossier_ctrl);
end
if isfile(fullfile(dossier_ctrl, 'compacite_assemblees.mat'))
    figure_compacite(dossier_ctrl);
end
end
