function ecrire_analysisALL(sessions, config)
%ECRIRE_ANALYSISALL  analysisALL.xlsx dans le dossier de l'étape 1 : une
%   ligne par session et par K déjà calculés (ancien exportdata).
%   Appelé après CHAQUE session de l'étape 1 : le fichier grossit au fur et à
%   mesure de l'analyse. Il est réécrit entièrement à chaque fois (pas de
%   doublons). Si le fichier est ouvert dans Excel, il n'est pas mis à jour
%   (simple avertissement) : il le sera après la session suivante.
dossier_etape1 = sessions.dossier;
if ismember('dossier_source', sessions.Properties.VariableNames), dossier_etape1 = sessions.dossier_source; end
racine = config.PathSave;
if isfield(config, 'PathEtape1') && ~isempty(config.PathEtape1), racine = config.PathEtape1; end
lignes = {};
for s = 1:height(sessions)
    for ck = config.liste_K
        fichier_K = fullfile(dossier_etape1{s}, sprintf('k%d', ck), 'assemblees.mat');
        if ~isfile(fichier_K), continue; end
        data = load(fichier_K, 'resume_K');
        if ~isfield(data, 'resume_K'), continue; end
        lignes{end+1, 1} = cell2table(struct2cell(data.resume_K)', ...
            'VariableNames', fieldnames(data.resume_K)'); %#ok<AGROW>
    end
end
if isempty(lignes), return; end
path_all = fullfile(racine, 'analysisALL.xlsx');
try
    T = vertcat(lignes{:});
    if isfile(path_all), delete(path_all); end
    writetable(T, path_all);
catch err
    warning('analysisALL.xlsx non mis à jour (ouvert dans Excel ?) : %s', err.message);
end
end
