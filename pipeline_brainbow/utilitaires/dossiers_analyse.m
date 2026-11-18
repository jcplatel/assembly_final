function config = dossiers_analyse(config)
%DOSSIERS_ANALYSE  Où sont lus et écrits les résultats.
%
%   Un dossier par réglage de l'étape 1 (config.PathSave). Les autres
%   analyses faites sur CETTE étape 1 (autre règle de choix du K, autres
%   tests...) sont rangées DEDANS, dans analyses\<nom> :
%
%   <PathSave>\                        ex. ...\resultats
%       <session>\pretraitement.mat    étape 1 (commune à toutes les analyses)
%       <session>\k<K>\assemblees.mat  étape 1
%       <session>\choix_K.mat          étape 2  \
%       <session>\clonal\, paires\     étapes 3, 3b > analyse principale
%       agregation\, paires_agregation\  étape 4  /
%       analyses\<nom>\                autre analyse, même étape 1 :
%           <session>\choix_K.mat, clonal\, paires\
%           agregation\, paires_agregation\
%
%   config.analyse = ''      analyse principale : tout dans PathSave
%   config.analyse = 'nom'   étapes 2, 3, 3b, 4 dans <PathSave>\analyses\nom ;
%                            l'étape 1 reste lue (et écrite) dans PathSave
%   Une analyse ne contient jamais d'étape 1 : c'est toujours celle du
%   dossier qui la contient. Pour la garder en changeant d'ordinateur, on
%   déplace tout le dossier PathSave.
%
%   Après l'appel :
%     config.PathEtape1  dossier de l'étape 1
%     config.PathSave    dossier des étapes 2 à 4 (et du journal)
%     config.PathSource  = PathEtape1 (lu par lire_table_sessions)
%   Un deuxième appel ne change rien.
%   Anciens scripts : config.PathSource donné sans config.analyse -> étape 1
%   lue dans PathSource, étapes 2 à 4 dans PathSave (comme avant).

if isfield(config, 'PathEtape1') && ~isempty(config.PathEtape1)
    return      % déjà fait
end
nom = '';
if isfield(config, 'analyse') && ~isempty(config.analyse)
    nom = strtrim(char(config.analyse));
end
% Version des assemblées (11 nov.) : la version 2 (cellules partagées entre
% assemblées, RACEORTHO_V2) est une propriété de l'analyse ; son nom est
% 'v2', commence par 'v2_' ou finit par '_v2' (voir est_v2). Avec
% config.version_assemblees = 'v2', '_v2' est ajouté au nom s'il n'y est
% pas ('v2' sans nom) ; inversement, choisir une analyse de ce nom donne la
% version 2.
nom_est_v2 = ~isempty(regexpi(nom, '(^v2$|^v2_|_v2$)', 'once'));
if isfield(config, 'version_assemblees') && strcmpi(char(config.version_assemblees), 'v2') && ~nom_est_v2
    if isempty(nom), nom = 'v2'; else, nom = [nom '_v2']; end
    nom_est_v2 = true;
end
if nom_est_v2, config.version_assemblees = 'v2'; else, config.version_assemblees = 'v1'; end
if ~isempty(nom)
    if ~isempty(regexp(nom, '[\\/:*?"<>|]', 'once')) || any(strcmp(nom, {'.', '..'}))
        error('Nom d''analyse invalide : "%s" (sans \\ / : * ? " < > |).', nom);
    end
    config.analyse = nom;
    config.PathEtape1 = config.PathSave;
    config.PathSave = fullfile(config.PathSave, 'analyses', nom);
elseif isfield(config, 'PathSource') && ~isempty(config.PathSource)
    config.analyse = '';
    config.PathEtape1 = config.PathSource;
else
    config.analyse = '';
    config.PathEtape1 = config.PathSave;
end
config.PathSource = config.PathEtape1;
end
