function f = fichier_assemblees(dossier_session_analyse, dossier_session_etape1, K)
%FICHIER_ASSEMBLEES  Fichier des assemblées d'une session pour un K
%   (11 nov.). K = 0 : consensus, dans le dossier de l'analyse. Sinon :
%   version 1 (analyses ordinaires) : <étape 1>\<session>\k<K>\assemblees.mat
%   version 2 (analyse dont le nom est 'v2', 'v2_...' ou '..._v2', voir
%   est_v2) : <étape 1>\<session>\k<K>\assemblees_v2.mat, calculé à la
%   demande à partir de la version 1 (creer_assemblees_v2, sans refaire le
%   k-means) et refait si la version 1 est plus récente.
%   Renvoie le chemin même s'il n'existe pas (version 1 absente).
if K == 0
    f = fullfile(dossier_session_analyse, 'consensus', 'assemblees.mat');
    return
end
f1 = fullfile(dossier_session_etape1, sprintf('k%d', K), 'assemblees.mat');
f = f1;
if ~est_v2(dossier_session_analyse) || ~isfile(f1), return; end
f = fullfile(dossier_session_etape1, sprintf('k%d', K), 'assemblees_v2.mat');
d1 = dir(f1);
d2 = dir(f);
if isempty(d2) || d2.datenum < d1.datenum
    creer_assemblees_v2(f1, f);
end
end
