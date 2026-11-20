function nom = dossier_sequences(opts)
%DOSSIER_SEQUENCES  Nom du dossier des résultats de l'étape 5 (20 nov. 2026).
%   'sequences' pour le calcium (opts.signal = 'calcium', comme avant) ;
%   'sequences_cascade' pour la probabilité de spike de CASCADE. Le dossier
%   de l'agrégation suit : 'sequences_agregation' ou
%   'sequences_agregation_cascade'.
nom = 'sequences';
if isfield(opts, 'signal') && strcmp(opts.signal, 'cascade')
    nom = 'sequences_cascade';
end
end
