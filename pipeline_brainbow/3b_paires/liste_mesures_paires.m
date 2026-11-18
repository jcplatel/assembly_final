function [mesures, dt_s] = liste_mesures_paires(opts)
%LISTE_MESURES_PAIRES  Mesures calculées par l'étape 3b : celles de
%   opts.mesures (la première est la principale), puis le STTC aux fenêtres
%   de sensibilité opts.dt_sensibilite_s (ex. 'sttc_100ms', 'sttc_500ms').
%   dt_s : fenêtre de chaque mesure STTC (NaN pour 'dff').
mesures = opts.mesures(:)';
dt_s = nan(1, numel(mesures));
dt_s(strcmp(mesures, 'sttc')) = opts.dt_sttc_s;
if isfield(opts, 'dt_sensibilite_s')
    for dt = opts.dt_sensibilite_s(:)'
        mesures{end+1} = sprintf('sttc_%dms', round(1000 * dt)); %#ok<AGROW>
        dt_s(end+1) = dt; %#ok<AGROW>
    end
end
end
