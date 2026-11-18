function bilan = creer_data_recording_liste(fichiers_nwb)
%% =========================================================================
%  CREER_DATA_RECORDING_LISTE (18 nov. 2026, sessions cues)
%  Comme create_data_recording_from_nwb.m (même lecture : openingnwb, même
%  contenu de data_recording), mais pour une LISTE de fichiers NWB, et :
%   - si un dossier contient plusieurs NWB de la liste, seul celui dont le
%     nombre de cellules (fluorescence_all_cells) = iscell du Fall.mat est gardé
%   - data_recording.mat existant avec le bon nombre de cellules : gardé tel quel
%   - data_recording.mat existant avec un AUTRE nombre de cellules que le
%     Fall.mat : renommé data_recording_ancien_<n>cellules.mat (rien n'est
%     effacé), puis refait
%  bilan : une ligne par dossier (état, nombre de cellules Fall / data_recording)
%  =========================================================================
fichiers_nwb = cellstr(fichiers_nwb(:));
dossiers = unique(cellfun(@fileparts, fichiers_nwb, 'UniformOutput', false), 'stable');
n = numel(dossiers);
bilan = table(dossiers, repmat({''}, n, 1), nan(n, 1), nan(n, 1), repmat({''}, n, 1), ...
    'VariableNames', {'dossier', 'etat', 'n_iscell', 'n_data_recording', 'nwb'});
for i = 1:n
    d = dossiers{i};
    try
        Fa = load(fullfile(d, 'Fall.mat'), 'F', 'iscell');
        garde = Fa.iscell(:, 1) > 0;
        bilan.n_iscell(i) = sum(garde);
        fdr = fullfile(d, 'data_recording.mat');
        if isfile(fdr)
            m = matfile(fdr); x = m.data_recording; nx = size(x.F, 1);
            if nx == sum(garde)
                bilan.etat{i} = 'existant, bon nombre de cellules';
                bilan.n_data_recording(i) = nx;
                continue
            end
            ancien = fullfile(d, sprintf('data_recording_ancien_%dcellules.mat', nx));
            movefile(fdr, ancien);
        end
        % NWB de ce dossier : celui qui a le même nombre de cellules que le Fall.mat
        cand = fichiers_nwb(strcmp(cellfun(@fileparts, fichiers_nwb, 'UniformOutput', false), d));
        cand = cand(cellfun(@isfile, cand));
        if isempty(cand)            % NWB de la liste absent : ceux du dossier
            l = dir(fullfile(d, '*.nwb')); l = l(~startsWith({l.name}, '._'));
            cand = fullfile(d, {l.name})';
        end
        choisi = '';
        if numel(cand) == 1
            choisi = cand{1};
        else
            for k = 1:numel(cand)
                nwb = nwbRead(cand{k}, 'ignorecache');
                o = nwb.processing.get('ophys');
                r = o.nwbdatainterface.get('fluorescence_all_cells').roiresponseseries.get('f_raw_suite2p');
                if r.data.dims(1) == sum(garde), choisi = cand{k}; break; end
            end
            if isempty(choisi), error('aucun NWB n''a %d cellules', sum(garde)); end
        end
        bilan.nwb{i} = choisi;
        [identifier, sampling_rate, speed, position] = lire_nwb(choisi, d); %#ok<ASGLU>
        data_recording = struct();
        data_recording.identifier    = identifier;
        data_recording.F             = double(Fa.F(garde, :));
        data_recording.sampling_rate = sampling_rate;
        data_recording.speed         = speed;
        data_recording.position      = position;
        data_recording.path          = d;
        save(fdr, 'data_recording', '-v7.3');
        bilan.n_data_recording(i) = size(data_recording.F, 1);
        if exist('ancien', 'var') && isfile(ancien)
            bilan.etat{i} = 'refait (ancien renommé)'; clear ancien
        else
            bilan.etat{i} = 'créé';
        end
    catch err
        bilan.etat{i} = ['ERREUR : ' err.message];
    end
    fprintf('%d/%d %s : %s\n', i, n, d, bilan.etat{i});
end
end


function [identifier, sampling_rate, speed, position] = lire_nwb(filename, ~)
% Mêmes lectures que le script openingnwb, mais nwbRead(..., 'ignorecache') :
% les classes matnwb ne sont pas régénérées (les fichiers 2025 font échouer
% la régénération : « Unrecognized property 'data_offset' »). Classes du
% schéma 2.6.0 à mettre en tête du chemin avant l'appel (generateCore).
read_nwb = nwbRead(filename, 'ignorecache');
identifier = read_nwb.identifier;
position = [];
try
    position = read_nwb.processing.get('behavior').nwbdatainterface.get('Position') ...
        .spatialseries.get('mouse_position').data.load;
catch
end
speed = read_nwb.processing.get('behavior').nwbdatainterface.get('BehavioralTimeSeries') ...
    .timeseries.get('running_speed_by_frame').data.load;
t = read_nwb.processing.get('ophys').nwbdatainterface.get('fluorescence_all_cells') ...
    .roiresponseseries.get('f_raw_suite2p').timestamps.load;
sampling_rate = length(t) / t(end);
end
