clear; clc;

%% =========================================================
% VÉRIFICATION DU REGROUPEMENT PAR ANIMAL_ID
%
% Ce script lit les sessions EXACTEMENT comme
% correlation_events_loader.m (même database_file, même fonction
% load_recording/get_recording_path), mais SANS charger Fall.mat ni
% lancer preprocessing_6 : il ne sert qu'à diagnostiquer comment
% Animal_ID est peuplé, avant de lancer l'analyse complète.
%
% Il compare Animal_ID (tel que lu dans data_recording) à un code
% animal DEVINÉ à partir de l'identifiant de session (en supposant la
% convention <prefixe>_<numéro_animal>_<AAMMJJ>_planeX, par exemple
% "P2M_411582_230320_plane0" -> animal deviné = "P2M_411582").
% Toute divergence, ou tout repli sur l'identifiant complet en
% l'absence d'Animal_ID, est signalée explicitement.
%% =========================================================

if ispc
    database_file = "E:\Data\Aurelie\data\nocues_data_recording_all.mat";
elseif ismac
    database_file = "/Users/platel/Desktop/exp/brainbow_hippocampus/analysis/final_analysis/bestnocuesold.mat";
else
    error('Système non configuré : définir database_file.');
end

DB = load(database_file);
if isfield(DB,'matfile')
    recording_files = DB.matfile;
elseif isfield(DB,'files')
    recording_files = DB.files;
elseif isfield(DB,'filelist')
    recording_files = DB.filelist;
else
    error('La base %s ne contient pas la variable matfile/files/filelist.',database_file);
end
if isstring(recording_files), recording_files = cellstr(recording_files); end
if ischar(recording_files), recording_files = cellstr(recording_files); end

n = numel(recording_files);
file_num       = (1:n)';
filename_col   = strings(n,1);
identifier_col = strings(n,1);
path_col       = strings(n,1);
animal_field_present = false(n,1);
animal_id_raw  = strings(n,1);   % Animal_ID tel que lu, "" si absent
animal_id_used = strings(n,1);   % ce que correlation_events_loader utiliserait réellement
animal_guess   = strings(n,1);   % deviné à partir de l'identifiant
guess_matches_raw = false(n,1);
status_col     = strings(n,1);

for k = 1:n
    filename_col(k) = string(recording_files{k});
    try
        [Session,recording_file] = load_recording(filename_col(k));
        pathload = get_recording_path(Session,recording_file);
        path_col(k) = string(pathload);

        if isfield(Session,'identifier') && ~isempty(Session.identifier)
            identifier_col(k) = string(Session.identifier);
        else
            identifier_col(k) = "IDENTIFIER_MANQUANT";
            status_col(k) = status_col(k) + "identifier absent; ";
        end

        if isfield(Session,'Animal_ID') && ~isempty(Session.Animal_ID) && string(Session.Animal_ID) ~= ""
            animal_field_present(k) = true;
            animal_id_raw(k) = string(Session.Animal_ID);
            animal_id_used(k) = animal_id_raw(k);
        else
            animal_field_present(k) = false;
            animal_id_raw(k) = "";
            animal_id_used(k) = identifier_col(k); % repli, comme dans le loader
            status_col(k) = status_col(k) + "Animal_ID ABSENT -> repli sur identifier; ";
        end

        animal_guess(k) = guess_animal_from_identifier(identifier_col(k));
        if animal_field_present(k)
            guess_matches_raw(k) = strcmpi(animal_guess(k),animal_id_raw(k)) || ...
                contains(animal_id_raw(k),animal_guess(k),'IgnoreCase',true) || ...
                contains(animal_guess(k),animal_id_raw(k),'IgnoreCase',true);
            if ~guess_matches_raw(k)
                status_col(k) = status_col(k) + sprintf( ...
                    "Animal_ID (''%s'') != animal deviné (''%s''); ", ...
                    animal_id_raw(k),animal_guess(k));
            end
        end

    catch ME
        status_col(k) = "ERREUR : " + string(ME.message);
    end
end

diag_table = table(file_num,filename_col,identifier_col,path_col, ...
    animal_field_present,animal_id_raw,animal_id_used,animal_guess, ...
    guess_matches_raw,status_col, ...
    'VariableNames',{'file_num','filename','identifier','path', ...
    'Animal_ID_present','Animal_ID_raw','Animal_ID_used','Animal_ID_guess', ...
    'guess_matches_raw','status'});

disp(diag_table)

%% =========================================================
% RÉSUMÉ PAR ANIMAL (tel que le pipeline le verrait réellement)
%% =========================================================
fprintf('\n=== RÉSUMÉ ===\n');
fprintf('Nombre total de sessions listées   : %d\n',n);

ok = animal_id_used ~= "" & animal_id_used ~= "IDENTIFIER_MANQUANT";
n_missing_animal = sum(~animal_field_present);
fprintf('Sessions AVEC Animal_ID renseigné  : %d\n',sum(animal_field_present));
fprintf('Sessions SANS Animal_ID (repli)    : %d\n',n_missing_animal);

animals_used = unique(animal_id_used(ok));
fprintf('Nombre d''"animaux" vus par le pipeline (Animal_ID_used) : %d\n', ...
    numel(animals_used));

animals_guess = unique(animal_guess(animal_guess~=""));
fprintf('Nombre d''animaux DEVINÉS depuis l''identifiant           : %d\n', ...
    numel(animals_guess));

if numel(animals_used) >= n*0.8 && n_missing_animal>0
    fprintf(2, ['\nALERTE : le nombre d''"animaux" utilisés par le pipeline ' ...
        '(%d) est proche du nombre de sessions (%d).\n' ...
        'Cela suggère fortement que Animal_ID est absent pour la ' ...
        'plupart des enregistrements, et que summary_animals_true ' ...
        'traite en réalité CHAQUE SESSION comme un animal distinct.\n'], ...
        numel(animals_used),n);
end

fprintf('\n--- Sessions par animal (Animal_ID_used, réellement utilisé) ---\n');
tab_used = groupcounts(animal_id_used(ok));
tbl_used = table(unique(animal_id_used(ok)),tab_used, ...
    'VariableNames',{'Animal_ID_used','n_sessions'});
tbl_used = sortrows(tbl_used,'n_sessions','descend');
disp(tbl_used)

fprintf('\n--- Sessions par animal DEVINÉ (depuis l''identifiant) ---\n');
tab_guess = groupcounts(animal_guess(animal_guess~=""));
tbl_guess = table(unique(animal_guess(animal_guess~="")),tab_guess, ...
    'VariableNames',{'Animal_ID_guess','n_sessions'});
tbl_guess = sortrows(tbl_guess,'n_sessions','descend');
disp(tbl_guess)

n_mismatch = sum(animal_field_present & ~guess_matches_raw);
fprintf(['\nSessions où Animal_ID renseigné NE correspond PAS à ' ...
    'l''animal deviné depuis l''identifiant : %d\n'],n_mismatch);
if n_mismatch>0
    disp(diag_table(animal_field_present & ~guess_matches_raw, ...
        {'file_num','identifier','Animal_ID_raw','Animal_ID_guess','status'}))
end

%% =========================================================
% FIGURE : nombre de sessions par animal, deviné vs utilisé
%% =========================================================
f = figure('Position',[100 100 1000 480],'Color','w');

subplot(1,2,1)
bar(tbl_used.n_sessions,'FaceColor',[0.75 0.35 0.30])
xticks(1:height(tbl_used));
xticklabels(tbl_used.Animal_ID_used);
xtickangle(60);
ylabel('Nombre de sessions');
title(sprintf('Animal\\_ID RÉELLEMENT utilisé par le pipeline (N=%d "animaux")', ...
    numel(animals_used)));
set(gca,'Box','off','TickDir','out','FontSize',8);

subplot(1,2,2)
bar(tbl_guess.n_sessions,'FaceColor',[0.30 0.45 0.75])
xticks(1:height(tbl_guess));
xticklabels(tbl_guess.Animal_ID_guess);
xtickangle(60);
ylabel('Nombre de sessions');
title(sprintf('Animal deviné depuis l''identifiant (N=%d animaux)', ...
    numel(animals_guess)));
set(gca,'Box','off','TickDir','out','FontSize',8);

sgtitle('Vérification du regroupement par animal','FontWeight','bold');

if ispc
    out_dir = "E:\Data\Aurelie\analysis\Final\correlation\";
elseif ismac
    out_dir = "/Users/platel/Desktop/exp/brainbow_hippocampus/analysis/nocues";
end
if ~exist(out_dir,'dir'), mkdir(out_dir); end

exportgraphics(f,fullfile(out_dir,'Check_Animal_Grouping.png'),'Resolution',200);
writetable(diag_table,fullfile(out_dir,'Check_Animal_Grouping_detail.csv'));
writetable(tbl_used,fullfile(out_dir,'Check_Animal_Grouping_used.csv'));
writetable(tbl_guess,fullfile(out_dir,'Check_Animal_Grouping_guess.csv'));

fprintf('\nFichiers écrits dans %s :\n',out_dir);
fprintf('  Check_Animal_Grouping.png\n');
fprintf('  Check_Animal_Grouping_detail.csv\n');
fprintf('  Check_Animal_Grouping_used.csv\n');
fprintf('  Check_Animal_Grouping_guess.csv\n');

%% =========================================================
% FONCTIONS LOCALES
%% =========================================================
function [Session,recording_file] = load_recording(entry)
% Identique à correlation_events_loader.m
entry = string(entry);
if isfolder(entry)
    candidates = [fullfile(entry,'data_recording.mat'),fullfile(entry,'brainbow.mat')];
    recording_file = "";
    for k=1:numel(candidates)
        if isfile(candidates(k)), recording_file=candidates(k); break; end
    end
    if recording_file=="", error('Aucun data_recording.mat/brainbow.mat dans %s.',entry); end
elseif isfile(entry)
    recording_file = entry;
else
    error('Entrée introuvable : %s',entry);
end

S = load(recording_file);
if isfield(S,'data_recording') && isstruct(S.data_recording)
    Session = S.data_recording;
elseif isfield(S,'Session') && isstruct(S.Session)
    Session = S.Session;
else
    Session = S;
end
end

function pathload = get_recording_path(Session,recording_file)
% Identique à correlation_events_loader.m
if isfield(Session,'path') && ~isempty(Session.path)
    pathload = char(string(Session.path));
elseif isfield(Session,'DossierSource') && ~isempty(Session.DossierSource)
    pathload = char(string(Session.DossierSource));
else
    pathload = char(fileparts(recording_file));
end
if ~isfolder(pathload)
    pathload = char(fileparts(recording_file)); % dégradé : ne bloque pas le diagnostic
end
end

function code = guess_animal_from_identifier(identifier)
% Devine le "code animal" à partir d'un identifiant supposé suivre la
% convention <prefixe>_<numéro>_<AAMMJJ>_planeX, ex :
%   "P2M_411582_230320_plane0" -> "P2M_411582"
% Heuristique : on cherche un bloc de 6 chiffres consécutifs (la date),
% et on garde tout ce qui précède, sans le underscore final.
code = "";
identifier = string(identifier);
if identifier=="" || identifier=="IDENTIFIER_MANQUANT"
    return
end

tok = regexp(identifier,'^(.*?)_?\d{6}_?','tokens','once');
if ~isempty(tok) && tok{1} ~= ""
    code = string(tok{1});
else
    % Repli : garde les deux premiers segments séparés par "_"
    parts = split(identifier,'_');
    if numel(parts)>=2
        code = strjoin(parts(1:2),'_');
    else
        code = identifier;
    end
end
end
