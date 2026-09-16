
%% ASSEMBLY_DATABASE_MANAGER - GUI for managing the imaging sessions database
% DESCRIPTION:
%   This script launches a two-tab graphical user interface (uifigure) that lets the user scan a directory of Suite2P/
% Assembly results and build/update a centralized MATLAB database (.mat file), then manually edit selected fields of 
% that database through an interactive table (uitable).
%
%   Tab 1 - "Scanner / Mettre à jour" (Scan / Update):
%       Recursively searches a source folder for 'results.mat' files located inside 'k6' subfolders. For each session found, it
%       creates or updates an entry (a 'Session' struct) in the database, loading only the variables that are still missing
%       ('identifier', 'F', 'iscell', 'sampling_rate', 'speed','position', 'path') to minimize loading time. It also searches
%       for the best available 'colorcell*.mat' file according to a predefined preference order.
%
%   Tab 2 - "Éditeur de Base de Données" (Database Editor):
%       Displays all sessions in the database as an editable table (Animal_ID, Colorcell_Version, Suite2P_Version, analysis flags
%       Assembly/Correlation/Sequence/PCs, Comment), allowing the user to flag sessions for permanent deletion and to save changes
%       back to the .mat file (including physical cleanup of the file after deletions).
%
% INPUT:
%   Upstream, the script relies on 'results.mat' files already generated
%   by the Suite2P/Assembly pipeline (inside 'k<N>' subfolders) and, when
%   available, on 'colorcell*.mat' files present in the session folder.
%
% OUTPUT:
% The script creates/updates a database file in .mat format ('Global_Database.mat', path defined in the global variable
%   'path_save'), containing one struct-type variable per session (named 'S_<identifier>'), with the fields:
%     Condition, DossierSource, DateExtraction, Suite2P_Version, Commentaire, Animal_ID, Analysis_Assembly, Analysis_Correlation,
%     Analysis_Sequence, Analysis_PCs, K, identifier, F, iscell, sampling_rate, speed, position, path, Colorcell_Version.
%   The file is updated incrementally ('-append') during scans, and fully rewritten ('-v7.3') after deleting sessions, to physically
%   free up disk space.
%%
clear; clc; close all;

%% 1. PARAMÈTRES GLOBAUX
global path_save
path_save = 'E:\Data\Aurelie\data\Global_Database.mat';

%% 2. CRÉATION DE L'INTERFACE GRAPHIQUE PRINCIPALE
f = uifigure('Name', 'Gestionnaire de Base de Données (Assembly)', 'Position', [50 50 1400 700]);

% Création du groupe d'onglets
tabgroup = uitabgroup(f, 'Position', [0 0 1400 700]);

% --- ONGLET 1 : MISE À JOUR / SCAN ---
tab1 = uitab(tabgroup, 'Title', '1. Scanner / Mettre à jour');
construire_onglet_scan(tab1);

% --- ONGLET 2 : ÉDITEUR INTERACTIF ---
tab2 = uitab(tabgroup, 'Title', '2. Éditeur de Base de Données');
construire_onglet_editeur(tab2);


%% ========================================================================
% FONCTIONS DE CONSTRUCTION DES ONGLETS
% =========================================================================

function construire_onglet_scan(parent_tab)
    g = uigridlayout(parent_tab, [7 3], 'RowHeight', {30, 30, 30, 30, 30, 50, '1x'}, 'ColumnWidth', {150, '1x', 100});
    
    % Dossier Source
    uilabel(g, 'Text', 'Dossier Source :');
    champ_dossier = uieditfield(g, 'text', 'Value', 'E:\Data\Aurelie\analysis\June2026\assembly\nocues\');
    champ_dossier.Layout.Column = 2;
    
    btn_browse = uibutton(g, 'Text', '📁 Parcourir...');
    btn_browse.Layout.Column = 3;
    btn_browse.ButtonPushedFcn = @(src, event) choisir_dossier(champ_dossier);
    
    % Condition
    lbl_cond = uilabel(g, 'Text', 'Condition :');
    lbl_cond.Layout.Column = 1;
    champ_cond = uieditfield(g, 'text', 'Value', 'nocues');
    champ_cond.Layout.Column = [2 3]; 
    
    % Version Suite2P
    lbl_s2p = uilabel(g, 'Text', 'Version Suite2P :');
    lbl_s2p.Layout.Column = 1;
    champ_s2p = uidropdown(g, 'Items', {'Aurélie', 'JC'});
    champ_s2p.Layout.Column = [2 3];
    
    % Commentaire global
    lbl_com = uilabel(g, 'Text', 'Commentaire global :');
    lbl_com.Layout.Column = 1;
    champ_comment = uieditfield(g, 'text', 'Value', 'Scan initial');
    champ_comment.Layout.Column = [2 3];
    
    % Options
    lbl_opt = uilabel(g, 'Text', 'Options :');
    lbl_opt.Layout.Column = 1;
    cb_forcer = uicheckbox(g, 'Text', 'Forcer la mise à jour des métadonnées des sessions existantes', 'Value', false);
    cb_forcer.Layout.Column = [2 3];
    
    % Bouton Lancer
    btn_scan = uibutton(g, 'Text', '🚀 Lancer le Scan et la Mise à Jour', 'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', [0.8 1 0.8]);
    btn_scan.Layout.Column = [1 3];
    
    % Zone de log
    txt_log = uitextarea(g, 'Editable', 'off', 'FontName', 'Consolas');
    txt_log.Layout.Column = [1 3];
    txt_log.Value = {'Prêt à scanner...'};
    
    btn_scan.ButtonPushedFcn = @(src, event) lancer_scan(champ_dossier.Value, champ_cond.Value, ...
        champ_s2p.Value, champ_comment.Value, cb_forcer.Value, txt_log);
end

function choisir_dossier(champ_texte)
    dossier_choisi = uigetdir(champ_texte.Value, 'Sélectionner le dossier source à scanner');
    if dossier_choisi ~= 0
        champ_texte.Value = dossier_choisi;
    end
end

function construire_onglet_editeur(parent_tab)
    g = uigridlayout(parent_tab, [3 1], 'RowHeight', {40, '1x', 50});
    
    btn_load = uibutton(g, 'Text', '🔄 Recharger la Base de Données', 'FontSize', 12);
    btn_load.Layout.Row = 1;
    
    uit = uitable(g);
    uit.Layout.Row = 2;
    
    btn_save = uibutton(g, 'Text', '💾 SAUVEGARDER (Modifications et Suppressions)', ...
               'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', [0.9 0.9 0.9]);
    btn_save.Layout.Row = 3;
    btn_save.Enable = 'off'; 
    
    btn_load.ButtonPushedFcn = @(src, event) charger_table_editeur(uit, btn_save);
    btn_save.ButtonPushedFcn = @(src, event) sauvegarder_et_supprimer(uit);
    
    global path_save;
    if isfile(path_save)
        charger_table_editeur(uit, btn_save);
    end
end


%% ========================================================================
% LOGIQUE DU SCAN (ONGLET 1)
% =========================================================================

function lancer_scan(dossier_source, condition, version_s2p, commentaire, forcer_mise_a_jour, txt_log)
    global path_save;
    
    if ~isfolder(dossier_source)
        uialert(gcbf, 'Le dossier source spécifié n''existe pas.', 'Erreur de dossier');
        return;
    end
    
    log_msg(txt_log, sprintf('Recherche globale dans %s...', dossier_source));
    
    % Recherche récursive des deux types de fichiers
    liste_rec = dir(fullfile(dossier_source, '**', 'data_recording.mat'));
    liste_res = dir(fullfile(dossier_source, '**', 'results.mat'));
    
    log_msg(txt_log, sprintf('Trouvé : %d "data_recording.mat" (bruts) et %d "results.mat" (analyses).', ...
        length(liste_rec), length(liste_res)));
    
    if isempty(liste_rec) && isempty(liste_res)
        log_msg(txt_log, '⚠️ Aucun fichier pertinent trouvé.');
        return;
    end
    
    if ~isfile(path_save)
        log_msg(txt_log, sprintf('Création de la base : %s', path_save));
        DB_Info = sprintf('Database créée le %s', datestr(now));
        save(path_save, 'DB_Info', '-v7.3');
        variables_existantes = {}; 
    else
        variables_existantes = who('-file', path_save);
    end
    
    nb_ajouts = 0; nb_ignores = 0; nb_mis_a_jour = 0;
    
    % Construire la liste de tous les éléments à traiter
    % Structure : chemin_fichier, type ('RAW' ou 'ANALYSIS')
    elements = struct('chemin', {}, 'type', {});
    
    % Priorité 1 : On intègre les manips brutes
    for r = 1:length(liste_rec)
        elements(end+1).chemin = fullfile(liste_rec(r).folder, liste_rec(r).name); %#ok<AGROW>
        elements(end).type = 'RAW';
    end
    % Priorité 2 : On intègre/met à jour avec les analyses
    for s = 1:length(liste_res)
        elements(end+1).chemin = fullfile(liste_res(s).folder, liste_res(s).name); %#ok<AGROW>
        elements(end).type = 'ANALYSIS';
    end
    
    total_elements = length(elements);
    h_wait = waitbar(0, 'Mise à jour de la base de données...');
    
    for i = 1:total_elements
        elem = elements(i);
        dossier_parent = fileparts(elem.chemin);
        
        try
            id_session = '';
            a_ete_modifiee = false;
            est_nouvelle = false;
            
            % =============================================================
            % CAS A : FICHIER BRUT (data_recording.mat) -> Non analysé ou brut
            % =============================================================
            if strcmp(elem.type, 'RAW')
                rec_data = load(elem.chemin);
                if isfield(rec_data, 'data_recording')
                    struct_rec = rec_data.data_recording;
                else
                    struct_rec = rec_data;
                end
                
                if isfield(struct_rec, 'identifier') && ~isempty(struct_rec.identifier)
                    id_session = char(struct_rec.identifier);
                else
                    [~, id_session] = fileparts(dossier_parent);
                end
                
                nom_variable_session = matlab.lang.makeValidName(['S_' id_session]);
                
                % Chargement ou initialisation
                if ismember(nom_variable_session, variables_existantes)
                    Session = load(path_save, nom_variable_session).(nom_variable_session);
                else
                    Session = struct();
                    est_nouvelle = true;
                    a_ete_modifiee = true;
                end
                
                % Remplissage métadonnées si nouvelle ou forcé
                if est_nouvelle || forcer_mise_a_jour
                    Session.Condition = condition;
                    Session.DateExtraction = datestr(now, 'dd/mm/yyyy HH:MM');
                    Session.Suite2P_Version = version_s2p;
                    if ~isfield(Session, 'Commentaire') || isempty(Session.Commentaire)
                        Session.Commentaire = commentaire;
                    end
                    a_ete_modifiee = true;
                end
                
                % Variables indispensables depuis data_recording
                Session.identifier = id_session;
                Session.DossierSource = dossier_parent;
                
                champs_rec = {'F', 'sampling_rate', 'speed', 'position', 'path'};
                for c = 1:length(champs_rec)
                    f_name = champs_rec{c};
                    if isfield(struct_rec, f_name) && (~isfield(Session, f_name) || isempty(Session.(f_name)) || forcer_mise_a_jour)
                        Session.(f_name) = struct_rec.(f_name);
                        a_ete_modifiee = true;
                    end
                end
                
                % Pour une session non analysée : on garantit des champs par défaut
                if ~isfield(Session, 'iscell'), Session.iscell = []; a_ete_modifiee = true; end
                if ~isfield(Session, 'K'), Session.K = NaN; a_ete_modifiee = true; end
                
                % Recherche Brainbow / Colorcell dans le dossier brut
                if ~isfield(Session, 'Colorcell_Version') || isempty(Session.Colorcell_Version) || forcer_mise_a_jour
                    Session.Colorcell_Version = trouver_meilleur_colorcell(dossier_parent);
                    a_ete_modifiee = true;
                end
                
            % =============================================================
            % CAS B : FICHIER ANALYSE (results.mat)
            % =============================================================
            elseif strcmp(elem.type, 'ANALYSIS')
                res_meta = load(elem.chemin, 'identifier', 'path', 'iscell');
                
                if isfield(res_meta, 'identifier') && ~isempty(res_meta.identifier)
                    id_session = char(res_meta.identifier);
                else
                    [~, nom_dossier_sup] = fileparts(fileparts(elem.chemin));
                    id_session = nom_dossier_sup;
                end
                
                nom_variable_session = matlab.lang.makeValidName(['S_' id_session]);
                
                if ismember(nom_variable_session, variables_existantes)
                    Session = load(path_save, nom_variable_session).(nom_variable_session);
                else
                    Session = struct();
                    Session.identifier = id_session;
                    est_nouvelle = true;
                    a_ete_modifiee = true;
                end
                
                if est_nouvelle || forcer_mise_a_jour
                    Session.Condition = condition;
                    Session.DateExtraction = datestr(now, 'dd/mm/yyyy HH:MM');
                    Session.Suite2P_Version = version_s2p;
                    if ~isfield(Session, 'Commentaire') || isempty(Session.Commentaire)
                        Session.Commentaire = commentaire;
                    end
                    a_ete_modifiee = true;
                end
                
                Session.DossierSource = dossier_parent; % pointe vers k6
                
                % Extraction K
                k_tokens = regexp(dossier_parent, '[\\/]k(\d+)', 'tokens');
                if ~isempty(k_tokens)
                    Session.K = str2double(k_tokens{end}{1});
                    a_ete_modifiee = true;
                end
                
                % Extraction iscell
                if isfield(res_meta, 'iscell')
                    Session.iscell = res_meta.iscell;
                    a_ete_modifiee = true;
                end
                
                % Récupération du dossier brut via path si F ou d'autres signaux manquent
                dossier_brut = '';
                if isfield(res_meta, 'path') && ~isempty(res_meta.path)
                    dossier_brut = char(res_meta.path);
                    Session.path = dossier_brut;
                elseif isfield(Session, 'path') && ~isempty(Session.path)
                    dossier_brut = char(Session.path);
                end
                
                % Compléter avec data_recording.mat si besoin
                if ~isempty(dossier_brut) && isfile(fullfile(dossier_brut, 'data_recording.mat'))
                    if ~isfield(Session, 'F') || isempty(Session.F) || forcer_mise_a_jour
                        tmp_brut = load(fullfile(dossier_brut, 'data_recording.mat'));
                        if isfield(tmp_brut, 'data_recording'), s_b = tmp_brut.data_recording; else, s_b = tmp_brut; end
                        
                        champs_rec = {'F', 'sampling_rate', 'speed', 'position'};
                        for c = 1:length(champs_rec)
                            fn = champs_rec{c};
                            if isfield(s_b, fn), Session.(fn) = s_b.(fn); a_ete_modifiee = true; end
                        end
                    end
                end
                
                % Recherche Colorcell (Dossier analyse puis Dossier brut)
                if ~isfield(Session, 'Colorcell_Version') || isempty(Session.Colorcell_Version) || forcer_mise_a_jour
                    cc_found = trouver_meilleur_colorcell(fileparts(dossier_parent));
                    if strcmp(cc_found, 'Non trouvé') && ~isempty(dossier_brut)
                        cc_found = trouver_meilleur_colorcell(dossier_brut);
                    end
                    Session.Colorcell_Version = cc_found;
                    a_ete_modifiee = true;
                end
            end
            
            % =============================================================
            % CHAMPS STANDARDS & ÉCRITURE BDD
            % =============================================================
            if ~isfield(Session, 'Animal_ID'), Session.Animal_ID = ''; a_ete_modifiee = true; end
            if ~isfield(Session, 'Analysis_Assembly'), Session.Analysis_Assembly = false; a_ete_modifiee = true; end
            if ~isfield(Session, 'Analysis_Correlation'), Session.Analysis_Correlation = false; a_ete_modifiee = true; end
            if ~isfield(Session, 'Analysis_Sequence'), Session.Analysis_Sequence = false; a_ete_modifiee = true; end
            if ~isfield(Session, 'Analysis_PCs'), Session.Analysis_PCs = false; a_ete_modifiee = true; end
            
            if a_ete_modifiee
                eval([nom_variable_session ' = Session;']);
                save(path_save, nom_variable_session, '-append');
                if est_nouvelle
                    variables_existantes{end+1} = nom_variable_session; 
                    nb_ajouts = nb_ajouts + 1;
                else
                    nb_mis_a_jour = nb_mis_a_jour + 1;
                end
            else
                nb_ignores = nb_ignores + 1;
            end
            
        catch exception
            log_msg(txt_log, sprintf('⚠️ Erreur sur %s : %s', elem.chemin, exception.message));
        end
        waitbar(i / total_elements, h_wait);
    end
    
    close(h_wait);
    log_msg(txt_log, sprintf('✅ Scan terminé : %d ajouts | %d MAJ | %d ignorés.', nb_ajouts, nb_mis_a_jour, nb_ignores));
    uialert(gcbf, 'Scan terminé avec succès ! Cliquez sur "Recharger" dans l''onglet Éditeur.', 'Succès');
end

function nom_fichier_trouve = trouver_meilleur_colorcell(dossier_cible)
    ordre_preference = {'colorcell_GT*.mat','colorcell_Maxwellnew*.mat','colorcell_Maxwell*.mat','colorcellnew*.mat','colorcell.mat'};
    nom_fichier_trouve = 'Non trouvé'; 
    if ~isfolder(dossier_cible), return; end
    
    for p = 1:length(ordre_preference)
        fichiers_trouves = dir(fullfile(dossier_cible, ordre_preference{p}));
        if ~isempty(fichiers_trouves)
            nom_fichier_trouve = fichiers_trouves(1).name; return;
        end
    end
    fichiers_presents = dir(fullfile(dossier_cible, 'colorcell*.mat'));
    if ~isempty(fichiers_presents), nom_fichier_trouve = fichiers_presents(1).name; end
end

function log_msg(txt_area, msg)
    txt_area.Value = [txt_area.Value; {msg}];
    scroll(txt_area, 'bottom');
    drawnow;
end

%% ========================================================================
% LOGIQUE DE L'ÉDITEUR (ONGLET 2)
% =========================================================================

function charger_table_editeur(uit, btn_save)
    global path_save;
    if ~isfile(path_save), uialert(gcbf, 'BDD introuvable. Faites un scan d''abord.', 'Erreur'); return; end
    
    h_wait = waitbar(0, 'Chargement de la base de données...');
    liste_variables = who('-file', path_save);
    noms_sessions = liste_variables(startsWith(liste_variables, 'S_'));
    nb_sessions = length(noms_sessions);
    
    % Pré-allocation
    Supprimer  = false(nb_sessions, 1);
    Variable_Name = string(noms_sessions);
    Animal_ID  = strings(nb_sessions, 1);
    Session_ID = strings(nb_sessions, 1);
    Condition  = strings(nb_sessions, 1);
    Nb_Tot_Cells = NaN(nb_sessions, 1);
    NCell      = NaN(nb_sessions, 1);
    Colorcell  = strings(nb_sessions, 1);
    S2P_Vers   = strings(nb_sessions, 1);
    Ana_Assembly = false(nb_sessions, 1);
    Ana_Correl   = false(nb_sessions, 1);
    Ana_Sequence = false(nb_sessions, 1);
    Ana_PCs      = false(nb_sessions, 1);
    Comment    = strings(nb_sessions, 1);
    
    for i = 1:nb_sessions
        nom_var = noms_sessions{i};
        data = load(path_save, nom_var);
        Session = data.(nom_var);
        
        Session_ID(i) = strrep(nom_var, 'S_', '');
        
        if isfield(Session, 'Animal_ID'), Animal_ID(i) = string(Session.Animal_ID); else Animal_ID(i) = ""; end
        if isfield(Session, 'Condition'), Condition(i) = string(Session.Condition); else Condition(i) = ""; end
        
        if isfield(Session, 'F') && ~isempty(Session.F)
            Nb_Tot_Cells(i) = size(Session.F, 1); 
        end
        
        % Extraction NCell (Somme de iscell(:,1) > 0)
        if isfield(Session, 'iscell') && ~isempty(Session.iscell)
            NCell(i) = sum(Session.iscell(:, 1) > 0); 
        end
        
        if isfield(Session, 'Colorcell_Version'), Colorcell(i) = string(Session.Colorcell_Version); else Colorcell(i) = ""; end
        if isfield(Session, 'Suite2P_Version'), S2P_Vers(i) = string(Session.Suite2P_Version); else S2P_Vers(i) = ""; end
        
        if isfield(Session, 'Analysis_Assembly') && Session.Analysis_Assembly == true, Ana_Assembly(i) = true; end
        if isfield(Session, 'Analysis_Correlation') && Session.Analysis_Correlation == true, Ana_Correl(i) = true; end
        if isfield(Session, 'Analysis_Sequence') && Session.Analysis_Sequence == true, Ana_Sequence(i) = true; end
        if isfield(Session, 'Analysis_PCs') && Session.Analysis_PCs == true, Ana_PCs(i) = true; end
        
        if isfield(Session, 'Commentaire'), Comment(i) = string(Session.Commentaire); else Comment(i) = ""; end
        waitbar(i / nb_sessions, h_wait);
    end
    close(h_wait);
    
    % Création de la table
    T = table(Supprimer, Variable_Name, Animal_ID, Session_ID, Condition, ...
              Nb_Tot_Cells, NCell, Colorcell, S2P_Vers, ...
              Ana_Assembly, Ana_Correl, Ana_Sequence, Ana_PCs, Comment);
    uit.Data = T;
    
    col_editable = [true, false, true, false, false, false, false, true, true, true, true, true, true, true]; 
    uit.ColumnEditable = col_editable;
    
    uit.ColumnWidth = {40, 0, 100, 160, 80, 50, 50, 150, 80, 60, 60, 60, 60, '1x'};
    
    uit.ColumnName = {'Suppr', 'VarName', 'Animal ID', 'Session ID', 'Condition', ...
                      'Total', 'NCell', 'Colorcell', 'Suite2P', ...
                      'Assembly', 'Correl', 'Sequence', 'PCs', 'Commentaire'};
                  
    btn_save.Enable = 'on';
end

function sauvegarder_et_supprimer(uit)
    global path_save;
    NouvelleTable = uit.Data;
    nb_sessions = height(NouvelleTable);
    
    idx_a_supprimer = find(NouvelleTable.Supprimer == true);
    nb_a_supprimer = length(idx_a_supprimer);
    
    if nb_a_supprimer > 0
        msg = sprintf('⚠️ Vous allez SUPPRIMER DÉFINITIVEMENT %d session(s).\nContinuer ?', nb_a_supprimer);
        reponse = uiconfirm(gcbf, msg, 'Confirmation', 'Options', {'Oui', 'Annuler'}, 'DefaultOption', 2);
        if strcmp(reponse, 'Annuler'), return; end
    end
    
    h_wait2 = waitbar(0, 'Sauvegarde...');
    database = matfile(path_save, 'Writable', true);
    nb_modif = 0;
    
    for i = 1:nb_sessions
        nom_var = char(NouvelleTable.Variable_Name(i));
        
        if NouvelleTable.Supprimer(i) == true
            database.(nom_var) = []; 
            continue;
        end
        
        Session_Actuelle = database.(nom_var);
        if isempty(Session_Actuelle), continue; end
        
        a_modifie = false;
        
        a_modifie = verif_update_string(Session_Actuelle, 'Animal_ID', char(NouvelleTable.Animal_ID(i))) || a_modifie;
        a_modifie = verif_update_string(Session_Actuelle, 'Colorcell_Version', char(NouvelleTable.Colorcell(i))) || a_modifie;
        a_modifie = verif_update_string(Session_Actuelle, 'Suite2P_Version', char(NouvelleTable.S2P_Vers(i))) || a_modifie;
        a_modifie = verif_update_string(Session_Actuelle, 'Commentaire', char(NouvelleTable.Comment(i))) || a_modifie;
        
        a_modifie = verif_update_bool(Session_Actuelle, 'Analysis_Assembly', NouvelleTable.Ana_Assembly(i)) || a_modifie;
        a_modifie = verif_update_bool(Session_Actuelle, 'Analysis_Correlation', NouvelleTable.Ana_Correl(i)) || a_modifie;
        a_modifie = verif_update_bool(Session_Actuelle, 'Analysis_Sequence', NouvelleTable.Ana_Sequence(i)) || a_modifie;
        a_modifie = verif_update_bool(Session_Actuelle, 'Analysis_PCs', NouvelleTable.Ana_PCs(i)) || a_modifie;
        
        if a_modifie
            Session_Actuelle.Animal_ID = char(NouvelleTable.Animal_ID(i));
            Session_Actuelle.Colorcell_Version = char(NouvelleTable.Colorcell(i));
            Session_Actuelle.Suite2P_Version = char(NouvelleTable.S2P_Vers(i));
            Session_Actuelle.Commentaire = char(NouvelleTable.Comment(i));
            
            Session_Actuelle.Analysis_Assembly = NouvelleTable.Ana_Assembly(i);
            Session_Actuelle.Analysis_Correlation = NouvelleTable.Ana_Correl(i);
            Session_Actuelle.Analysis_Sequence = NouvelleTable.Ana_Sequence(i);
            Session_Actuelle.Analysis_PCs = NouvelleTable.Ana_PCs(i);
            
            database.(nom_var) = Session_Actuelle;
            nb_modif = nb_modif + 1;
        end
        waitbar(i / nb_sessions, h_wait2);
    end
    close(h_wait2);
    
    if nb_a_supprimer > 0
        h_net = waitbar(0, 'Nettoyage physique du fichier .mat...');
        d = load(path_save);
        champs = fieldnames(d);
        garder = {}; for c=1:length(champs), if ~isempty(d.(champs{c})), garder{end+1}=champs{c}; end, end
        d_clean = rmfield(d, setdiff(champs, garder));
        save(path_save, '-struct', 'd_clean', '-v7.3');
        close(h_net);
    end
    
    uialert(gcbf, sprintf('✅ Opération terminée !\n- Modifiées : %d\n- Supprimées : %d', nb_modif, nb_a_supprimer), 'Succès');
    btn_save = gcbo; charger_table_editeur(uit, btn_save);
end

function est_modifie = verif_update_string(Struct, champ, new_val)
    est_modifie = false;
    old_val = ''; if isfield(Struct, champ), old_val = Struct.(champ); end
    if ~strcmp(old_val, new_val), est_modifie = true; end
end
function est_modifie = verif_update_bool(Struct, champ, new_val)
    est_modifie = false;
    old_val = false; if isfield(Struct, champ) && Struct.(champ)==true, old_val = true; end
    if old_val ~= new_val, est_modifie = true; end
end