clear; clc; close all;

%% 1. CHEMIN DE LA BASE DE DONNÉES
path_save = 'E:\Data\Aurelie\data\Global_Database.mat';

if ~isfile(path_save)
    error('Le fichier de base de données est introuvable : %s', path_save);
end

%% 2. CHARGEMENT ET RÉSUMÉ
fprintf('Lecture des variables dans la base de données...\n');
liste_variables = who('-file', path_save);
noms_sessions = liste_variables(startsWith(liste_variables, 'S_'));
nb_sessions = length(noms_sessions);

fprintf('%d sessions trouvées. Création du tableau interactif...\n', nb_sessions);

% Pré-allocation
Supprimer  = false(nb_sessions, 1); % NOUVEAU : Colonne logique (checkbox)
Session_ID = strings(nb_sessions, 1);
Condition  = strings(nb_sessions, 1);
Nb_Cells   = NaN(nb_sessions, 1);
Colorcell  = strings(nb_sessions, 1);
S2P_Vers   = strings(nb_sessions, 1);
Comment    = strings(nb_sessions, 1);

h_wait = waitbar(0, 'Chargement...', 'Name', 'Lecture BDD');
for i = 1:nb_sessions
    nom_var = noms_sessions{i};
    data = load(path_save, nom_var);
    Session = data.(nom_var);
    
    Session_ID(i) = string(strrep(nom_var, 'S_', ''));
    if isfield(Session, 'Condition'), Condition(i) = string(Session.Condition); else Condition(i) = ""; end
    if isfield(Session, 'F'), Nb_Cells(i) = size(Session.F, 1); end
    if isfield(Session, 'Colorcell_Version'), Colorcell(i) = string(Session.Colorcell_Version); else Colorcell(i) = ""; end
    if isfield(Session, 'Suite2P_Version'), S2P_Vers(i) = string(Session.Suite2P_Version); else S2P_Vers(i) = ""; end
    if isfield(Session, 'Commentaire'), Comment(i) = string(Session.Commentaire); else Comment(i) = ""; end
    
    waitbar(i / nb_sessions, h_wait);
end
close(h_wait);

% Création de la table avec la colonne "Supprimer" en premier
T = table(Supprimer, Session_ID, Condition, Nb_Cells, Colorcell, S2P_Vers, Comment);

%% 3. CRÉATION DE L'INTERFACE GRAPHIQUE (UI)
f = uifigure('Name', 'Éditeur de Base de Données', 'Position', [100 100 1100 600]);

g = uigridlayout(f, [2 1], 'RowHeight', {'1x', 50});

% Le tableau interactif
% Colonnes : 1=Supprimer (éditable), 2,3,4=Fixes, 5,6,7=Texte (éditables)
col_editable = [true, false, false, false, true, true, true]; 
uit = uitable(g, 'Data', T, 'ColumnEditable', col_editable);
uit.Layout.Row = 1;
uit.Layout.Column = 1;
% Ajustement de la largeur des colonnes pour que la case à cocher soit petite
uit.ColumnWidth = {70, 200, 100, 80, 150, 150, '1x'};

btn = uibutton(g, 'Text', '💾 SAUVEGARDER (Modifications et Suppressions)', ...
               'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', [0.9 0.9 0.9]);
btn.Layout.Row = 2;
btn.Layout.Column = 1;

% Callback qui lance la sauvegarde
btn.ButtonPushedFcn = @(btn,event) sauvegarder_et_supprimer(uit, noms_sessions, path_save);

%% 4. FONCTION DE SAUVEGARDE ET SUPPRESSION
function sauvegarder_et_supprimer(uit, noms_sessions_originaux, path_save)
    NouvelleTable = uit.Data;
    nb_sessions = height(NouvelleTable);
    
    % --- 1. Identifier les sessions à supprimer ---
    idx_a_supprimer = find(NouvelleTable.Supprimer == true);
    nb_a_supprimer = length(idx_a_supprimer);
    
    % Demander confirmation si l'utilisateur veut supprimer des choses
    if nb_a_supprimer > 0
        msg = sprintf('⚠️ Vous êtes sur le point de SUPPRIMER DÉFINITIVEMENT %d session(s) de la base de données.\n\nContinuer ?', nb_a_supprimer);
        reponse = uiconfirm(gcbf, msg, 'Confirmation de suppression', ...
                            'Options', {'Oui, supprimer', 'Annuler'}, ...
                            'DefaultOption', 2, 'CancelOption', 2, 'Icon', 'warning');
        if strcmp(reponse, 'Annuler')
            return; % On annule tout
        end
    end
    
    h_wait2 = waitbar(0, 'Traitement en cours...', 'Name', 'Écriture BDD');
    database = matfile(path_save, 'Writable', true);
    nb_modif = 0;
    
    % --- 2. Boucle de traitement (Modifications + Suppressions) ---
    for i = 1:nb_sessions
        nom_var = noms_sessions_originaux{i};
        
        % A. SUPPRESSION
        if NouvelleTable.Supprimer(i) == true
            % Il n'existe pas de commande "delete variable from matfile".
            % L'astuce est de charger la BDD sans cette variable et de recréer le fichier,
            % OU dans Matlab moderne, utiliser l'option de nettoyage par clear (complexe sur matfile).
            % La méthode la plus sûre (et officielle) :
            % on met la variable existante à 'clear' depuis l'espace de travail, et on re-sauvegardera
            % tout à la fin, MAIS ça demande de tout charger.
            
            % ASTUCE RAPIDE MATLAB : On ne peut pas facilement effacer via 'matfile'.
            % On assigne simplement un tableau vide pour vider le poids du fichier.
            % (La variable existera toujours mais pèsera 0 octet et sera ignorée au prochain chargement)
            database.(nom_var) = []; 
            continue;
        end
        
        % B. MODIFICATIONS
        new_colorcell = char(NouvelleTable.Colorcell(i));
        new_s2p = char(NouvelleTable.S2P_Vers(i));
        new_comment = char(NouvelleTable.Comment(i));
        
        Session_Actuelle = database.(nom_var);
        
        % Si la session avait été vidée précédemment, on l'ignore
        if isempty(Session_Actuelle)
            continue;
        end
        
        a_ete_modifie = false;
        
        old_colorcell = ''; if isfield(Session_Actuelle, 'Colorcell_Version'), old_colorcell = Session_Actuelle.Colorcell_Version; end
        if ~strcmp(old_colorcell, new_colorcell)
            Session_Actuelle.Colorcell_Version = new_colorcell; a_ete_modifie = true;
        end
        
        old_s2p = ''; if isfield(Session_Actuelle, 'Suite2P_Version'), old_s2p = Session_Actuelle.Suite2P_Version; end
        if ~strcmp(old_s2p, new_s2p)
            Session_Actuelle.Suite2P_Version = new_s2p; a_ete_modifie = true;
        end
        
        old_comment = ''; if isfield(Session_Actuelle, 'Commentaire'), old_comment = Session_Actuelle.Commentaire; end
        if ~strcmp(old_comment, new_comment)
            Session_Actuelle.Commentaire = new_comment; a_ete_modifie = true;
        end
        
        if a_ete_modifie
            database.(nom_var) = Session_Actuelle;
            nb_modif = nb_modif + 1;
        end
        
        waitbar(i / nb_sessions, h_wait2);
    end
    
    close(h_wait2);
    
    % --- 3. Nettoyage final du fichier ---
    if nb_a_supprimer > 0
        % Pour vraiment effacer les variables (pas juste les vider), on doit utiliser load puis save
        nettoyer_fichier_physique(path_save);
    end
    
    uialert(gcbf, sprintf('✅ Opération terminée !\n- %d session(s) modifiée(s)\n- %d session(s) supprimée(s)', nb_modif, nb_a_supprimer), 'Succès', 'Icon', 'success');
    
    % Optionnel : fermer l'UI après sauvegarde
    % delete(gcbf); 
end

%% Fonction utilitaire pour nettoyer physiquement le fichier .mat
function nettoyer_fichier_physique(path_save)
    h = waitbar(0, 'Nettoyage du fichier en cours (cela peut prendre 1 minute)...');
    % On lit tout
    data_all = load(path_save);
    champs = fieldnames(data_all);
    
    % On garde uniquement ceux qui ne sont PAS vides (ceux qu'on a mis à [] plus haut)
    champs_a_garder = {};
    for c = 1:length(champs)
        if ~isempty(data_all.(champs{c})) || strcmp(champs{c}, 'DB_Info')
            champs_a_garder{end+1} = champs{c};
        end
    end
    
    % On crée une nouvelle structure propre
    data_clean = rmfield(data_all, setdiff(champs, champs_a_garder));
    
    % On écrase l'ancien fichier
    save(path_save, '-struct', 'data_clean', '-v7.3');
    close(h);
end