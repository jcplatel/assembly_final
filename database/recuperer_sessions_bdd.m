function sessions_a_traiter = recuperer_sessions_bdd(path_database, filtre_condition, filtre_analyse)
% RECUPERER_SESSIONS_BDD : Interroge la base de données globale pour
% extraire les noms des variables des sessions correspondant aux critères.
%
% UTILISATION :
%   sessions = recuperer_sessions_bdd(path_BDD, 'nocues', 'Analysis_Correlation');
%   sessions = recuperer_sessions_bdd(path_BDD, 'cues', 'Analysis_Assembly');
%   sessions = recuperer_sessions_bdd(path_BDD, 'nocues', ''); % Pas de filtre d'analyse
%   sessions = recuperer_sessions_bdd(path_BDD, '', ''); % Toutes les sessions

    fprintf('\nInterrogation de la base de données globale...\n');
    if ~isfile(path_database)
        error('Fichier BDD introuvable : %s', path_database);
    end

    % 1. Lister toutes les variables
    liste_vars = who('-file', path_database);
    sessions_noms = liste_vars(startsWith(liste_vars, 'S_'));
    
    sessions_a_traiter = {};
    
    % 2. Filtrer en boucle
    for i = 1:length(sessions_noms)
        nom_var = sessions_noms{i};
        
        % On charge juste la structure sans charger ses sous-champs
        % (ça va extrêmement vite, même si F est lourd, who ne lit que l'en-tête)
        S_test = load(path_database, nom_var);
        S = S_test.(nom_var);
        
        % --- Vérification Condition ---
        if nargin >= 2 && ~isempty(filtre_condition)
            is_bonne_condition = isfield(S, 'Condition') && strcmp(string(S.Condition), string(filtre_condition));
        else
            is_bonne_condition = true; % Si on ne filtre pas, c'est OK
        end
        
        % --- Vérification Checkbox d'Analyse ---
        if nargin >= 3 && ~isempty(filtre_analyse)
            % On vérifie si le champ existe et s'il vaut bien 'true'
            is_analyse_checked = isfield(S, filtre_analyse) && S.(filtre_analyse) == true;
        else
            is_analyse_checked = true; % Si on ne filtre pas, c'est OK
        end
        
        % --- Validation Finale ---
        if is_bonne_condition && is_analyse_checked
            sessions_a_traiter{end+1} = nom_var;
        end
    end
    
    % 3. Résumé console
    str_cond = 'toutes conditions';
    if nargin >= 2 && ~isempty(filtre_condition), str_cond = sprintf('condition "%s"', filtre_condition); end
    
    str_ana = 'aucun filtre d''analyse';
    if nargin >= 3 && ~isempty(filtre_analyse), str_ana = sprintf('case "%s" cochée', strrep(filtre_analyse, '_', ' ')); end
    
    fprintf('-> %d sessions trouvées pour [%s] et [%s].\n\n', length(sessions_a_traiter), str_cond, str_ana);
    
    if isempty(sessions_a_traiter)
        warning('Aucune session ne correspond à vos critères.');
    end
end