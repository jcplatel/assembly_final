function C = consensus_assemblees(dossier_etape1, liste_K, NCell, opts)
%CONSENSUS_ASSEMBLEES  Assemblées de consensus entre plusieurs K (étape 2,
%   règle 'consensus'). Ne refait aucun clustering : lit les assemblées
%   validées de l'étape 1 (k<K>\assemblees.mat) pour chaque K de liste_K.
%
%   Repris de assembly_core_analysis.m (co-appartenance des cellules aux
%   assemblées validées, arbre 'average', coupe au seuil, co-appartenance
%   moyenne interne >= seuil), avec une autre normalisation, choisie sur
%   simulation (opts.normalisation) :
%     'union' : nb de K où i et j sont dans la même assemblée
%               / nb de K où i OU j est dans une assemblée
%     'tous'  : ... / nb de K (comme le script d'origine). Une assemblée qui
%               n'apparaît qu'aux grands K n'atteint alors jamais le seuil.
%   Simulation (22 sessions de 489 cellules et 370 SCE, K de 4 à 12 ;
%   assemblées retrouvées = Jaccard >= 0,5 ; fausses = Jaccard < 0,2 avec
%   toutes les vraies), par session :
%                             4 implantées   6 implantées   8 implantées   aucune
%     'union' 0,7, min 2 K    2,8 (0 fausse) 4,3 (0)        6,0 (0)        0
%     'tous'  0,7             1,8 (0)        2,5 (0)        4,5 (0,2)      0
%     un seul K ('exces')     2,3 (1,7)      4,8 (1,2)      5,8 (1,0)      0
%     meilleur K a posteriori 2,7 (0,5)      5,0 (0,8)      6,5 (1,0)      0
%   Une paire ensemble dans moins de opts.min_K K compte 0 (une seule
%   rencontre ne suffit pas).
%
%   opts : seuil (fraction, ex. 0.7), min_K, normalisation, taille_min
%   C : assemblees (cellules de chaque assemblée, de la plus grande à la
%       plus petite), consensus_moyen (co-appartenance moyenne interne),
%       K_presence (K où une assemblée de l'étape 1 ressemble à celle-ci :
%       Jaccard >= 0,5), K_utilises

if ~isfield(opts, 'normalisation') || isempty(opts.normalisation), opts.normalisation = 'union'; end
nK = numel(liste_K);
ens = zeros(NCell, NCell);
assign = false(nK, NCell);
asm_K = cell(1, nK);
K_utilises = [];
n = 0;
% (13 nov.) opts.assemblees_K{iK} : assemblées du K n° iK données directement
% (cell de listes de cellules ; [] = ce K n'a pas été calculé), au lieu
% d'être lues dans les fichiers. Utilisé par screening_reglages.
en_memoire = isfield(opts, 'assemblees_K') && ~isempty(opts.assemblees_K);
for iK = 1:nK
    if en_memoire
        if ~iscell(opts.assemblees_K{iK}), continue; end
        S = struct();
        S.assemblyortho = opts.assemblees_K{iK};
    else
        if isfield(opts, 'dossier_analyse') && ~isempty(opts.dossier_analyse)
            f = fichier_assemblees(opts.dossier_analyse, dossier_etape1, liste_K(iK));   % v1 ou v2
        else
            f = fullfile(dossier_etape1, sprintf('k%d', liste_K(iK)), 'assemblees.mat');
        end
        if ~isfile(f), continue; end
        S = load(f, 'assemblyortho');
    end
    n = n + 1;
    K_utilises(end+1) = liste_K(iK); %#ok<AGROW>
    asm_K{n} = {};
    ens_K = false(NCell, NCell);     % une paire compte une fois par K, même si deux
    for a = 1:numel(S.assemblyortho)  % assemblées la partagent (version 2, 11 nov.)
        cel = unique(S.assemblyortho{a}(:));
        cel = cel(cel >= 1 & cel <= NCell);
        if isempty(cel), continue; end
        ens_K(cel, cel) = true;
        assign(n, cel) = true;
        asm_K{n}{end+1} = cel;
    end
    ens = ens + ens_K;
end
asm_K = asm_K(1:n);
assign = double(assign(1:n, :));
C = struct('assemblees', {cell(1, 0)}, 'consensus_moyen', zeros(1, 0), ...
    'K_presence', {cell(1, 0)}, 'K_utilises', K_utilises, 'opts', opts);
if n == 0, return; end

switch opts.normalisation
    case 'tous'
        den = n * ones(NCell, NCell);
    case 'union'
        nb = sum(assign, 1);
        den = nb' + nb - assign' * assign;
    case 'both'
        den = assign' * assign;
    otherwise
        error('consensus : normalisation ''union'', ''tous'' ou ''both''.');
end
M = zeros(NCell, NCell);
ok = den > 0;
M(ok) = ens(ok) ./ den(ok);
M(ens < opts.min_K) = 0;
M(1:NCell + 1:end) = 1;

% Arbre (lien moyen) sur les cellules qui ont au moins un partenaire
autre = M;
autre(1:NCell + 1:end) = 0;
actives = find(any(autre > 0, 2));
if numel(actives) < opts.taille_min, return; end
D = 1 - M(actives, actives);
D = (D + D') / 2;
D(1:numel(actives) + 1:end) = 0;
Z = linkage(squareform(D, 'tovector'), 'average');
lab = cluster(Z, 'Cutoff', 1 - opts.seuil + 1e-9, 'Criterion', 'distance');

for c = unique(lab)'
    cel = actives(lab == c);
    if numel(cel) < opts.taille_min, continue; end
    sm = M(cel, cel);
    moy = mean(sm(triu(true(numel(cel)), 1)));
    if moy < opts.seuil, continue; end
    C.assemblees{end+1} = cel(:)';
    C.consensus_moyen(end+1) = moy;
    Kp = [];
    for iK = 1:n
        for a = 1:numel(asm_K{iK})
            b = asm_K{iK}{a};
            if numel(intersect(cel, b)) / numel(union(cel, b)) >= 0.5
                Kp(end+1) = K_utilises(iK); %#ok<AGROW>
                break
            end
        end
    end
    C.K_presence{end+1} = Kp;
end
% de la plus grande à la plus petite
[~, ordre] = sort(cellfun(@numel, C.assemblees), 'descend');
C.assemblees = C.assemblees(ordre);
C.consensus_moyen = C.consensus_moyen(ordre);
C.K_presence = C.K_presence(ordre);
end
