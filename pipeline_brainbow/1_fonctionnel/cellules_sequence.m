function [dans_sequence, info] = cellules_sequence(F, speed, opts_seq)
%CELLULES_SEQUENCE  Cellules activées en séquence pendant la course (15 nov.).
%   Sert à config.pre.population = 'sequence' / 'hors_sequence' : chercher
%   les SCE et les assemblées dans une partie des cellules seulement
%   (idée de l'utilisateur, 15 nov.).
%
%   Même calcul que l'étape 5 (etape5_sequences, parties 1 et 2, puis
%   decalage_cellules sur TOUTES les courses : c'est « sequence_tout »),
%   recopié sans changement : F divisé par un polynôme de degré 3, filtre de
%   Savitzky-Golay (3, 5), z-score ; courses = vitesse lissée > seuil ;
%   corrélation croisée de chaque cellule avec les débuts de course ;
%   cellule de la séquence = maximum > seuil d'Otsu et décalage < 0.
%   N'utilise que l'activité et la vitesse, jamais les couleurs.
%
%   F : cellules x trames (data_recording.F) ; speed : vitesse ;
%   opts_seq : réglages de options_sequences (défaut : ceux de l'étape 5).
%   dans_sequence : 1 x NCell logique.
%   info : n_courses, debuts, fins, corr_max, decalage, seuil, ordre, statut
%          ('' si tout va bien ; sinon la raison, et dans_sequence est vide).

if nargin < 3 || isempty(opts_seq), opts_seq = options_sequences(struct()); end

%% 1. Activité et vitesse (comme etape5_sequences)
F = double(F);
speed = double(speed(:)');
[NCell, Nz] = size(F);
DF = F;
ws = warning('off', 'all');
for i = 1:NCell
    p0 = polyfit(1:Nz, DF(i, :), 3);
    DF(i, :) = DF(i, :) ./ polyval(p0, 1:Nz);
end
warning(ws);
DF = sgolayfilt(DF', 3, 5)';
DF = normalize(DF, 2, 'zscore');
clear F

% vitesse : même longueur que la fluorescence
if numel(speed) >= Nz
    speed = speed(1:Nz);
else
    speed(end+1:Nz) = NaN;
end

%% 2. Courses (comme etape5_sequences)
vitesse_lissee = smoothdata(speed, 'gaussian', opts_seq.lissage_vitesse);
en_course = double(vitesse_lissee > opts_seq.seuil_vitesse);
changements = diff([0, en_course, 0]);
debuts = find(changements == 1) - 1;
fins = find(changements == -1) - 1;
garder = debuts >= 1;
debuts = debuts(garder);
fins = fins(garder);
n_courses = numel(debuts);

info = struct('n_courses', n_courses, 'debuts', debuts, 'fins', fins, 'corr_max', [], 'decalage', [], ...
    'seuil', NaN, 'ordre', [], 'statut', '');
dans_sequence = false(1, 0);
if n_courses < opts_seq.min_courses
    info.statut = sprintf('%d course(s) : trop peu (minimum %d)', n_courses, opts_seq.min_courses);
    return
end

%% 3. Cellules de la séquence, sur toutes les courses
[info.corr_max, info.decalage, info.seuil, dans_sequence, info.ordre] = ...
    decalage_cellules(DF, debuts, fins, 1:n_courses, opts_seq.max_lag);
dans_sequence = logical(dans_sequence(:)');
end
