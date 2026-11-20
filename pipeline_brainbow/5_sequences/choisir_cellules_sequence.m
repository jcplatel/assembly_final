function [corr_max, decalage, critere, dans_sequence, ordre] = choisir_cellules_sequence(DF, debuts, fins, courses, ...
    opts, trames_course)
%CHOISIR_CELLULES_SEQUENCE  Cellules de la séquence sur les courses choisies, selon
%   opts.selection (20 nov. 2026) :
%     'otsu'        : decalage_cellules (code d'origine : trames de course
%                     collées, seuil d'Otsu) ; critere = seuil d'Otsu ;
%     'faux_debuts' : decalage_cellules_test (signal continu, test par
%                     cellule contre des faux débuts tirés pendant la
%                     course, FDR) ; critere = p de chaque cellule.
%   opts.max_lag_trames : décalage maximal en trames (fixé par
%   etape5_sequences selon le mode).
if isfield(opts, 'selection') && strcmp(opts.selection, 'faux_debuts')
    [corr_max, decalage, critere, dans_sequence, ordre] = decalage_cellules_test(DF, debuts, courses, ...
        trames_course, opts.max_lag_trames, opts.n_faux_debuts, opts.fdr_sequence);
else
    [corr_max, decalage, critere, dans_sequence, ordre] = decalage_cellules(DF, debuts, fins, courses, ...
        opts.max_lag_trames);
end
end
