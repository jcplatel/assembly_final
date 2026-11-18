function resume_session = export_data_brainbow2(resume_K, res, dist, p_detection, p_F0, source_couleur, seuilstat)
% = export_data_brainbow : une ligne de résumé par session (K retenu).
% Changements :
%  - les anciennes colonnes couleur (successrate, color_assembly, p
%    hypergéométriques non corrigées) sont remplacées par les résultats du
%    test principal (p_minp corrigé, hasard attendu, test du champ)
%  - plus d'ajout ligne à ligne : analysis_brainbow_ALL.xlsx est réécrit
%    entièrement par l'agrégation à partir de ces lignes (avant, chaque
%    relance et chaque K ajoutaient des lignes)
%  - analysis_assemblies_pvalues.mat (une ligne par assemblée, ajoutée pour
%    CHAQUE K) est remplacé par agregation\<variante>\resultats_biais_clonal.csv
%  - distances en µm (distance_calculation2)
% resume_K : la ligne d'exportdata2 pour le K retenu (mêmes colonnes qu'avant)

resume_session = resume_K;
testable = res.testable;
resume_session.source_couleur           = source_couleur;
resume_session.assemblees_testables     = sum(testable);
resume_session.assemblees_biaisees      = sum(res.p_minp(testable) < seuilstat);
resume_session.attendu_par_hasard       = sum(res.taux_nul(testable));
resume_session.p_champ                  = res.p_champ(1);
resume_session.z_champ                  = res.z_champ(1);
resume_session.p_detection              = p_detection;
resume_session.p_F0                     = p_F0;
resume_session.mean_distance_assemblies_um = dist.mean_distance_assembly;
resume_session.mean_distance_all_cells_um  = dist.meandistance_gcamp;
noms_couleur = {'R', 'G', 'B', 'RG', 'RB', 'GB'};
for n = 1:6
    resume_session.(['p_regroupement_' noms_couleur{n}]) = dist.pvalue_distance(n);
end
end
