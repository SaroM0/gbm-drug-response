# Resultados actuales generados

Generado por `scripts/update_results_docs.py` desde salidas con checksums verificados. Estas ejecuciones son exploratorias; el holdout ya se había observado. Ninguna cifra acredita eficacia de formulaciones o validación clínica.

## Química y priorización

Se resolvieron 243 filas DRUG_ID y 234 nombres, con 232 entidades normalizadas distintas. Hay 152 filas en el informe de revisión; un sinónimo puede ser correcto y seguir requiriendo verificación.

T98G: 214 ajustes por encima del máximo, 0 por debajo del mínimo y 66 dentro del rango.

| stage | drug_id_rows | distinct_nonmissing_compound_ids |
| --- | --- | --- |
| t98g_qc_assays | 280 | 217 |
| structure_accepted | 223 | 213 |
| all_assays_in_range | 50 | 49 |
| within_chosen_loading_scenario | 42 | 41 |
| structural_route_tier_1_to_4 | 24 | 23 |

Lista corta: 23 moléculas provisionales; revisiones químicas aprobadas: 0. Protocolo comparativo listo: False. Los umbrales de carga son escenarios elegidos, no capacidades experimentales universales.

| rank | drug_name | tier | n_assays | ic50_min_um | ic50_max_um | required_loading_pct |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | Mycophenolic acid | 1 | 1 | 3.097493 | 3.097493 | 3.969016 |
| 2 | OSI-027 | 1 | 1 | 3.568838 | 3.568838 | 5.802159 |
| 3 | Staurosporine | 3 | 1 | 0.009286 | 0.009286 | 0.017328 |
| 4 | Tanespimycin | 3 | 1 | 0.110534 | 0.110534 | 0.258957 |
| 5 | AZD7762 | 3 | 1 | 0.450612 | 0.450612 | 0.653261 |
| 6 | LMP744 | 3 | 1 | 0.619885 | 0.619885 | 1.1219 |
| 7 | Docetaxel | 4 | 2 | 0.005148 | 0.008231 | 0.026598 |
| 8 | Luminespib | 4 | 1 | 0.026188 | 0.026188 | 0.048768 |
| 9 | Vinorelbine | 4 | 1 | 0.028617 | 0.028617 | 0.089165 |
| 10 | Dinaciclib | 4 | 1 | 0.09581 | 0.09581 | 0.151954 |
| 11 | Vinblastine | 4 | 1 | 0.064951 | 0.064951 | 0.210699 |
| 12 | Gemcitabine | 4 | 1 | 0.348779 | 0.348779 | 0.367195 |
| 13 | Dasatinib | 4 | 1 | 0.205206 | 0.205206 | 0.400575 |
| 14 | Schweinfurthin A | 4 | 1 | 0.191978 | 0.191978 | 0.421369 |
| 15 | Lestaurtinib | 4 | 1 | 0.551422 | 0.551422 | 0.969336 |
| 16 | AZD8055 | 4 | 1 | 0.542214 | 0.542214 | 1.009719 |
| 17 | MK-1775 | 4 | 1 | 1.077443 | 1.077443 | 2.157503 |
| 18 | Telomerase Inhibitor IX | 4 | 1 | 2.70506 | 2.70506 | 4.115543 |
| 19 | Topotecan | 4 | 1 | 2.490768 | 2.490768 | 4.198966 |
| 20 | Pevonedistat | 4 | 1 | 3.131446 | 3.131446 | 5.555549 |
| 21 | Sabutoclax | 4 | 1 | 2.376465 | 2.376465 | 6.661592 |
| 22 | 5-azacytidine | 4 | 1 | 6.852798 | 6.852798 | 6.694005 |
| 23 | Teniposide | 4 | 1 | 2.603722 | 2.603722 | 6.839061 |

La potencia representativa usa el máximo LN_IC50 entre todos los ensayos, conservando IDs, discrepancias y ambos extremos de rango. Ver `reports/candidates/assays.csv`, `dose_scenarios.csv`, `scenario_summary.csv` y `formulation_readiness.csv`.

## Metadatos

| model | n | mae | prediction_rmse | r2 | spearman | macro_drug_mae | macro_line_mae |
| --- | --- | --- | --- | --- | --- | --- | --- |
| global_mean | 2677 | 2.080341 | 2.773233 | -2.3e-05 |  | 2.063203 | 2.102689 |
| drug_mean | 2677 | 0.93891 | 1.286945 | 0.784644 | 0.82109 | 0.939356 | 0.981006 |
| elastic_net_metadata | 2677 | 1.02058 | 1.354041 | 0.761603 | 0.822455 | 1.024406 | 1.058283 |
| extra_trees_metadata | 2677 | 0.95708 | 1.308395 | 0.777405 | 0.815552 | 0.955816 | 0.99939 |

## Expresión real

Matriz preparada: 50 líneas × 37602 genes. Líneas sin perfil: SIDM01021. Se verificó escala TPM cruda mediante sus sumas por perfil y se aplicó log2(TPM+1) una vez. Se conservaron los ausentes; no se imputó antes de CV.

Estado molecular: `trained`. Los fármacos se eligieron por cobertura en desarrollo, no por error de test.

| DRUG_ID | DRUG_NAME | model | n | mae | baseline_mae | delta_ci95_low | delta_ci95_high |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1003 | Camptothecin | pca_elastic_net | 11 | 1.114042 | 1.193835 | -0.28088 | 0.184039 |
| 1007 | Docetaxel | pca_elastic_net | 11 | 1.179656 | 1.43281 | -0.579434 | -0.009575 |
| 1011 | Navitoclax | pca_extra_trees | 11 | 0.920111 | 0.951191 | -0.145029 | 0.078411 |

Se muestran los modelos elegidos por CV, aunque otro tenga menor error de test. Los intervalos bootstrap son condicionales al entrenamiento fijo; tres fármacos y once líneas de test por fármaco no constituyen validación general.

Los desgloses `outcome_diagnostics.csv` informan métricas con el mismo modelo en ajustes dentro/por encima/por debajo de rango. No se selecciona modelo con esos subgrupos ni se convierten extrapolaciones en cotas de censura medidas.

## Ejecuciones

| kind | run_id | status |
| --- | --- | --- |
| chemistry | 20261029T220646789117Z-0264ec65 | complete |
| candidates | 20261029T224440566967Z-2c38f65f | complete |
| metadata_baseline | 20261029T221151543898Z-00c7f9eb | complete |
| molecular | 20261029T221250672376Z-0b91f133 | complete |
| expression_preparation | 20261029T220647968006Z-d34a9c67 | complete |
| audit | 20261029T221101927170Z-bc52ddfa | complete |
