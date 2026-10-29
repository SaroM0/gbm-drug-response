# Primer resultado ejecutado

> Registro histórico de la primera ejecución. El estado vigente está en [resultados_actuales.md](resultados_actuales.md) y [correcciones de auditoría](correcciones_auditoria_2026-10-29.md); las menciones a expresión ausente o salidas sobrescritas describen aquella versión.

La explicación detallada de todas las métricas, el inventario de implementación y las instrucciones de reproducción están en el [informe completo](informe_completo.md).

Ejecución con los cuatro adjuntos originales y los parámetros iniciales de config.json. Los notebooks conservan los outputs. Los resultados completos y el hash de entrada están en reports/metadata_baseline/.

## Auditoría

GDSC2: 242.036 curvas, 969 SIDM, 295 DRUG_ID y 286 nombres de fármaco. Subconjunto GBM/glioma: 12.001 curvas y 51 SIDM (37 glioblastoma, 14 glioma). T98G: 280 curvas de un solo SIDM, ya incluido en el subconjunto. CSV/XLSX del subconjunto coinciden después de normalizar tipos y vacíos; T98G y GBM/glioma están contenidos en GDSC2 con valores coincidentes. No se detectaron duplicados de fila, curva ni par SIDM–DRUG_ID.

El filtro principal RMSE ≤ 0.3 conserva las 12.001 curvas. Los filtros descriptivos 0.15/0.20/0.25 conservarían 11.327/11.839/11.968, respectivamente. No se entrenaron variantes para seleccionar esos cortes.

Se marcaron 9.181 curvas (76,50%) con IC50 ajustada por encima de la concentración máxima ensayada y 36 (0,30%) por debajo del mínimo. Son extrapolaciones conservadas como valores puntuales por el modelo inicial; un buen RMSE de curva no elimina esa incertidumbre. Su tratamiento debe revisarse metodológicamente sin escoger reglas por las métricas de test.

## Referencia de metadatos

40 líneas en desarrollo y 11 reservadas para test. Test contiene 2.677 curvas; todos sus DRUG_ID aparecen en desarrollo. La búsqueda por cuatro folds agrupados seleccionó la media por fármaco antes de evaluar test.

| Modelo | MAE CV de desarrollo | MAE test | RMSE predictivo test | R² test agrupado |
|---|---:|---:|---:|---:|
| Media global | 2.0400 | 2.0803 | 2.7732 | -0.0000 |
| Media por fármaco | 0.9188 | 0.9389 | 1.2869 | 0.7846 |
| Elastic Net con metadatos | 0.9945 | 1.0206 | 1.3540 | 0.7616 |
| Extra Trees con metadatos | 0.9380 | 0.9571 | 1.3084 | 0.7774 |

Los errores están en escala LN_IC50. El MAE dando igual peso a cada línea para la media por fármaco fue 0.9810, con intervalo bootstrap exploratorio de 0.7523–1.3557. El MAE por curva es diferente porque las líneas tienen distintas cantidades de curvas.

El R² agrupado de 0.7846 se alcanza con una predicción constante para cada medicamento. Por tanto, no demuestra capacidad de predecir diferencias individuales entre líneas. Los modelos con metadatos no mejoraron la referencia en esta ejecución. Estos resultados no justifican cambiar el holdout o buscar parámetros sobre él.

### Cómo interpretar las cifras

- MAE 0,9389 es distancia media en logaritmo natural de IC50, no µM ni porcentaje de error. Su exponencial, aproximadamente 2,56, es un factor de discrepancia geométrico medio; no describe cada predicción ni proporciona un intervalo individual.
- RMSE predictivo 1,2869 mide error al predecir LN_IC50. Es distinto del RMSE experimental usado para filtrar con 0,3.
- R² 0,7846 no significa «78,46% de aciertos». El resultado agrupado captura diferencias entre fármacos; la media por fármaco no distingue líneas dentro de un mismo medicamento.
- Elastic Net tuvo un MAE de test 8,70% mayor que la referencia, y Extra Trees 1,94% mayor. Son comparaciones descriptivas, sin prueba de superioridad estadística.
- Spearman agrupado de Elastic Net fue ligeramente mayor (0,8225 frente a 0,8211), pero su error absoluto también fue mayor. La selección seguía MAE, definido antes de evaluar test.
- El intervalo 0,7523–1,3557 corresponde al MAE macro por línea de 0,9810. La diferencia bootstrap del ganador contra la referencia es cero porque son el mismo modelo, no porque carezca de incertidumbre.

T98G quedó en test con 280 curvas y MAE 0,6899. No participó en el ajuste, pero sigue siendo evaluación interna. NMC-G1 registró MAE 2,6945; el error varía considerablemente entre líneas. No se excluyó ninguna línea por tener un resultado desfavorable. La distinta cobertura de fármacos impide atribuir toda esa variación exclusivamente a biología celular.

## Estado molecular y comprobaciones

Notebook 03 se ejecutó hasta la comprobación de datos y reportó `missing_expression`. El pipeline molecular está implementado, pero no se entrenó con expresión biológica real porque no está en los adjuntos. Su ruta de entrenamiento, evaluación y guardado se verificó con una fixture sintética aislada en los tests; esas métricas no forman parte de este informe.

Los tres notebooks completaron su ejecución sin errores de celdas. Nueve pruebas automatizadas verifican límites de QC, vacíos, separación y persistencia de grupos, media aprendida en entrenamiento, transformaciones sin reajuste en test, contrato molecular, ejecución molecular y bootstrap por línea.

El siguiente insumo para evaluar una mejora biológica es la matriz de expresión basal emparejada por SIDM, con escala y procedencia verificadas. Ver docs/expresion.md.
