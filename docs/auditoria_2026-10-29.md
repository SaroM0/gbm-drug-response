# Auditoría científica y técnica — 29 de octubre de 2026

**Dictamen:** el flujo reproduce sus resultados y sirve como exploración de respuesta a fármacos libres. La priorización para Fe₃O₄@SiO₂-APTES necesita correcciones antes de cerrar candidatos y dosis. Los resultados actuales no demuestran mejora de una nanoformulación ni predicción molecular individualizada. El propio proyecto reconoce correctamente estos dos límites.

Se revisaron los cuatro originales, código, configuración, pruebas, notebooks guardados, resultados, documentación y fuentes primarias. Se reprodujeron auditoría, entrenamiento de metadatos, construcción química desde caché y ranking en una carpeta temporal. No se modificaron código científico, datos, modelos, partición ni resultados existentes. La evidencia numérica queda en [evidencia_auditoria_2026-10-29.json](evidencia_auditoria_2026-10-29.json).

## 1. Hallazgos que afectan decisiones

### A01 · Alta · Pevonedistat se prioriza a partir de un grupo funcional mal identificado

**Ubicación:** `src/gbm_drug_response/chemistry.py:35`, `:199` y `:232`; `reports/candidates/shortlist.csv`.

El patrón de amina primaria excluye ciertos entornos carbonílicos, pero admite nitrógeno unido a azufre sulfonilo. En la estructura real de pevonedistat, la única coincidencia de `primary_aliphatic_amine` es el N de `–O–S(=O)₂–NH₂`: un sulfamato, no una amina alifática. Esto le asigna ruta de amina, nivel 3, posición 7 y una carga aproximada positiva que no se deduce de ese grupo.

La química de desarrollo del compuesto identifica explícitamente el sulfamato y su importancia farmacológica. No debe proponerse como amina intercambiable para acoplamiento. Véase el [trabajo de sus desarrolladores](https://pubs.acs.org/doi/10.1021/bk-2016-1240.ch001) y el [estudio mecanístico original](https://pubmed.ncbi.nlm.nih.gov/20129059/).

**Corrección:** excluir sulfamatos/sulfonamidas de estos patrones, revisar sus coincidencias átomo por átomo y recalcular rutas y banderas. Una posible ruta mediante otro grupo del compuesto requiere revisión independiente; esta auditoría no determina que sea imposible formularlo. Añadir un caso de regresión con la estructura real.

### A02 · Alta · Los umbrales de carga se generalizan desde materiales y mecanismos diferentes

**Ubicación:** `scripts/rank_candidates.py:29`; `docs/quimica_candidatos.md:64`.

La conversión de unidades de `required_loading_pct()` es correcta. La justificación de usar 5% como carga reportada y 10% como frontera de factibilidad es insuficiente para todos los candidatos.

La fuente de 51,86 µg/mg estudia doxorrubicina adsorbida sobre ferrita **dopada con Gd**, sílice y APTES; no establece capacidad universal para otros compuestos o unión covalente. La fuente de 34,6 mg/g estudia una superficie **carboxilada Fe₃O₄@SiO₂-Glu**, cuya interacción con doxorrubicina difiere de una superficie aminada. Los artículos respaldan ejemplos, no un máximo general ni que 10% sea «el doble de lo mejor publicado». [Estudio APTES/Gd](https://pubs.rsc.org/en/content/articlehtml/2026/ra/d6ra00215c), [estudio de superficie carboxilada](https://researchnow.flinders.edu.au/en/publications/modified-green-synthesis-of-fesub3subosub4subsiosub2sub-nanoparti/).

Con liberación parcial, la carga necesaria para entregar la misma concentración libre aumentaría aproximadamente como `carga_actual / fracción_liberada`, bajo los demás supuestos del cribado. Una conjugación estable puede requerir otra interpretación de la exposición. El límite actual también puede descartar formulaciones que mejoren captación o actúen por otro mecanismo: no identifica imposibilidad física general.

**Corrección:** presentar 5% y 10% como escenarios de cribado elegidos por el equipo; analizar estabilidad del ranking frente a carga y liberación; reemplazar los escenarios por mediciones de cada formulación/lote cuando existan. Reservar «viable» para un estado que especifique qué se ha comprobado.

### A03 · Alta · Alcanzar la IC50 a la dosis máxima no garantiza poder estimarla con 5–25 ppm

**Ubicación:** `scripts/rank_candidates.py:47`; `docs/plan_investigacion_t98g_nanoparticulas.md:20`; `docs/quimica_candidatos.md:127`.

El ranking sólo pregunta si el máximo de nanopartícula puede entregar la IC50 libre. Para estimar una curva también hace falta cubrir su transición. Los siguientes son cálculos de escenario, **no dosis prescritas ni predicciones experimentales**: carga de 5% respecto a masa de nanopartícula, liberación completa, PM del archivo y misma potencia libre de GDSC2.

| Candidato | IC50 libre de referencia (µM) | Fármaco equivalente a 5 ppm NP (µM) | A 25 ppm NP (µM) | Consecuencia bajo ese escenario |
|---|---:|---:|---:|---|
| Docetaxel, ID 1007 | 0,00515 | 0,309 | 1,547 | La dosis mínima ya sería unas 60 veces la IC50 |
| Gemcitabina | 0,349 | 0,950 | 4,749 | La dosis mínima ya superaría la IC50 unas 2,7 veces |
| OSI-027 | 3,569 | 0,615 | 3,076 | El máximo no alcanzaría la IC50 |

Por tanto, un «gran margen» de potencia puede dejar toda la serie por encima de la transición. En sentido contrario, OSI-027 depende de una carga superior al escenario de 5%. La carga efectiva y liberación cambiarán estos cálculos, y MTT puede producir otra potencia: precisamente por eso hace falta un piloto de rango por formulación.

**Corrección:** definir concentraciones a partir de la carga medida, verificar que la respuesta rodee 50% y comparar libre/cargado a dosis equivalentes de fármaco. Si no hay cruce, informar IC50 no estimable o límite justificado, sin forzar un valor. Preespecificar tiempo y contraste principal; 24/48/72 h y varios candidatos no son una única comparación.

### A04 · Media · El informe confunde conteos marginales, niveles y moléculas

**Ubicación:** `docs/quimica_candidatos.md:99` y `:107`; `scripts/rank_candidates.py:163`.

En T98G hay **214/280 = 76,43%** de ajustes por encima del máximo, y **66/280** dentro del rango; no 174/280. Los 174 son los registros asignados al nivel 8 después de que el nivel 9 haya absorbido los casos de estructura no resuelta/ambigua. Son preguntas diferentes.

El embudo publicado tampoco es acumulativo: 66 y 55 son conteos marginales. Aplicando los criterios conjuntamente, el flujo real es:

| Paso acumulativo del algoritmo actual | Registros DRUG_ID | `compound_id` no nulos distintos |
|---|---:|---:|
| T98G | 280 | 218 |
| Estructura aceptada por `structure_ok` | 226 | 216 |
| Además, IC50 dentro del rango en estos datos | 52 | 51 |
| Además, carga calculada ≤10% | 43 | 42 |
| Además, ruta de niveles 1–4 | 25 | 24 |

La última deduplicación produce la lista de **24 moléculas**. Las 243 estructuras de la tabla química son filas por DRUG_ID con estructura resuelta, no 243 moléculas distintas. Los `compound_id` de esta tabla conservan además las limitaciones de A07.

**Corrección:** separar explícitamente conteos marginales de un embudo acumulativo y mantener la unidad de conteo en cada columna. Generar las cifras documentales desde las salidas para evitar desincronización.

### A05 · Media · La deduplicación conserva sistemáticamente el ensayo más favorable

**Ubicación:** `scripts/rank_candidates.py:130` y `:145`.

Se ordena por menor carga requerida y luego se aplica `drop_duplicates(compound_id, keep='first')`. Para docetaxel se conservan 5,15 nM del ID 1007; el ID 1819 da 8,23 nM. La razón es 1,60. Las cargas respectivas son 0,0166% y 0,0266%. Ambos valores son compatibles con alta potencia, pero conservar siempre el menor induce optimismo al resumir ensayos repetidos.

`all_gdsc_drug_ids` también se calcula después de filtrar niveles: no garantiza conservar todos los ensayos de una identidad, incluidos los discrepantes o rechazados. En los datos hay diferencias mayores para otros fármacos; los valores extrapolados requieren especial cuidado.

**Corrección:** conservar todos los ensayos por identidad, rangos, QC y procedencia; predefinir cómo se obtiene la potencia representativa o usar un análisis conservador de sensibilidad. No promediar sin considerar rangos y extrapolación. Registrar los IDs excluidos con sus motivos.

### A06 · Media · «Dentro de rango» ignora la concentración mínima

**Ubicación:** `scripts/rank_candidates.py:65`, `:75`, `:80` y `:86`.

El cálculo usa `~above_range`. Una IC50 por debajo del mínimo también queda declarada dentro de rango. Una fixture aislada con `LN_IC50=-10`, `MIN_CONC=0,001` y `MAX_CONC=10` reproduce `t98g_in_range=True` y nivel 1 si cumple estructura/dosis.

**Impacto comprobado:** T98G no contiene ajustes por debajo del mínimo, de modo que este error no cambia sus 66 registros actuales ni demuestra por sí solo un cambio de la lista corta. Sí afecta las fracciones del panel GBM/glioma, que contiene **36** ajustes por debajo del mínimo, y futuras entradas legítimas.

**Corrección:** usar conjuntamente ambas banderas de extrapolación, conservar los tres estados y añadir cobertura de este caso a las pruebas del ranking.

### A07 · Media · La identidad y los grupos funcionales dependen de cómo se representa la sal

**Ubicación:** `src/gbm_drug_response/chemistry.py:33`, `:35`, `:142` y `:274`; `tests/test_chemistry.py`.

El primer bloque del InChIKey de la estructura completa no garantiza agrupar una sal con su molécula parental. La detección de grupos también depende del estado de protonación. Reproducciones con RDKit:

| Par de representaciones | Primer bloque InChIKey | Ruta asignada actualmente |
|---|---|---|
| Ácido acético / acetato de sodio | `QTBSBXVTEAMEQO` / `VMHLLURERBWHNL` | COOH / sin asa |
| Butilamina / cloruro de butilamonio | `HQABUPZFAYXKJW` / `ICXXXLGATNSZAV` | Amina / sin asa |

El test de alias utiliza InChIKeys inventados con el mismo primer bloque y una sal representada como `CCCCN.Cl`; por ello no detecta el fallo con `CCCC[NH3+].[Cl-]`. El primer bloque también omite información estereoquímica relevante para no fusionar isómeros indiscriminadamente. Véanse la [especificación de InChI](https://www.inchi-trust.org/technical-faq/) y las [herramientas de estandarización de RDKit](https://www.rdkit.org/docs/source/rdkit.Chem.MolStandardize.rdMolStandardize.html).

**Corrección:** distinguir identidad de la forma adquirida, entidad activa y familia de conectividad; estandarizar fragmentos/cargas con reglas documentadas que preserven estereoquímica y excepciones para complejos metálicos. Calcular el PM sobre la base de masa que realmente se dosifica. Este contraejemplo demuestra un defecto del método, no que las 24 identidades actuales estén todas mal.

### A08 · Alta para la siguiente fase · La mayor parte del outcome molecular sigue siendo extrapolación

**Ubicación:** `src/gbm_drug_response/data.py:107`; `src/gbm_drug_response/modeling.py:148`; `docs/metodologia.md`.

En GBM/glioma hay **9.181/12.001 = 76,50%** de IC50 ajustadas por encima del máximo, 36 por debajo y 2.784 dentro del rango. El entrenamiento las trata como números puntuales. Es una limitación ya reconocida, pero dominante: reproducir las cifras ajustadas no prueba predecir una IC50 experimental bien identificada. Tampoco debe llamarse automáticamente «cota» al valor numérico extrapolado del archivo.

Sanger confirma que los valores son ajustes y que GDSC2 usa CellTiter-Glo a 72 h. El RMSE experimental bajo evalúa ajuste a los puntos disponibles y no demuestra precisión del cruce extrapolado. La AUC también depende del intervalo de concentraciones. [Documentación primaria de GDSC](https://depmap.sanger.ac.uk/documentation/datasets/drug-sensitivity/).

**Corrección:** definir antes de la siguiente evaluación un análisis principal y sensibilidades que representen observación/extrapolación e incertidumbre; considerar curvas originales y un tratamiento de censura cuando esté justificado. No resolverlo eliminando valores fuera de rango para maximizar el test. El holdout ya observado conserva carácter exploratorio para desarrollos guiados por él.

### A09 · Media · Faltan criterios operativos para atribuir mejora y transporte

**Ubicación:** `docs/plan_investigacion_t98g_nanoparticulas.md:18`, `:301` y `:336`.

El plan reconoce controles de interferencia MTT, dosis equivalentes, liberación, réplicas independientes y el posible brazo libre+NP. Son decisiones pendientes, no evidencia de que ya estén aplicadas. Faltan datos crudos del antecedente de oxaliplatino para auditar sus 15 ppm y una definición operativa del éxito: tiempo principal, magnitud de mejora relevante, incertidumbre aceptable y criterio endotelial.

El riesgo de interferencia requiere verificación sobre **estas** partículas y concentraciones: el estudio citado utilizó otro recubrimiento, línea y rango de exposición; no cuantifica el sesgo de este experimento. [Estudio primario de interferencia](https://pmc.ncbi.nlm.nih.gov/articles/PMC10845919/).

Para Transwell falta especificar qué se mide como transporte —fármaco libre, fármaco total, partícula o marcador— y cómo se distingue paso con barrera íntegra de fuga por daño endotelial. El trabajo original con hCMEC/D3 evalúa propiedades de permeabilidad y transporte; nombrar la línea no valida por sí solo la barrera. [Estudio de caracterización](https://pubmed.ncbi.nlm.nih.gov/19013850/).

**Corrección:** cerrar estos campos en el protocolo del piloto, incorporar controles pertinentes, criterios de integridad/recuperación y unidades equivalentes. Predefinir análisis de réplicas biológicas, placa/lote y resultados no estimables. La mezcla libre+NP responde una pregunta adicional; su necesidad depende de la atribución mecanística que se quiera afirmar.

### A10 · Media · La trazabilidad de nuevas ejecuciones y de la química es incompleta

**Ubicación:** `.gitignore`; `scripts/rank_candidates.py:155`; `src/gbm_drug_response/molecular.py:120`; `docs/expresion.md`.

Los resultados se sobrescriben en rutas comunes. `reports/candidates/run.json` no registra hashes de entradas, tabla química, código o parámetros completos de anotación; tampoco versión de RDKit. Una repetición de la química depende de una caché excluida de Git. La ausencia de expresión tras una corrida molecular previa puede dejar modelos/métricas antiguos junto a un nuevo `missing_expression`, porque esa rama sólo actualiza el estado.

Además, la documentación que afirma que el ZIP no se descargó ya está desactualizada: hay un ZIP de 940.707.972 bytes con matriz TPM cuyos encabezados incluyen **50/51** SIDM. La línea ausente en esos encabezados es `SIDM01021`. Esto no certifica integridad completa del ZIP, escala, calidad de mediciones o cobertura de genes: sólo se inspeccionaron índice y cabeceras. Faltan el procesamiento y los metadatos canónicos que exige el pipeline.

La tabla local `model_list_20260921.csv` sí permite una comprobación adicional: las 51 líneas tienen 51 `patient_id` distintos. Hay dos modelos con padres anotados; ninguno de esos padres está en la cohorte actual. No se encontró cruce por `patient_id` ni por esas dos relaciones en la partición. Es evidencia favorable dentro de esa versión local, cuya autenticidad de descarga no se verificó.

**Corrección:** guardar ejecuciones identificadas con código, configuración, inputs, caché y versiones; asociar todos los archivos a un estado de ejecución para impedir interpretar salidas antiguas como actuales. Actualizar inventario y preparar la expresión con selección de perfiles, transformación y procedencia explícitas. Repetir la auditoría de relaciones al ampliar a pan-cáncer.

## 2. Resultados y controles que sí se sostienen

- Los cuatro hashes originales coinciden con el manifiesto local. CSV y XLSX de GBM/glioma contienen la misma información; GBM/glioma pertenece al GDSC2 aportado y T98G al subconjunto. No se comprobó identidad de bytes frente a una descarga oficial.
- Se reprodujeron 242.036 curvas/969 líneas en el archivo completo y 12.001 curvas/51 líneas en GBM/glioma. Las tablas no se concatenan como cohortes independientes.
- La partición mantiene 40 líneas de desarrollo y 11 de test; las líneas no se cruzan entre conjuntos ni folds. Las transformaciones moleculares se aprenden dentro del pipeline de CV.
- Se reprodujeron búsqueda CV, métricas, predicciones y desgloses por línea/fármaco. La media por fármaco gana entre los candidatos probados: **MAE 0,93891**, RMSE 1,28695 y R² agrupado 0,78464 sobre 2.677 respuestas de test. Ese R² refleja en gran parte diferencias entre medicamentos y no acredita individualización.
- Como comprobación adicional se evaluó la mediana por fármaco: MAE CV 0,92993 frente a 0,91884 de la media; MAE test 0,94465 frente a 0,93891. No cambia la conclusión. Es una comprobación exploratoria de auditoría, no una selección confirmatoria nueva.
- Se reprodujeron la tabla química, el ranking y las 24 identidades de la lista corta usando la caché existente. Esto prueba reproducibilidad computacional, no corrección de los supuestos químicos.
- **23 pruebas pasan**. Los tres notebooks guardados contienen 16 celdas de código ejecutadas y ninguna salida de error. No se reejecutaron los notebooks: se ejecutaron directamente los flujos de código en aislamiento. El entrenamiento molecular científico continúa pendiente; el test sintético comprueba software.

## 3. Decisión sobre los tres candidatos propuestos

| Candidato | Lectura defendible hoy | Condición pendiente principal |
|---|---|---|
| Docetaxel | Potencia libre alta en dos ensayos de T98G; candidato para estudio de formulación | Integrar ambos ensayos y ajustar rango a carga/liberación reales; no asumir que 5–25 ppm permite estimar IC50 |
| Gemcitabina | Potencia libre dentro de rango y carga teórica pequeña | Verificar derivatización, liberación, actividad y rango; ausencia de banderas del script no acredita ausencia de riesgos |
| OSI-027 | Carboxilo detectado y potencia libre dentro de rango | Factibilidad de carga, mecanismo de acción/liberación tras unión y conservación de actividad; no se ha demostrado que el COOH sea su único grupo ionizable ni que medie la unión a la diana |

La frase «tres rutas químicas distintas» del documento no se cumple con esta terna: docetaxel y gemcitabina comparten la categoría de derivatización; OSI-027 aporta la ruta directa por COOH. Pueden diferir sus detalles, pero éstos todavía no se documentaron. Tampoco puede inferirse de una búsqueda fallida en PubChem que los 52 nombres no resueltos sean todos no divulgados o imposibles de comprar.

## 4. Orden de trabajo recomendado

1. **Corregir el cribado:** SMARTS de sulfamatos/sulfonamidas, estados de rango, identidad normalizada, resumen de múltiples ensayos y conteos documentales. Verificar identidades/rutas de los finalistas contra fuentes y forma comercial.
2. **Cerrar el piloto:** contraste principal, carga y liberación medidas, dosis equivalentes y rango, controles de ensayo, réplicas independientes, criterio endotelial y definición de transporte. Registrar éxito/fallo antes de observar resultados.
3. **Congelar una versión de priorización:** escenarios de carga/liberación, ranking y exclusiones, hashes y protocolo. Confirmar prospectivamente los candidatos; con tres moléculas el objetivo inmediato es aprender sobre formulación.
4. **Continuar expresión sólo si responde la pregunta científica elegida:** preparar la matriz existente, auditar correspondencias y extrapolación, y evaluar aporte frente a la media. La expresión no reemplaza las etiquetas de respuesta de nanopartículas ni justifica por sí sola deep learning.

El resultado defendible para la investigación actual es una **priorización exploratoria auditada seguida de validación experimental**. Las probabilidades de éxito de formulaciones, superioridad del cargado y transporte con barrera íntegra siguen siendo hipótesis por contrastar.
