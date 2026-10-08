# Documentación del proyecto y explicación de resultados

## 1. Qué se construyó y qué se obtuvo

El repositorio `gbm-drug-response` contiene un flujo reproducible para auditar cuatro archivos de respuesta a fármacos, construir su diccionario de datos, aplicar control de calidad y evaluar regresores de `LN_IC50` reservando líneas celulares completas.

La primera ejecución real entrenó referencias simples, Elastic Net y Extra Trees con metadatos. La **media por fármaco** obtuvo el menor error en validación de desarrollo y fue el modelo seleccionado. En test alcanzó MAE **0,9389** y R² agrupado **0,7846**. Los algoritmos más complejos no mejoraron esa referencia con las entradas actuales.

También se implementó el entrenamiento molecular con expresión, PCA y Elastic Net/árboles. **No se entrenó con expresión biológica real:** esa matriz no está en los adjuntos. Su funcionamiento se comprobó con datos sintéticos aislados en pruebas de software; esos datos no se incorporaron a los notebooks científicos ni a las métricas publicadas aquí.

La conclusión respaldada es que existe una referencia funcional que captura diferencias medias entre medicamentos. Todavía no se ha demostrado capacidad de predecir la respuesta individual a partir de características moleculares.

## 2. Datos de entrada y procedencia

Los cuatro adjuntos del usuario se copiaron a `data/raw/`, conservando sus bytes. Se verificó igualdad de SHA-256 entre cada adjunto y su copia. El manifiesto registra nombres, tamaños y hashes. Esta comprobación acredita conservación de los adjuntos; no verifica que sus bytes sean idénticos a una descarga oficial de GDSC.

| Archivo | Filas | Columnas | SIDM distintos | DRUG_ID distintos | Nombres de fármaco |
|---|---:|---:|---:|---:|---:|
| GDSC2_fitted_dose_response_27Oct23.xlsx | 242.036 | 16 | 969 | 295 | 286 |
| GBMGliomafiltrado.csv | 12.001 | 16 | 51 | 295 | 286 |
| GBMGliomafiltrado.xlsx | 12.001 | 16 | 51 | 295 | 286 |
| Aislados T98G - Hoja 1.csv | 280 | 16 | 1 | 280 | 271 |

Cada fila contiene resultados de una curva ajustada para una combinación de modelo celular y fármaco. No representa un paciente nuevo ni una dosis individual. El subconjunto contiene 37 modelos etiquetados Glioblastoma y 14 Glioma: 8.796 y 3.205 curvas, respectivamente.

Los 51 SIDM son modelos celulares distintos identificables en estos archivos. No se verificó que procedan de 51 donantes diferentes. Para establecer independencia por paciente hay que incorporar relaciones entre modelos, muestras y donantes; Sanger mantiene esa información en Cell Model Passports [3].

### Relaciones verificadas

La comparación completa por `NLME_CURVE_ID`, normalizando tipos y espacios vacíos y con tolerancia numérica de 1e-12, confirmó:

1. Las versiones CSV y XLSX de GBM/glioma contienen las mismas 12.001 curvas y valores.
2. Todas esas curvas están presentes en GDSC2.
3. Las 280 curvas de T98G están presentes en GBM/glioma.

Por tanto, se usa un solo archivo GBM/glioma para entrenar. El GDSC2 completo se auditó, pero **no se usaron sus 969 líneas para entrenar esta primera ejecución**. No se concatenaron los archivos para aumentar artificialmente el tamaño muestral.

### Duplicados, faltantes y cobertura

No se encontraron duplicados de fila completa, de `NLME_CURVE_ID` ni del par `SANGER_MODEL_ID + DRUG_ID` dentro de cada archivo. Hay nueve nombres de fármaco asociados a dos IDs cada uno: Docetaxel, Selumetinib, Oxaliplatin, Fulvestrant, Uprosertib, GSK343, Acetalax, Dactinomycin y Ulixertinib. Se conservaron los IDs separados; un nombre compartido no autoriza fusionar experimentos [1].

Después de convertir vacíos y espacios en ausentes, sólo `PUTATIVE_TARGET` contiene faltantes:

| Dataset | Faltantes | Porcentaje de filas |
|---|---:|---:|
| GDSC2 | 27.872 | 11,52% |
| GBM/glioma, en ambos formatos | 1.193 | 9,94% |
| T98G | 38 | 13,57% |

No faltan valores de `LN_IC50`, AUC, RMSE o Z_SCORE. La ausencia de una diana farmacológica no equivale a ausencia de la etiqueta de entrenamiento, y esa anotación no se usó como predictor en el primer modelo.

La cobertura no es uniforme: cada DRUG_ID tiene entre 16 y 51 líneas en GBM/glioma, con mediana de 47. Las combinaciones sin ensayo no se rellenaron con respuestas inventadas.

## 3. Diccionario y elección de X e y

Cada diccionario CSV incluye nombre, definición, rol, unidad, tipo leído, ejemplos, número de valores distintos y faltantes. Las definiciones se apoyan en la documentación de Sanger [1].

| Campo | Significado y uso en este proyecto |
|---|---|
| DATASET | Colección de origen; procedencia. |
| NLME_RESULT_ID | Identificador del conjunto de ajustes; no predictor. |
| NLME_CURVE_ID | Identificador de curva; verificación y trazabilidad. |
| CELL_LINE_NAME | Nombre de línea; presentación, no predictor. |
| SANGER_MODEL_ID | Identificador de modelo; unión con expresión y agrupación de validación. |
| CANCER_TYPE | Tipo tumoral; predictor categórico de la referencia con metadatos. |
| DRUG_ID | Identificador de fármaco; predictor para medicamentos ya conocidos. |
| DRUG_NAME | Nombre para presentación; no sustituye al ID. |
| PUTATIVE_TARGET | Diana farmacológica propuesta; no es la etiqueta y. |
| PATHWAY_NAME | Vía anotada para la diana; no actividad de vías medida en una célula. |
| MIN_CONC / MAX_CONC | Concentraciones ensayadas en µM; comprobación de rango. |
| LN_IC50 | Logaritmo natural de la IC50 ajustada en µM; outcome primario. |
| AUC | Resultado alternativo de respuesta; excluido de X. |
| RMSE | Calidad del ajuste experimental; filtro, excluido de X. |
| Z_SCORE | Estandarización de LN_IC50 por fármaco; excluido de X por derivarse del outcome. |

En el notebook 02: `X = [DRUG_ID, CANCER_TYPE]`; `y = LN_IC50`. Es regresión, no clasificación. No se transforma otra vez el logaritmo ni se define sensibilidad con un corte global arbitrario.

En el notebook 03: `X = expresión basal de genes` y `y = LN_IC50`, con un modelo separado por DRUG_ID. Así cada línea aparece una vez en cada problema de entrenamiento, sin repetir su expresión cientos de veces al ajustar PCA.

Los identificadores de línea se utilizan para separar conjuntos y unir datos. No se entregan al algoritmo como si fueran características biológicas. Las concentraciones se excluyeron de X por una decisión conservadora del protocolo: evitar que el diseño del ensayo domine la predicción.

## 4. Control de calidad preespecificado

El criterio principal es `0 ≤ RMSE ≤ 0,3`. La documentación de GDSC indica que excluye curvas por encima de 0,3 antes de publicarlas [1]. Se añaden comprobaciones de outcome finito, AUC en [0, 1], concentraciones positivas y ordenadas e identificadores esenciales completos. Las 12.001 curvas pasan esos controles.

| RMSE máximo | Curvas conservadas | Curvas excluidas | Líneas conservadas |
|---|---:|---:|---:|
| 0,15 | 11.327 | 674 | 51 |
| 0,20 | 11.839 | 162 | 51 |
| 0,25 | 11.968 | 33 | 51 |
| 0,30, principal | 12.001 | 0 | 51 |

Las filas de esta tabla son un análisis descriptivo de retención. No se entrenaron cuatro variantes para escoger el umbral con mejor resultado en test. El control principal no mejora por sí solo la calidad: confirma que estos adjuntos cumplen los límites definidos.

### IC50 extrapoladas: una limitación importante

Se comparó `LN_IC50` con `ln(MIN_CONC)` y `ln(MAX_CONC)`:

- 9.181 curvas, **76,50%**, tienen IC50 ajustada por encima del máximo ensayado.
- 36, **0,30%**, quedan por debajo del mínimo.
- 2.784, **23,20%**, quedan dentro del intervalo, incluyendo sus límites.

Estos son valores ajustados fuera del rango experimental, no mediciones directas de una IC50 alcanzada dentro de ese rango. El pipeline los marca, pero no los elimina ni recorta automáticamente; hacerlo podría cambiar sustancialmente la distribución de respuestas.

El modelo actual los trata como valores puntuales de entrenamiento. Las métricas evalúan reproducción de esos ajustes publicados, no precisión contra una medición biológica sin incertidumbre. La alta proporción de extrapolaciones hace prioritaria una revisión de incertidumbre/censura y, eventualmente, del AUC como outcome alternativo bajo un protocolo separado. Un RMSE experimental bajo no resuelve esta limitación.

## 5. División de datos y prevención de fuga

La partición se generó sobre SIDM únicos, sin usar sus valores de respuesta, con `GroupShuffleSplit`, semilla 42 y fracción de test 0,20. El redondeo produce 11 líneas reservadas de 51; no exactamente el 20% de filas. La asignación se conserva en `data/processed/line_split.csv` y se verifica con una firma de IDs, semilla y proporción.

| Partición | Curvas | Líneas | Glioblastoma | Glioma | DRUG_ID |
|---|---:|---:|---:|---:|---:|
| Desarrollo | 9.324 | 40 | 28 líneas | 12 líneas | 295 |
| Test | 2.677 | 11 | 9 líneas | 2 líneas | 295 |

Todas las curvas de una línea permanecen en el mismo conjunto. Los 295 DRUG_ID del test aparecen en desarrollo: esta evaluación comprueba nuevas líneas para fármacos conocidos, no predicción de moléculas nunca observadas.

Dentro de desarrollo se utilizaron cuatro folds `GroupKFold` [4]. Cada fold reserva 10 líneas; sus cantidades de curvas son 2.373, 2.356, 2.285 y 2.310. El entrenamiento de cada fold usa las otras 30 líneas. El modelo se selecciona por MAE promedio entre folds; después se reajusta sobre las 40 líneas de desarrollo y se evalúa en las 11 reservadas.

La codificación categórica, la media por fármaco y, en la ruta molecular, la imputación, selección de genes, escalado y PCA se aprenden exclusivamente en entrenamiento dentro de cada fold [5]. El test no se utiliza para elegir hiperparámetros.

La partición no fue estratificada por tipo tumoral y sólo incluye dos líneas Glioma en test. No permite conclusiones sólidas sobre desempeño en ese subtipo. Reservar SIDM distintos tampoco garantiza independencia por donante si hay modelos emparentados.

**T98G, SIDM01171, quedó en test** con sus 280 curvas. No se incluyó en el ajuste del modelo, pero es una evaluación interna del mismo dataset. No constituye validación externa de otra institución, ensayo o cohorte.

## 6. Modelos entrenados y selección

| Modelo | Qué aprende | Búsqueda realizada |
|---|---|---|
| Media global | Un único promedio de LN_IC50 para todas las curvas de entrenamiento. | Sin hiperparámetros. |
| Media por fármaco | Un promedio de LN_IC50 por DRUG_ID. | Sin hiperparámetros; se recalcula en cada fold. |
| Elastic Net con metadatos | Regresión regularizada sobre DRUG_ID y tipo tumoral codificados one-hot. | alpha: 0,001 / 0,01; l1_ratio: 0,1 / 0,5. |
| Extra Trees con metadatos | Árboles aleatorizados sobre la misma codificación, con posibles interacciones. | 120 árboles; min_samples_leaf: 5 / 15; max_features: 1,0. |

Elastic Net se ejecutó con hasta 20.000 iteraciones y tolerancia 1e-4. Extra Trees usa semilla 42 y dos trabajos. La búsqueda de hiperparámetros se ejecuta con un trabajo y limita los pools numéricos a dos hilos. La ejecución fue en CPU; no requirió GPU, base vectorial ni servidor de base de datos.

Elastic Net eligió alpha 0,001 y l1_ratio 0,5. Extra Trees eligió cinco muestras mínimas por hoja. Esas elecciones corresponden a las opciones ensayadas; no son óptimos universales.

| Modelo | MAE medio en CV | Desviación entre folds |
|---|---:|---:|
| Media por fármaco | **0,9188** | 0,0562 |
| Extra Trees con metadatos | 0,9380 | 0,0609 |
| Elastic Net con metadatos | 0,9945 | 0,0568 |
| Media global | 2,0400 | 0,0257 |

La desviación describe variación entre cuatro folds; no es un intervalo de confianza. Se seleccionó la media por fármaco antes de mirar las métricas de test. El archivo `models/metadata_baseline.joblib` contiene **ese modelo ganador**, sus features declaradas y la configuración de ejecución. Los otros modelos se entrenaron y evaluaron, pero sus objetos ajustados no quedaron persistidos individualmente; sus predicciones y resultados sí se guardaron.

## 7. Resultados de test y su significado

Todas las siguientes métricas se calcularon sobre las mismas 2.677 curvas de las 11 líneas reservadas:

| Modelo | MAE | RMSE predictivo | R² agrupado | Spearman agrupado |
|---|---:|---:|---:|---:|
| Media global | 2,0803 | 2,7732 | -0,000023 | No definido |
| Media por fármaco | **0,9389** | **1,2869** | **0,7846** | 0,8211 |
| Elastic Net con metadatos | 1,0206 | 1,3540 | 0,7616 | 0,8225 |
| Extra Trees con metadatos | 0,9571 | 1,3084 | 0,7774 | 0,8156 |

### MAE: distancia media al valor ajustado

MAE promedia `abs(LN_IC50 observado − LN_IC50 predicho)`. El valor 0,9389 está en escala logarítmica natural: no son 0,9389 µM ni un porcentaje de error.

Para orientar la escala, `exp(0,9389) ≈ 2,56`. Matemáticamente, equivale a la media geométrica del factor de discrepancia simétrico `max(IC50_predicha/IC50_observada, IC50_observada/IC50_predicha)`. No significa que cada predicción tenga un error de 2,56 veces ni da un intervalo de predicción individual.

Frente a la media por fármaco, el MAE de Elastic Net fue **8,70% mayor** y el de Extra Trees **1,94% mayor**. Son diferencias descriptivas de esta partición y búsqueda; no se realizó una prueba de superioridad estadística entre todos los modelos. No demuestran que esos algoritmos carezcan de utilidad al incorporar expresión u otras variables.

### Los dos RMSE no se comparan

El campo `RMSE` original mide ajuste de la curva dosis–respuesta y se usa para QC. El `prediction_rmse` de 1,2869 mide errores del predictor de LN_IC50. Tienen escalas y funciones diferentes: no se debe comparar 1,2869 con el umbral 0,3. El RMSE predictivo penaliza errores grandes más intensamente que MAE.

### R² no es porcentaje de aciertos

R² compara la suma de errores cuadrados con predecir la media de los valores observados en test. Un R² de 0,7846 indica una reducción del 78,46% en esa suma frente a esa referencia matemática, **no que se acierte en el 78,46% de pacientes o medicamentos**.

Aquí se consigue usando sólo un promedio por medicamento: la predicción de un fármaco es idéntica para todas las líneas. Por tanto, gran parte del resultado agrupado recoge diferencias entre fármacos. No acredita que el modelo sepa qué línea será más sensible a un fármaco concreto.

El R² de la media global es ligeramente negativo porque la media aprendida en desarrollo no coincide exactamente con la media del test. No se trata de un error del cálculo.

### Spearman: orden, no calibración

Spearman mide concordancia de rangos. Elastic Net obtuvo una correlación agrupada ligeramente mayor que la media por fármaco, pero un MAE mayor. Ordenar respuestas y aproximar sus valores numéricos son objetivos diferentes; el criterio preespecificado fue MAE, por eso no se cambió de ganador según la correlación.

Spearman de la media global no está definido porque todas sus predicciones son iguales. Dentro de cada uno de los 295 fármacos, Spearman de la media por fármaco tampoco está definido: no hay variación de las predicciones entre líneas. Es una evidencia directa de que esa referencia no distingue sensibilidad individual dentro del mismo medicamento.

### Promedios por curva, por fármaco y por línea

Para la media por fármaco:

| Agregación | MAE | Qué recibe el mismo peso |
|---|---:|---|
| Por curva | 0,9389 | Cada una de las 2.677 curvas. |
| Macro por fármaco | 0,9394 | Cada uno de los 295 DRUG_ID. |
| Macro por línea | 0,9810 | Cada una de las 11 líneas. |

Los valores difieren porque la cobertura no es uniforme. En test cada fármaco tiene respuestas de entre 3 y 11 líneas, lo que limita la precisión de sus métricas individuales.

El bootstrap usa 1.000 remuestreos de las 11 líneas con reemplazo y calcula el MAE macro por línea. Su intervalo percentil exploratorio del 95% es **0,7523–1,3557**, correspondiente a 0,9810, no al MAE por curva 0,9389. Se condiciona al modelo ya ajustado; no incluye la variación de repetir la selección y el entrenamiento.

En `run.json`, la diferencia bootstrap del ganador contra la media por fármaco es cero, con intervalo [0, 0], porque el ganador es esa misma referencia. No significa ausencia de incertidumbre ni una comparación independiente favorable al modelo.

### Heterogeneidad entre líneas reservadas

| Línea | Tipo | Curvas | MAE de la media por fármaco |
|---|---|---:|---:|
| GAMG | Glioblastoma | 280 | 0,6538 |
| M059J | Glioblastoma | 277 | 0,6550 |
| SNB75 | Glioblastoma | 180 | 0,6814 |
| T98G | Glioblastoma | 280 | 0,6899 |
| YH-13 | Glioblastoma | 280 | 0,7666 |
| Hs-683 | Glioma | 179 | 0,8503 |
| D-247MG | Glioblastoma | 180 | 0,9048 |
| 42-MG-BA | Glioblastoma | 279 | 0,9105 |
| AM-38 | Glioblastoma | 279 | 0,9396 |
| LN-405 | Glioblastoma | 279 | 1,0447 |
| NMC-G1 | Glioma | 184 | 2,6945 |

NMC-G1 presenta un error mucho mayor. Esto merece una investigación descriptiva de respuesta, cobertura y procedencia; **no justifica excluirla por empeorar el resultado**. Tampoco demuestra por sí mismo un error de anotación ni explica su causa biológica. La composición de fármacos también difiere entre líneas, así que esta tabla no aísla un efecto biológico puro de cada línea.

## 8. Ruta molecular implementada y pendiente

La entrada requerida consiste en `data/external/expression.csv` y `expression_metadata.json`. El CSV debe contener SIDM únicos y genes numéricos con prefijo `gene_`; los metadatos registran fuente, versión, escala, sistema de IDs, procesamiento y SHA-256. Se admite RMA log2 o log2(TPM+1); no counts crudos, infinitos, muestras totalmente vacías ni IDs repetidos. El formato completo está en `docs/expresion.md`.

El flujo implementado es:

1. Unir expresión y respuestas por SIDM con validación many-to-one; registrar líneas sin coincidencia.
2. Reutilizar la partición de líneas existente, sin reasignar según cobertura molecular.
3. Seleccionar hasta tres DRUG_ID por mayor cobertura en desarrollo, con desempate por ID, o usar una lista preespecificada.
4. Requerir al menos 25 líneas de desarrollo y 5 de test por fármaco. Son mínimos operativos; no un cálculo de potencia.
5. En cada fold, imputar medianas, escoger hasta 2.000 genes variables, estandarizar, calcular 5/10 componentes PCA y estandarizar esos componentes.
6. Comparar la media del fármaco, Elastic Net y Extra Trees; seleccionar por CV y evaluar en test.

Elastic Net molecular explora alpha 0,1/1,0 y l1_ratio 0,1/0,5. Los árboles moleculares usan 150 estimadores y hojas mínimas de 3/5 muestras. PCA limita su dimensión con las muestras y genes disponibles en el fold. Una selección de genes de alta varianza es una reducción estadística, no identificación confirmada de biomarcadores.

Se implementó **PCA**, no actividad de vías. Para trabajar con vías hacen falta un conjunto de genes por vía, una versión y un método de puntuación reproducible. No se utiliza `PATHWAY_NAME` como sustituto.

El estado real guardado es `missing_expression`. Se localizó una distribución histórica oficial de RNA-seq de 940.707.972 bytes, pero no se descargó ni transformó. Aún hay que escoger una fuente/versión, verificar escala y correspondencia de perfiles. La documentación de expresión de Sanger sirve como punto de partida [2].

## 9. Organización del repositorio

| Ubicación | Función |
|---|---|
| `README.md` | Entrada al proyecto y comandos iniciales. |
| `config.json` | Outcome, grupo, QC, semilla, tamaños de partición y límites moleculares. |
| `pyproject.toml`, `uv.lock` | Dependencias declaradas y resolución reproducible. |
| `notebooks/01_auditoria_diccionario.ipynb` | Auditoría, diccionarios, QC, cobertura y partición. |
| `notebooks/02_baseline_metadatos.ipynb` | Entrenamiento y evaluación reales con los adjuntos. |
| `notebooks/03_expresion_pca_modelos.ipynb` | Validación de entrada molecular y entrenamiento condicionado a su disponibilidad. |
| `src/gbm_drug_response/data.py` | Esquema, lectura, diccionario, comparaciones y QC. |
| `src/gbm_drug_response/modeling.py` | Particiones persistidas, folds, modelos de metadatos, métricas y bootstrap. |
| `src/gbm_drug_response/molecular.py` | Contrato de expresión, reducción PCA y modelos por fármaco. |
| `scripts/run_notebooks.py` | Ejecución secuencial, kernel temporal y conservación de outputs. |
| `tests/test_protocol.py` | Nueve comprobaciones del protocolo y su implementación. |
| `data/raw/` | Copias de los cuatro adjuntos. |
| `data/processed/` | Tablas con banderas/QC y partición persistida. |
| `data/external/` | Ubicación prevista de expresión real. |
| `reports/` | Evidencias de auditoría, búsqueda, evaluación y estado molecular. |
| `models/metadata_baseline.joblib` | Media por fármaco ajustada sólo con desarrollo. |
| `docs/` | Este informe, protocolo, resultado resumido y contrato molecular. |

Los notebooks contienen las explicaciones y sus outputs. La lógica reutilizable se mantiene en módulos para compartir exactamente los mismos controles entre ejecución interactiva, ejecución automática y pruebas.

## 10. Cómo reproducir y consultar resultados

En el equipo actual:

```bash
cd /home/saro/research/gbm-drug-response
source .venv/bin/activate
jupyter lab
```

Para repetir la ejecución completa sin interfaz:

```bash
.venv/bin/python scripts/run_notebooks.py
.venv/bin/python -m pytest -q
```

En otro entorno con uv y Python compatible, `uv sync --frozen` reconstruye dependencias desde el lock. También hay que transferir los cuatro archivos originales a `data/raw/`: Git no los incluye. La ejecución comprobada usó Python 3.11.15, pandas 2.3.3, NumPy 2.4.6 y scikit-learn 1.7.2; las versiones de cada corrida de metadatos se guardan en `run.json`.

El orden recomendado es 01 → 02 → 03. El notebook 01 produce los diccionarios y tablas de auditoría; 02 realiza el entrenamiento; 03 completa su ejecución mostrando el estado de ausencia cuando no existe expresión. Ese final normal no equivale a haber entrenado un modelo molecular.

El ejecutor usa un kernel temporal del Python activo, admite nombres de notebooks como argumentos y conserva outputs incluso si una celda falla. Cada celda tiene un límite de 900 segundos. Se puede abrir Jupyter desde la raíz o desde `notebooks/`.

### Evidencia para cada cifra

| Archivo de salida | Pregunta que permite responder |
|---|---|
| `reports/audit/manifest.json` | ¿Qué archivos exactos se analizaron y cuáles son sus hashes? |
| `reports/audit/inventory.csv` | ¿Cuántas curvas, líneas, fármacos y duplicados existen? |
| `reports/audit/dictionary_*.csv` | ¿Qué significa cada columna y cuántos faltantes tiene? |
| `reports/audit/relationships.json` | ¿Se comprobaron las relaciones entre adjuntos? |
| `reports/audit/quality_counts.csv` | ¿Qué criterios QC fallan y cuántas IC50 están fuera de rango? |
| `reports/audit/rmse_sensitivity.csv` | ¿Cuántas curvas conservaría cada umbral descriptivo? |
| `reports/audit/drug_coverage.csv`, `cell_lines.csv` | ¿Cómo se distribuye la cobertura? |
| `data/processed/line_split.csv` | ¿Qué SIDM quedaron en cada conjunto? |
| `reports/metadata_baseline/cv_assignments.csv` | ¿Qué curvas y líneas pertenecen a cada fold? |
| `reports/metadata_baseline/cv_search.csv`, `cv_summary.csv` | ¿Qué candidatos se probaron y cómo se eligió el modelo? |
| `reports/metadata_baseline/test_predictions.csv` | ¿Qué predicción produjo cada modelo para cada curva de test? |
| `reports/metadata_baseline/test_metrics.csv`, `test_by_drug.csv`, `test_by_line.csv` | ¿Cuáles son las métricas y su heterogeneidad? |
| `reports/metadata_baseline/run.json` | ¿Qué configuración, versión, entrada y bootstrap se usaron? |
| `reports/molecular/status.json` | ¿Se entrenó el modelo molecular o falta un insumo? |

Los datos, modelos, reportes y entorno están excluidos mediante `.gitignore`; los notebooks con outputs, código, configuración y documentación pueden versionarse. Los resultados locales no están archivados automáticamente en Git. Durante la revisión se observó el commit `36849c7` (`wip`) con la implementación. Esta ampliación documental no creó commits ni publicó cambios remotos.

Como respaldo documental versionable, `docs/evidencia_primera_ejecucion.json` conserva una instantánea de hashes de entrada, inventario, QC, relaciones, asignación de las 51 líneas, configuración, versiones, búsqueda resumida y métricas de esta primera ejecución. No reemplaza las predicciones individuales ni los archivos originales, que siguen en sus rutas locales y deben conservarse por separado.

Las ejecuciones vuelven a escribir outputs en las mismas rutas; no hay un gestor de experimentos con carpetas únicas. Antes de variar el protocolo, conservar una copia identificada de configuración, reportes, partición y modelos. El código detecta cambios incompatibles de semilla/IDs/proporción en la partición persistida. No cambiar la semilla después de observar test para buscar una cifra mejor.

## 11. Verificación realizada y sus límites

Se ejecutaron los tres notebooks completos y se comprobaron sus 16 celdas de código ejecutadas, sin outputs de error: seis en auditoría, cinco en metadatos y cinco en el notebook molecular. La ruta molecular científica se detuvo funcionalmente en la falta de expresión, informada como estado y no como excepción.

Las nueve pruebas automatizadas verificaron:

1. Límite de RMSE, outcome no finito, concentraciones/AUC inválidas y conservación de extrapolaciones marcadas.
2. Normalización de anotaciones vacías y separación entre diana ausente y outcome disponible.
3. Persistencia de la partición, separación de líneas en folds y detección de cambios manuales de asignación.
4. Media por fármaco aprendida sólo con entrenamiento y fallback para IDs desconocidos.
5. Imputación y PCA sin reajuste al predecir datos de test extremos.
6. Ausencia de métricas/modelo molecular cuando falta expresión.
7. Validación del checksum y rechazo de SIDM repetidos en expresión.
8. Flujo molecular completo con una fixture sintética temporal, incluidos selección, predicción y guardado.
9. Bootstrap basado en medias de líneas, sin dar más peso a una línea sólo por tener más curvas.

También se verificaron compatibilidad de dependencias y conservación de los cuatro originales. Estas pruebas respaldan aspectos concretos del software; no sustituyen validación biológica, cálculo de potencia, auditoría por donante ni validación clínica. No cubren todas las fuentes posibles de sesgo o errores en datos externos.

## 12. Qué falta y cómo avanzar

| Prioridad | Trabajo pendiente | Criterio de finalización |
|---|---|---|
| 1 | Obtener expresión basal y correspondencias oficiales. | Matriz real validada, escala documentada y reporte de SIDM emparejados. |
| 2 | Revisar relaciones de modelos/donantes y origen experimental. | Determinar si es necesario agrupar SIDM relacionados para evitar dependencia entre conjuntos. |
| 3 | Ejecutar el protocolo molecular predefinido. | Comparaciones por fármaco frente a la media sobre las mismas líneas evaluables. |
| 4 | Analizar las IC50 extrapoladas y la incertidumbre del outcome. | Protocolo justificado de tratamiento de extrapolación o análisis alternativo, sin seleccionar por test. |
| 5 | Ampliar datos y validar fuera de la cohorte actual. | Evaluación independiente por modelos, procedencia y, si el objetivo lo requiere, pacientes. |

El holdout de metadatos ya fue observado y se usó para describir errores. Se puede reutilizar para una comparación exploratoria preespecificada con expresión, pero no presentar desarrollos guiados por estos resultados como una nueva validación confirmatoria intacta. Elegir genes, excluir líneas o ajustar umbrales según este test requeriría una evaluación adicional independiente.

No existe aún evidencia de que el modelo identifique tratamientos eficaces para pacientes, generalice a nuevos compuestos o que 51 líneas sean suficientes para un modelo transcriptómico robusto. Tampoco se ha demostrado que los datos moleculares mejorarán las métricas: esa es la siguiente hipótesis a evaluar.

## 13. Fuentes metodológicas

Las cifras de este informe proceden de los adjuntos y de los reportes locales; no de cifras generales publicadas en páginas web. Las decisiones y definiciones se apoyan en:

1. [Sanger: Drug Sensitivity, columnas, QC y múltiples IDs](https://depmap.sanger.ac.uk/documentation/datasets/drug-sensitivity/).
2. [Sanger: Expression, procesamiento y escalas](https://depmap.sanger.ac.uk/documentation/datasets/expression/).
3. [Sanger: Model Relationships, modelos y donantes relacionados](https://depmap.sanger.ac.uk/documentation/cell-models/model-relationships/).
4. [Scikit-learn: GroupKFold](https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.GroupKFold.html).
5. [Scikit-learn: prevención de fuga de información](https://scikit-learn.org/stable/common_pitfalls.html).
6. [Scikit-learn: Elastic Net](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.ElasticNet.html).
7. [Sanger: finalidad del recurso GDSC](https://www.sanger.ac.uk/tool/gdsc-genomics-drug-sensitivity-cancer/).

El protocolo inicial está en `docs/metodologia.md`; los detalles del formato de expresión en `docs/expresion.md`. Las cifras aquí fijadas documentan la primera ejecución, por lo que futuras corridas con cambios de datos o configuración deben generar una interpretación nueva e identificar el experimento correspondiente.
