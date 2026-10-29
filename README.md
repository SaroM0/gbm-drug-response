# gbm-drug-response

Auditoría reproducible de GDSC2, predicción exploratoria de `LN_IC50` con expresión basal y priorización de fármacos para estudiar formulaciones en T98G. CPU; sin deep learning. Las fuentes externas y sus límites se registran por ejecución.

## Documentación

Empezar por los [resultados actuales generados](docs/resultados_actuales.md) y el [seguimiento de las correcciones de auditoría](docs/correcciones_auditoria_2026-10-29.md). El [informe de la primera ejecución](docs/informe_completo.md) queda como registro histórico.

| Documento | Contenido |
|---|---|
| [Plan de investigación T98G y nanopartículas](docs/plan_investigacion_t98g_nanoparticulas.md) | Definición experimental acordada con el laboratorio (Fe₃O₄@SiO₂-APTES, MTT, 4 brazos), inventario, datos faltantes, inputs/outputs y plan de validación probabilística, ML y deep learning. |
| [Química y priorización de candidatos](docs/quimica_candidatos.md) | Identidad química desde PubChem, anotación de conjugabilidad a la superficie APTES, factibilidad de dosis en ppm y lista corta de candidatos para el ensayo. |
| [Informe completo](docs/informe_completo.md) | Trabajo realizado y lectura detallada de MAE, RMSE, R², Spearman, extrapolaciones y variación entre líneas. |
| [Primer resultado](docs/primer_resultado.md) | Resumen de las métricas y conclusión de la primera ejecución. |
| [Metodología](docs/metodologia.md) | Protocolo preespecificado, validación, QC y fuentes. |
| [Entrada de expresión](docs/expresion.md) | Archivos y metadatos necesarios para ejecutar el entrenamiento molecular real. |
| [Evidencia de la primera ejecución](docs/evidencia_primera_ejecucion.json) | Instantánea versionable de hashes, configuración, partición, QC y métricas. |

Resultado principal: la media por fármaco fue seleccionada en CV y obtuvo MAE 0,9389 en test; los modelos con metadatos no la superaron. El R² agrupado de 0,7846 no acredita predicción individual molecular. El 76,50% de las IC50 ajustadas de GBM/glioma está por encima del máximo ensayado, una limitación relevante del outcome. La expresión real ya se prepara desde el ZIP local y se ejecuta el protocolo molecular; consultar sus resultados exploratorios actuales en la documentación generada.

## Abrir los notebooks

```bash
cd /home/saro/research/gbm-drug-response
source .venv/bin/activate
jupyter lab
```

El entorno ya está instalado. En otro equipo, usar Python 3.11–3.13 y `uv sync --frozen` (o `python3.11 -m venv .venv` y `.venv/bin/pip install -e .` para resolver dependencias sin lock). En VS Code seleccionar `.venv/bin/python` como kernel.

Ejecutar en orden:

| Notebook | Resultado |
|---|---|
| `01_auditoria_diccionario.ipynb` | Diccionario de 16 columnas, dimensiones, faltantes, duplicados, relaciones, cobertura, QC y partición persistida. |
| `02_baseline_metadatos.ipynb` | Entrenamiento real con los adjuntos: referencias de media, Elastic Net y Extra Trees; CV y test por líneas completas. |
| `03_expresion_pca_modelos.ipynb` | Entrenamiento por fármaco con expresión, PCA y Elastic Net/árboles. Requiere la entrada descrita en `docs/expresion.md`; sin ella muestra el motivo y no inventa resultados. |

Para ejecutar los tres sin interfaz:

```bash
.venv/bin/python scripts/run_notebooks.py
.venv/bin/python -m pytest -q
```

El ejecutor guarda outputs dentro de cada notebook, conserva los outputs completados si una celda falla y devuelve error. Se puede indicar uno o varios nombres de notebook como argumentos. Abrir Jupyter desde la raíz o desde `notebooks/` funciona; no se necesitan rutas personales dentro de los notebooks.

## Priorización de candidatos para el ensayo

Dos scripts independientes de los notebooks, para la fase de nanopartículas:

```bash
.venv/bin/python scripts/build_compounds.py   # entidad activa normalizada y caché archivada
.venv/bin/python scripts/rank_candidates.py   # todos los ensayos y escenarios de carga/liberación
```

El primero resuelve los nombres de fármaco contra PubChem y anota grupos funcionales con RDKit; guarda cada respuesta cruda en `data/external/pubchem_cache/`, de modo que repetirlo no vuelve a consultar la red. El segundo cruza esa tabla con la potencia medida en T98G y en el panel de GDSC2, y aplica el embudo de cribado. El método, sus supuestos y la lista corta están en [docs/quimica_candidatos.md](docs/quimica_candidatos.md).

## Datos y resultados

Los cuatro adjuntos se copiaron sin modificar a `data/raw/`. El GDSC2 completo se usa para auditoría; el primer entrenamiento utiliza las 51 líneas del subconjunto GBM/glioma. No se concatenan la versión CSV/XLSX ni T98G.

`config.json` preespecifica `LN_IC50`, RMSE máximo 0.3, semilla 42, holdout de 20% de líneas y cuatro folds de desarrollo. La partición se comparte entre los notebooks de entrenamiento. `reports/audit/manifest.json` conserva checksums de originales. Los archivos grandes, datos derivados, modelos y entorno están excluidos de Git; los notebooks, código, documentación y configuración sí pueden versionarse. Los datos originales deben transferirse por separado al reproducir el repositorio.

Resultados locales:

- `reports/audit/`: diccionarios CSV, inventario, retención por RMSE, cobertura y relaciones verificadas.
- `data/processed/`: tablas canónicas en Parquet y `line_split.csv`.
- `reports/metadata_baseline/`: búsqueda CV, predicciones, métricas por línea/fármaco, bootstrap y protocolo de ejecución.
- `reports/molecular/`: disponibilidad molecular o métricas reales cuando se incorpora expresión.
- `reports/chemistry/`: resolución de identidad química y filas que requieren revisión manual.
- `reports/candidates/`: ranking completo por `DRUG_ID`, lista corta por molécula y parámetros del cribado.
- `models/`: modelos ajustados sólo con desarrollo, más lista de features y configuración.

## Interpretación

La referencia de metadatos cuantifica cuánto puede predecirse sin medir biología molecular. Un R² global alto puede provenir de diferencias entre fármacos; revisar errores por fármaco y la mejora frente a su media. La diana `PUTATIVE_TARGET` no es la etiqueta `y`, y `PATHWAY_NAME` no es una medición celular de vías.

El notebook molecular implementa PCA, no puntuación de vías. Toda selección de genes, imputación, escala y PCA se aprende en entrenamiento dentro de CV. Los mínimos de cobertura son requisitos de ejecución, no un cálculo de tamaño muestral. Se verifican relaciones de modelos/donantes con la tabla local y se impiden cruces entre familias conocidas. La autenticidad de descarga de esos metadatos no está verificada.

Protocolo, fuentes, limitaciones y decisiones: `docs/metodologia.md`. Formato molecular y adquisición: `docs/expresion.md`. Primera ejecución y sus métricas: `docs/primer_resultado.md`. Los resultados son exploratorios in vitro; no validan predicciones clínicas para pacientes.


## Flujo corregido y conservación de ejecuciones

```bash
.venv/bin/python scripts/prepare_expression.py
.venv/bin/python scripts/run_notebooks.py
.venv/bin/python scripts/build_compounds.py
.venv/bin/python scripts/rank_candidates.py
.venv/bin/python scripts/update_results_docs.py
.venv/bin/python -m pytest -q
```

Los escenarios están en `screening.json`; las decisiones pendientes del laboratorio en `pilot_protocol.json`. Identidad/ruta revisadas y mediciones reales se registran en `data/candidate_reviews.csv` y `data/formulation_measurements.csv`. Las listas son provisionales; un campo vacío no se convierte en un supuesto favorable.

Cada etapa guarda un directorio distinto `reports/runs/<etapa>/<run_id>/` con snapshots y manifiesto. Las rutas cortas y modelos son enlaces a la corrida actual. Se preservan salidas anteriores; un fallo no vuelve a exponer automáticamente resultados viejos. Para reproducir en otro equipo hay que conservar también originales, ZIP y corridas, excluidos de Git.
