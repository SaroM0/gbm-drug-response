# gbm-drug-response

Primer acercamiento reproducible para auditar GDSC2 y predecir `LN_IC50` en líneas celulares GBM/glioma reservadas. CPU; sin deep learning ni bases de datos adicionales.

## Documentación

Empezar por el [informe completo del proyecto y explicación de resultados](docs/informe_completo.md). Describe los datos, el diccionario, el código, las decisiones, las pruebas, la reproducción y los límites científicos de la primera ejecución.

| Documento | Contenido |
|---|---|
| [Informe completo](docs/informe_completo.md) | Trabajo realizado y lectura detallada de MAE, RMSE, R², Spearman, extrapolaciones y variación entre líneas. |
| [Primer resultado](docs/primer_resultado.md) | Resumen de las métricas y conclusión de la primera ejecución. |
| [Metodología](docs/metodologia.md) | Protocolo preespecificado, validación, QC y fuentes. |
| [Entrada de expresión](docs/expresion.md) | Archivos y metadatos necesarios para ejecutar el entrenamiento molecular real. |
| [Evidencia de la primera ejecución](docs/evidencia_primera_ejecucion.json) | Instantánea versionable de hashes, configuración, partición, QC y métricas. |

Resultado principal: la media por fármaco fue seleccionada en CV y obtuvo MAE 0,9389 en test; los modelos con metadatos no la superaron. El R² agrupado de 0,7846 no acredita predicción individual molecular. El 76,50% de las IC50 ajustadas de GBM/glioma está por encima del máximo ensayado, una limitación relevante del outcome. El modelo molecular sigue pendiente de expresión real.

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

## Datos y resultados

Los cuatro adjuntos se copiaron sin modificar a `data/raw/`. El GDSC2 completo se usa para auditoría; el primer entrenamiento utiliza las 51 líneas del subconjunto GBM/glioma. No se concatenan la versión CSV/XLSX ni T98G.

`config.json` preespecifica `LN_IC50`, RMSE máximo 0.3, semilla 42, holdout de 20% de líneas y cuatro folds de desarrollo. La partición se comparte entre los notebooks de entrenamiento. `reports/audit/manifest.json` conserva checksums de originales. Los archivos grandes, datos derivados, modelos y entorno están excluidos de Git; los notebooks, código, documentación y configuración sí pueden versionarse. Los datos originales deben transferirse por separado al reproducir el repositorio.

Resultados locales:

- `reports/audit/`: diccionarios CSV, inventario, retención por RMSE, cobertura y relaciones verificadas.
- `data/processed/`: tablas canónicas en Parquet y `line_split.csv`.
- `reports/metadata_baseline/`: búsqueda CV, predicciones, métricas por línea/fármaco, bootstrap y protocolo de ejecución.
- `reports/molecular/`: disponibilidad molecular o métricas reales cuando se incorpora expresión.
- `models/`: modelos ajustados sólo con desarrollo, más lista de features y configuración.

## Interpretación

La referencia de metadatos cuantifica cuánto puede predecirse sin medir biología molecular. Un R² global alto puede provenir de diferencias entre fármacos; revisar errores por fármaco y la mejora frente a su media. La diana `PUTATIVE_TARGET` no es la etiqueta `y`, y `PATHWAY_NAME` no es una medición celular de vías.

El notebook molecular implementa PCA, no puntuación de vías. Toda selección de genes, imputación, escala y PCA se aprende en entrenamiento dentro de CV. Los mínimos de cobertura son requisitos de ejecución, no un cálculo de tamaño muestral. Las líneas se agrupan por SIDM; aún no se ha demostrado independencia por donante.

Protocolo, fuentes, limitaciones y decisiones: `docs/metodologia.md`. Formato molecular y adquisición: `docs/expresion.md`. Primera ejecución y sus métricas: `docs/primer_resultado.md`. Los resultados son exploratorios in vitro; no validan predicciones clínicas para pacientes.
