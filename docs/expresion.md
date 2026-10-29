# Expresión real: preparación y contrato

El ZIP `data/external/rnaseq_all_20220624.zip` está disponible localmente. `scripts/prepare_expression.py` convierte directamente su miembro `rnaseq_tpm_20220624.csv` al contrato del notebook molecular, sin extraer la tabla larga de 5,32 GB. Las cifras y el resultado del entrenamiento vigente están en [resultados_actuales.md](resultados_actuales.md).

## Qué se verifica

- Formato histórico: cabecera SIDM y filas `model_name`, `dataset_name`, `data_source`, `gene_id`.
- Un perfil único por SIDM; sin coincidencias aproximadas de nombres, duplicados o selección por outcome.
- Genes SIDG no repetidos, filas completas, valores no negativos y sin infinitos.
- Escala de TPM sin transformar: las sumas de mediciones por perfil deben estar dentro de 1e6 ±2%. El archivo local cumple este control. Esta comprobación de consistencia no certifica autenticidad ni procesamiento original.
- Aplicación de `log2(TPM+1)` exactamente una vez, transposición y selección de SIDM con respuesta. No se centran/escalan perfiles ni se imputan valores en esta preparación.
- Se lee el miembro TPM completo, verificando su CRC. Se registra además el SHA-256 del ZIP; no se afirma igualdad de bytes con una descarga oficial.

La documentación actual de Sanger describe una versión más reciente con TPM transformado; no se aplica automáticamente esa escala a la distribución histórica. [Fuente primaria](https://depmap.sanger.ac.uk/documentation/datasets/expression/).

## Resultado local y ausencias

Se emparejan **50 de 51** líneas y **37.602 genes**. `SIDM01021` no aparece en el encabezado de esta matriz; no se inventa un perfil. Cada línea emparejada contiene 339 genes ausentes en este archivo. Se conservan como NaN; el imputer del pipeline aprende sólo dentro del entrenamiento de cada fold. Los genes completamente vacíos o constantes en entrenamiento no superan la selección de varianza.

`reports/expression_preparation/` conserva `profiles.csv` con origen Sanger/Broad por línea, `gene_annotation.csv` con símbolos, `unmatched_lines.csv`, matriz, metadatos y manifiesto. Los features se llaman `gene_SIDG...`: son identificadores Sanger, no Ensembl simulados.

## Contrato de entrada

`data/external/expression.csv` tiene una fila por `SANGER_MODEL_ID` y columnas `gene_<identificador>`. `expression_metadata.json` declara `source_url`, `release`, `scale`, `gene_id_type`, `preprocessing` y `sha256`, además del origen y checksums registrados por la preparación. El loader comprueba el checksum, esquema, IDs únicos, ausencia de infinitos y escala admitida (`log2_tpm_plus1` o `rma_log2`). No aceptar counts crudos como si fueran TPM, ni repetir el logaritmo sobre datos ya transformados.

```bash
.venv/bin/python scripts/prepare_expression.py
.venv/bin/python scripts/run_notebooks.py 03_expresion_pca_modelos.ipynb
```

Las rutas canónicas son enlaces a la preparación actual. Si una nueva preparación falla, no se recurre silenciosamente a la matriz anterior. Cada preparación y entrenamiento conserva su directorio con run_id y estado. Si se incorpora otra distribución, su escala, perfiles y procedencia requieren un preparador documentado y pruebas correspondientes.

El modelo de expresión predice ajustes de fármaco libre en líneas reservadas. No predice carga, liberación, transporte ni respuesta de nanopartículas. El holdout ya observado se usa como evaluación exploratoria, y la extrapolación del outcome sigue siendo una limitación científica.
