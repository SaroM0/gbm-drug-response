# Entrada molecular para el notebook 03

Los cuatro adjuntos no incluyen expresión. El notebook informa `missing_expression` hasta incorporar dos archivos reales en `data/external/`. No genera genes sintéticos ni interpreta las anotaciones de fármacos como mediciones moleculares.

## Contrato

`expression.csv`: filas = modelos celulares únicos; primera columna `SANGER_MODEL_ID`; restantes columnas `gene_<identificador>` con números. Ejemplo de cabecera, sin datos ficticios:

```csv
SANGER_MODEL_ID,gene_ENSG00000141510,gene_ENSG00000146648
```

Escala admitida: `log2_tpm_plus1` o `rma_log2`. No pasar counts crudos ni datos centrados/escalados con todas las líneas. TPM sin log requiere `log2(TPM+1)` una sola vez. Para matrices ya transformadas no repetirlo. Los genes deben usar el mismo sistema de identificadores y no repetirse. Se permiten NaN parciales; la imputación se aprende sólo en entrenamiento. Los IDs COSMIC/ACH necesitan una tabla de correspondencia oficial; no usar coincidencia aproximada de nombres.

`expression_metadata.json` debe incluir:

```json
{
  "source_url": "URL real de descarga",
  "release": "versión real",
  "scale": "log2_tpm_plus1",
  "gene_id_type": "Ensembl gene ID sin versión, u otra definición real",
  "preprocessing": "Describir orientación, selección de perfil por SIDM y transformaciones realizadas",
  "sha256": "SHA-256 del expression.csv final"
}
```

El hash puede obtenerse con `sha256sum data/external/expression.csv`. No se incluyen archivos vacíos con estos nombres: el notebook debe detectar la ausencia, no confundir una plantilla con datos.

## Fuentes y preparación

La documentación de Sanger enlaza las descargas oficiales de Cell Model Passports:
https://depmap.sanger.ac.uk/documentation/datasets/expression/
https://cellmodelpassports.sanger.ac.uk/downloads

La matriz ancha oficial puede estar orientada genes × muestras. Identificar filas de metadatos y transponer únicamente las mediciones. Verificar la escala contra el README de la versión descargada; la documentación actual reporta TPM transformado. Mantener trazabilidad entre SIDM, perfil y plataforma; si hay varios perfiles de un SIDM, preespecificar cuál conservar y documentarlo.

También se verificó que la descarga histórica `https://cog.sanger.ac.uk/cmp/download/rnaseq_all_20220624.zip` existe (940.707.972 bytes). No se descargó automáticamente: falta elegir y documentar versión, perfil y escala. El notebook admite la matriz canónica, independientemente de la distribución elegida. La antigua URL RMA de cancerrxgene.org devolvió HTTP 410 en esta sesión.

Tras incorporar los archivos, ejecutar notebook 03 completo. Se registra cobertura de unión y SIDM sin expresión; no se imputa un perfil molecular entero para las líneas no emparejadas. La selección de hasta tres fármacos se basa sólo en cobertura de desarrollo. Para entrenar más, cambiar `molecular_max_drugs` o `molecular_drug_ids` antes de examinar test y volver a ejecutar como experimento exploratorio.
