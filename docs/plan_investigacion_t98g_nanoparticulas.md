# Plan de investigación: respuesta de T98G a fármacos y nanopartículas de óxido de hierro

Fecha del inventario: 7 de octubre de 2026.

Estado: propuesta metodológica; no representa un modelo de nanopartículas entrenado ni una eficacia experimental demostrada. El inventario describe los archivos locales. Los datos que pueda tener el equipo fuera de este repositorio deben incorporarse y verificarse.

## 1. Objetivo y alcance

Desarrollar y validar un sistema que utilice la estructura química del fármaco, su potencia previa y su interacción con nanopartículas de óxido de hierro para priorizar formulaciones capaces de reducir la viabilidad o provocar muerte celular en T98G, bajo condiciones experimentales definidas.

El resultado práctico será una lista de candidatos con respuesta esperada, incertidumbre y probabilidad de alcanzar un objetivo experimental. La validación será in vitro en T98G; extender las conclusiones a otras líneas, tumores o pacientes requerirá estudios adicionales.

### Preguntas que debemos separar

| Pregunta | Entradas principales | Resultado a predecir | Qué hace falta |
|---|---|---|---|
| A. ¿Qué moléculas tienen actividad como fármaco libre en T98G? | Estructura química; contexto celular si se usan otras líneas | LN_IC50 libre o curva de viabilidad | Respuestas al fármaco libre y estructuras verificadas |
| B. ¿Qué fármacos interactúan con una nanopartícula concreta? | Estructura; superficie de la nanopartícula; condiciones | Afinidad, carga o liberación, cada una por separado | Mediciones de interacción o simulaciones validadas |
| C. ¿Qué formulaciones tienen actividad en T98G? | Potencia libre, estructura, interacción, formulación, dosis y tiempo | Viabilidad o muerte celular de la formulación | Ensayos de fármaco cargado y controles |
| D. ¿La formulación mejora el efecto del fármaco libre? | Datos pareados a dosis equivalentes de fármaco | Diferencia de respuesta o de potencia | Comparación experimental libre/cargado |

La pregunta C es el objetivo principal propuesto. A y B aportan información y permiten avanzar antes de disponer de suficientes formulaciones. D es un objetivo secundario: una formulación puede ser activa sin superar al fármaco libre.

## 2. Qué significa «matar células»

LN_IC50 describe una concentración asociada a una reducción del 50 % en el ensayo correspondiente. En GDSC se estima a partir de curvas de viabilidad; una reducción del crecimiento o de la señal metabólica no demuestra por sí sola muerte celular. Véase la [documentación de GDSC/Sanger](https://depmap.sanger.ac.uk/documentation/datasets/drug-sensitivity/).

Proponemos mantener dos resultados separados:

1. **Respuesta de viabilidad:** fracción de viabilidad normalizada respecto al control, a dosis y tiempo definidos. Es el resultado continuo inicial.
2. **Muerte celular:** fracción de células muertas medida por un ensayo adecuado, por ejemplo Annexin V/PI con análisis y controles definidos. Es necesaria para sostener la afirmación de muerte celular.

Como ejemplo de evento para priorizar candidatos:

\[
P(V_{72h}<0.20\mid\text{estructura, potencia libre, interacción, formulación, dosis})
\]

El tiempo de 72 h y el umbral de 20 % son propuestas para discutir, no criterios aprobados. Esta probabilidad se refiere al resultado de un experimento bajo las condiciones especificadas, no a la probabilidad de muerte de una célula individual. Para muerte celular se definiría otro evento, basado en la medición de muerte.

La dosis debe estar fijada o formar parte de la entrada. No existe una etiqueta universal «mata/no mata» independiente de dosis, tiempo y formulación.

## 3. Datos que tenemos

### Inventario local

| Recurso | Contenido disponible | Uso y límites |
|---|---|---|
| `data/raw/GDSC2_fitted_dose_response_27Oct23.xlsx` | Según la auditoría existente: 242.036 curvas, 969 modelos celulares y 295 DRUG_ID | Respuesta a fármacos; no describe formulaciones con nanopartículas |
| `data/raw/GBMGliomafiltrado.csv` y `.xlsx` | 12.001 curvas, 51 modelos GBM/glioma; ambas versiones representan la misma tabla | Desarrollo en líneas relacionadas; no concatenar CSV y XLSX |
| `data/raw/Aislados T98G - Hoja 1.csv` | 280 curvas, 280 DRUG_ID y 271 nombres de fármacos; SIDM01171 | Respuestas de una sola línea; registros incluidos en GBM/glioma y GDSC2 |
| `data/processed/glioma_flagged.parquet` y `glioma_qc.parquet` | Tablas derivadas del pipeline actual | Reutilizables con trazabilidad al origen |
| `data/processed/line_split.csv` | Partición existente de 51 líneas: 40 en desarrollo y 11 en test | Protocolo de generalización a líneas reservadas; no evalúa química nueva |
| `data/external/model_list_20260921.csv` | 2.266 registros y 98 columnas de metadatos | Correspondencia de identificadores y relaciones; no es expresión génica |
| `data/external/rnaseq_all_20220624.zip` | Archivo descargado de unos 941 MB; contiene tablas de RNA-seq, TPM, FPKM y counts | Potencial fuente de expresión; falta preparar escala, orientación, unión y metadatos |
| Notebooks, modelos e informes existentes | Auditoría, referencias de metadatos y pipeline de expresión | Base reproducible; todavía no implementan el objetivo de nanopartículas |

Los conteos del GDSC2 completo y de líneas proceden de la [auditoría existente](primer_resultado.md). En esta revisión se comprobaron las dimensiones de los CSV, los identificadores de T98G, la existencia de archivos y el directorio interno del ZIP; no se volvió a auditar íntegramente el XLSX ni se procesó la matriz de RNA-seq.

### Variables presentes en las tablas de respuesta

Las tablas tienen 16 columnas: `DATASET`, `NLME_RESULT_ID`, `NLME_CURVE_ID`, `CELL_LINE_NAME`, `SANGER_MODEL_ID`, `CANCER_TYPE`, `DRUG_ID`, `DRUG_NAME`, `PUTATIVE_TARGET`, `PATHWAY_NAME`, `MIN_CONC`, `MAX_CONC`, `LN_IC50`, `AUC`, `RMSE` y `Z_SCORE`.

- `LN_IC50`: resultado ajustado de sensibilidad al fármaco. No aplicar otro logaritmo; en este proyecto la concentración de referencia está en µM.
- `AUC`: resumen de la curva en el rango ensayado. No mezclar directamente AUC de rangos o protocolos distintos.
- `MIN_CONC`/`MAX_CONC`: límites del ensayo; no son una serie de dosis con sus respuestas crudas.
- `RMSE`: calidad del ajuste experimental; no es error del modelo predictivo.
- `PUTATIVE_TARGET`/`PATHWAY_NAME`: anotaciones del fármaco; no son mediciones de actividad biológica de T98G.

### Datos que no están preparados o no se encontraron

| Dato | Estado local |
|---|---|
| SMILES e identificadores químicos verificados | No encontrados en las tablas disponibles |
| Afinidad fármaco–nanopartícula | No encontrada |
| Caracterización de nanopartículas y formulaciones | No encontrada |
| Carga efectiva y liberación del fármaco | No encontradas |
| Viabilidad o muerte de T98G expuesta a formulaciones | No encontradas |
| Curvas crudas por dosis y réplica | No encontradas en los adjuntos de respuesta ajustada |
| Expresión en formato admitido por el pipeline | No existen `expression.csv` ni `expression_metadata.json`; sí existe el ZIP de origen |

La documentación anterior describe expresión pendiente. El inventario actual precisa que el archivo de origen ya está descargado, pero todavía no equivale a una entrada molecular validada. El ZIP contiene una tabla de aproximadamente 5,32 GB sin comprimir y otras menores; conviene preparar la tabla adecuada sin cargar todo indiscriminadamente.

### Límites científicos ya detectados

- Las 280 curvas de T98G no son 280 pacientes, 280 formulaciones ni 280 réplicas biológicas.
- Los 271 nombres tampoco garantizan 271 estructuras distintas. Hay que resolver sales, alias y múltiples IDs antes de contar moléculas.
- En T98G, 214 de 280 IC50 ajustadas están por encima de la concentración máxima ensayada: **76,43 %**. En el conjunto GBM/glioma, la auditoría documenta **76,50 %**.
- Una IC50 extrapolada no debe tratarse como una medición precisa. Si se dispone de datos crudos, se revisará el ajuste; cuando el ensayo solo sustente un límite, se considerará una observación censurada. No convertir automáticamente toda extrapolación publicada en censura sin revisar su evidencia.
- La referencia actual de media por fármaco obtuvo MAE test 0,9389 y R² agrupado 0,7846. Ese R² no demuestra capacidad de predicción molecular ni eficacia de nanopartículas. [Resultados del proyecto](primer_resultado.md).

## 4. Datos que debemos buscar o generar

| Prioridad | Datos | Fuente o acción propuesta | Verificación necesaria |
|---|---|---|---|
| P0 | Identidad química y SMILES | Anotaciones oficiales de compuestos GDSC; contrastar con PubChem/ChEMBL o proveedor | Unir por identificadores; revisar isómeros, sales, mezclas, alias y discrepancias |
| P0 | Definición de nanopartícula objetivo | Equipo experimental | Óxido, recubrimiento, superficie, tamaño y condiciones de uso |
| P0 | Definición operacional de afinidad | Equipo experimental y literatura primaria de la misma superficie | Magnitud, método, unidades, pH, temperatura y medio |
| P0 | Respuesta de T98G a libre/cargado | Generar un piloto experimental | Dosis equivalentes de fármaco, controles, réplicas y lotes |
| P1 | Carga y liberación | Mediciones de las formulaciones | Fracción realmente cargada y curva de liberación en medios pertinentes |
| P1 | Curvas crudas de respuesta | Buscar datos originales de GDSC y generar los del piloto | Dosis, tiempo, réplica, controles y trazabilidad del ensayo |
| P1 | Expresión basal | Preparar el ZIP disponible y cruzar con metadatos | Escala, versión, genes, SIDM, cobertura y origen; contrato en [expresion.md](expresion.md) |
| P1 | Datos publicados de formulaciones comparables | Literatura primaria y suplementos | Separar estudios por línea, superficie, ensayo y condiciones; no asumir transferibilidad |
| P2 | Internalización celular y mecanismo de muerte | Ensayos adicionales | Mediciones independientes y comparables |
| P2 | Selectividad | Ensayos con un modelo no tumoral pertinente | Comparaciones en condiciones equivalentes |

La expresión aporta información si se aprende entre distintas líneas celulares. Si todas las observaciones corresponden a T98G con el mismo perfil basal, ese perfil es constante y no explica diferencias entre fármacos por sí solo.

### Afinidad, carga y liberación son variables distintas

Debemos elegir qué significa «afinidad». Puede ser una constante de adsorción, un Kd o una energía de interacción calculada, pero esas magnitudes no son intercambiables. El porcentaje de carga es una propiedad de formulación dependiente de condiciones; no es automáticamente una constante de afinidad.

Registrar por separado:

- Medición o cálculo de interacción, con unidades y condiciones.
- Cantidad de fármaco cargado por masa de nanopartícula y eficiencia de carga.
- Fracción liberada a lo largo del tiempo, con pH, medio y temperatura.

Una interacción fuerte podría favorecer la carga y dificultar la liberación. No imponer una relación monotónica «más afinidad → más muerte». Si se usan cálculos de interacción, validar su relación con mediciones experimentales para la superficie concreta; un cálculo aislado no constituye una etiqueta de citotoxicidad.

## 5. Inputs y contratos de datos propuestos

Esta es una propuesta de tablas futuras, no archivos ya implementados. Conservar datos originales y procedencia; no rellenar ausencias con números ficticios.

| Tabla | Campos mínimos propuestos | Unidad de registro |
|---|---|---|
| `compounds.csv` | `compound_id`, `gdsc_drug_id`, `drug_name`, `canonical_smiles`, `inchikey`, `chemical_form`, `source`, `release` | Correspondencia de ID de origen con identidad química validada |
| `nanoparticles.csv` | `np_id`, `np_batch_id`, `oxide_type`, `coating`, `core_size_nm`, `hydrodynamic_size_nm`, `pdi`, `zeta_mv`, condiciones y métodos | Material y lote caracterizado |
| `interaction_measurements.csv` | `measurement_id`, `compound_id`, `np_id`, `np_batch_id`, `metric_type`, `value`, `unit`, `method`, `ph`, `temperature_c`, `medium`, `replicate_id`, `source` | Una medición o cálculo de interacción |
| `formulations.csv` | `formulation_id`, `compound_id`, `np_id`, `np_batch_id`, `formulation_batch_id`, `loading_mass_ratio`, `loading_efficiency`, `preparation_protocol` | Formulación y lote |
| `release_measurements.csv` | `formulation_id`, `time_h`, `released_fraction`, `ph`, `temperature_c`, `medium`, `replicate_id` | Medición de liberación |
| `responses.csv` | `experiment_id`, `cell_model_id`, `compound_id`, `formulation_id`, `treatment_type`, `drug_conc_um`, `np_conc_ug_ml`, `time_h`, `assay`, `plate_id`, `biological_replicate_id`, `technical_replicate_id`, `raw_readout`, `normalized_response`, `response_type`, `qc_status` | Una observación experimental por pozo o muestra |
| `free_drug_potency.csv` | `compound_id`, `cell_model_id`, `ln_ic50_free`, rango de dosis, método, tiempo, fuente, incertidumbre o estado de censura | Potencia libre en un contexto definido |

Se necesita además un registro de controles y del procedimiento de normalización. Las concentraciones propuestas son unidades canónicas: documentar conversiones y definir si la dosis de nanopartícula representa masa total o masa de hierro. Para nanopartícula vacía o vehículo, permitir claves nulas según el tipo de tratamiento y reglas explícitas.

Las liberaciones y caracterizaciones deben corresponder al lote y protocolo usados. No incorporar mediciones que solo se conocerán después del resultado que se pretende predecir.

### Entradas del predictor final

1. Estructura: fingerprint Morgan/ECFP y descriptores químicos, o representación aprendida del grafo molecular.
2. Potencia libre: LN_IC50 previa en T98G, con contexto e incertidumbre, cuando esté disponible.
3. Interacción: medida definida de afinidad; carga y liberación disponibles antes del ensayo de respuesta.
4. Nanopartícula/formulación: propiedades de superficie y tamaño, entre otras caracterizaciones pertinentes.
5. Exposición: dosis de fármaco, dosis de nanopartícula, tiempo y condiciones relevantes.

Debe declararse el escenario de uso:

- **Candidato caracterizado:** tiene potencia libre y caracterización experimental. Podemos usarlas para predecir la respuesta de la formulación.
- **Molécula nueva sin potencia medida:** necesitamos otro predictor para estimar potencia, o un modelo que opere sin ella. Evaluar ese escenario por separado y propagar la incertidumbre; no asumir rendimiento equivalente.

Si la nanopartícula y el tiempo se mantienen fijos, pueden ser condiciones del protocolo en lugar de columnas variables. El modelo solo estará validado para ese contexto.

## 6. Outputs esperados

| Salida | Significado |
|---|---|
| Viabilidad esperada | Respuesta continua para la combinación de compuesto, formulación, dosis y tiempo |
| Intervalo predictivo | Rango de posibles resultados experimentales; especificar nivel y método |
| Probabilidad de alcanzar el umbral | Probabilidad predictiva del evento previamente definido |
| Curva dosis–respuesta | Respuesta esperada en el rango validado, cuando los datos permitan ajustarla |
| Mejora frente al fármaco libre | Diferencia pareada de respuesta, con incertidumbre, a dosis equivalentes |
| Dominio de aplicación | Indicación de semejanza química y soporte de condiciones de formulación |
| Ranking para experimentos | Priorización que considere respuesta, incertidumbre y factibilidad |

Para el objetivo secundario puede definirse `delta_viability = viability_loaded - viability_free`: valores negativos indican menor viabilidad con la formulación. No etiquetar esa diferencia como sinergia sin un diseño y modelo de referencia específicos.

Archivo de predicciones propuesto:

```text
compound_id, formulation_id, drug_conc_um, np_conc_ug_ml, time_h,
predicted_viability, prediction_lower, prediction_upper,
probability_below_threshold, threshold, interval_level,
domain_status, model_version
```

La salida podrá indicar «sin soporte suficiente» cuando falten entradas esenciales o la combinación esté fuera del dominio validado. Los intervalos predictivos de un experimento y los intervalos de la media esperada deben distinguirse.

## 7. Aproximación probabilística

### Modelo inicial recomendado

Comparar un modelo bayesiano jerárquico de dosis–respuesta con referencias sencillas. Modelar la respuesta según concentración y tiempo, usando regularización para las características químicas y efectos de lote/experimento cuando haya datos suficientes para estimarlos.

Una estructura conceptual es:

```text
Química + potencia libre + interacción + formulación
                    ↓
       Parámetros de respuesta a la dosis
                    ↓
       Dosis y tiempo de exposición
                    ↓
Distribución predictiva de viabilidad o muerte
                    ↓
  Probabilidad de alcanzar el objetivo experimental
```

La forma de la curva y la distribución de observación se elegirán según el ensayo. No asumir que toda formulación admite una sigmoide simple. Si la respuesta está realmente acotada entre 0 y 1, considerar distribuciones compatibles; si la normalización produce valores fuera del rango, revisar el proceso y modelar el readout pertinente en lugar de truncar automáticamente.

Un proceso gaussiano sobre una representación química manejable es otra alternativa a comparar con pocos datos independientes. La jerarquía y la incertidumbre no eliminan la necesidad de datos representativos ni garantizan calibración fuera del dominio observado.

La probabilidad final se obtendrá de la distribución predictiva del resultado futuro, incorporando ruido experimental y la incertidumbre de las entradas estimadas. La cobertura de intervalos debe verificarse en datos reservados.

## 8. Machine learning y deep learning

| Nivel | Modelos propuestos | Propósito |
|---|---|---|
| Referencias | Media de desarrollo; potencia libre como predictor; modelo simple de dosis y tiempo | Determinar cuánto aporta la complejidad |
| ML regularizado | Elastic Net; regresión logística si la etiqueta binaria es necesaria | Referencia interpretable con descriptores |
| ML no lineal | Random Forest/Extra Trees y gradient boosting | Capturar interacciones entre estructura y variables numéricas |
| Probabilístico | Modelo bayesiano jerárquico; proceso gaussiano | Distribución de respuesta, incertidumbre y probabilidades |
| Deep learning | D-MPNN u otra red de grafos, fusionada con entradas experimentales | Aprender una representación molecular para el objetivo definido |

Para deep learning, convertir SMILES en un grafo de átomos y enlaces; combinar el embedding molecular con potencia libre, interacción, formulación y exposición mediante una cabeza de predicción. Evaluar preentrenamiento o transferencia solo si aportan mejora en la validación correspondiente. Véase [Yang y colaboradores](https://dspace.mit.edu/entities/publication/c07288ef-8bfa-4557-abb7-b387e1d1e67e).

No existe un número universal de filas que justifique una red neuronal. Contar compuestos, familias químicas, formulaciones, lotes y experimentos independientes; evaluar curvas de aprendizaje. MoleculeNet documenta que modelos clásicos pueden ser más robustos bajo escasez de datos: [estudio original](https://doi.org/10.1039/C7SC02664A).

Una red entrenada solo en GDSC libre puede apoyar la tarea A; no queda validada para la tarea C por añadir una columna de afinidad. Esta última necesita etiquetas de formulaciones.

## 9. Validación y prevención de fuga de información

### Partición según la pregunta

| Escenario | Diseño de evaluación |
|---|---|
| Moléculas nuevas en T98G | Separación por compuesto y familia estructural; todos los IDs, dosis y réplicas de una identidad química viajan juntos |
| Nuevas formulaciones de compuestos conocidos | Reservar formulaciones; controlar qué lotes y compuestos se comparten según el uso previsto |
| Nuevo lote o experimento | Reservar lotes o campañas experimentales completas |
| Nueva línea celular | Reservar líneas completas y, si corresponde, donantes; es una pregunta distinta |
| Confirmación prospectiva | Predecir antes de generar respuestas de nuevos candidatos y congelar modelo y criterios |

La separación estructural permite evaluar generalización química de manera más exigente que una separación aleatoria de filas. [MoleculeNet](https://doi.org/10.1039/C7SC02664A).

### Reglas

- No usar LN_IC50 de T98G como entrada para predecir esa misma LN_IC50 o una etiqueta definida directamente a partir de ella.
- No usar AUC, Z_SCORE u otras salidas de la misma curva como predictores de esa respuesta.
- Para la tarea de formulación, la potencia libre puede ser entrada si está disponible en el uso previsto y proviene de un ensayo independiente.
- Ajustar imputación, selección, escalado, PCA, hiperparámetros y calibración exclusivamente dentro del desarrollo.
- No separar réplicas técnicas entre entrenamiento y test. Registrar dependencia entre réplicas biológicas, placas y lotes.
- Si una entrada proviene de otro modelo, generar valores fuera de muestra para entrenar la etapa siguiente; evitar introducir predicciones ajustadas sobre sus propias etiquetas.
- Mantener GDSC2, su subconjunto y T98G como vistas de una fuente, no como tres cohortes independientes.
- El test actual ya fue examinado. Reutilizarlo tras adaptar la estrategia a sus resultados es evaluación exploratoria; una nueva semilla no lo vuelve independiente.

### Métricas y comparación

- Respuesta continua: MAE, RMSE predictivo y errores por compuesto/formulación y rango de dosis.
- Incertidumbre: cobertura y anchura de intervalos; comprobar que la incertidumbre sea informativa.
- Probabilidad: Brier score, log loss y curva de calibración; PR-AUC cuando el evento positivo sea escaso.
- Priorización: precisión entre los primeros K candidatos, sensibilidad y enriquecimiento respecto a una referencia, con K fijado según capacidad experimental.
- Comparación de modelos: mismas particiones y evaluación de diferencias con incertidumbre agrupada por unidades independientes.
- Ablaciones: dosis/tiempo; añadir potencia libre; añadir estructura; añadir afinidad; añadir caracterización y liberación. Medir qué aporta cada bloque.

Si se estima la probabilidad a partir de una distribución continua, conservar la etiqueta continua y comprobar también la calibración del evento. El presupuesto experimental, K y los criterios de aceptación se fijarán antes de evaluar el test final.

## 10. Diseño del piloto experimental

Seleccionar candidatos que cubran diversidad química y combinaciones de potencia e interacción. Evitar ensayar únicamente los candidatos de mayor puntuación: incluir controles y candidatos informativos que permitan comprobar supuestos.

Comparar, bajo condiciones documentadas:

1. Control sin tratamiento y vehículo.
2. Fármaco libre.
3. Nanopartícula vacía.
4. Fármaco cargado en nanopartícula, con cuantificación de carga y dosis equivalente de fármaco.
5. Cuando sea pertinente, mezcla física libre + nanopartícula para distinguirla de la formulación cargada.

Usar una serie de dosis y tiempos preespecificados, aleatorización de posiciones, controles de placa y réplicas biológicas en experimentos independientes. El número de candidatos y réplicas se definirá a partir de variabilidad del piloto, precisión requerida y presupuesto; las réplicas técnicas no sustituyen a experimentos independientes.

Validar el ensayo en presencia de las nanopartículas, con controles sin células y una medición complementaria adecuada. Se ha demostrado interferencia de óxido de hierro con ensayos de viabilidad: [estudio experimental](https://pmc.ncbi.nlm.nih.gov/articles/PMC10845919/). Annexin V/PI también requiere sus controles y compensación correspondientes.

## 11. Plan de trabajo y entregables

| Etapa | Acciones | Entregable | Criterio de avance |
|---|---|---|---|
| 0. Definir el experimento | Acordar respuesta principal, superficie, afinidad, dosis, tiempos y escenario de predicción | Protocolo preespecificado y diccionario | Inputs medibles y etiqueta inequívoca |
| 1. Preparar química | Resolver IDs, obtener y validar SMILES, estandarizar y calcular descriptores | Tabla de compuestos, cobertura y conflictos | Identidades verificadas; exclusiones documentadas |
| 2. Preparar respuesta libre | Revisar extrapolaciones; buscar curvas crudas; conservar procedencia | Dataset de potencia con incertidumbre y QC | Definición explícita de valores fiables, límites y ausencias |
| 3. Preparar expresión si aplica | Procesar ZIP, cruzar SIDM y verificar escala | Matriz y metadatos del contrato existente | Unión y escala validadas; sin perfiles inventados |
| 4. Construir datos de formulación | Caracterizar interacción, carga y liberación; ejecutar piloto | Tablas experimentales completas | Controles válidos y respuestas reproducibles |
| 5. Entrenar referencias y probabilístico | Comparar modelos clásicos y bayesianos; medir aportes de entradas | Modelos, predicciones fuera de muestra e informe | Evaluación bajo el escenario de uso; incertidumbre comprobada |
| 6. Evaluar deep learning | Comparar GNN con el mejor modelo clásico | Experimento reproducible y curvas de aprendizaje | Mejora verificable; si no mejora, conservar el modelo más útil |
| 7. Validar prospectivamente | Congelar modelo y criterios; ensayar nuevos candidatos | Informe de predicciones vs. resultados | Cumplimiento de criterios predefinidos y dominio declarado |

Las etapas de preparación y búsqueda pueden avanzar en paralelo como actividades del equipo. La etapa 5 para formulaciones depende de la etapa 4; disponer de SMILES y afinidad sin respuestas de formulación no habilita un predictor validado de citotoxicidad cargada.

Este documento propone una ampliación del proyecto. No modifica el target, las exclusiones de variables ni la partición del [protocolo actual](metodologia.md). La inclusión de dosis y tiempo aquí corresponde a una nueva tarea de respuesta condicionada a exposición; requiere configuración y protocolo propios.

## 12. Decisiones pendientes y siguiente paso concreto

- ¿El objetivo principal será viabilidad, muerte celular o ambos con uno declarado primario?
- ¿Se predice una molécula nueva, una formulación nueva o un nuevo lote de una formulación conocida?
- ¿Qué óxido, recubrimiento y condiciones tendrán las nanopartículas?
- ¿La afinidad está medida, calculada o aún debe generarse? ¿Con qué magnitud y unidades?
- ¿Qué dosis, tiempos y ensayo se usarán? ¿Hay datos fuera del repositorio?
- ¿Se exige mejorar frente al fármaco libre o basta alcanzar una actividad mínima?
- ¿Cuántos candidatos puede confirmar el laboratorio y con qué precisión se requiere el ranking?

**Siguiente paso recomendado:** completar la definición experimental y construir la tabla de identidad química de los fármacos de T98G. Con ello se puede desarrollar un primer predictor de respuesta libre y seleccionar un piloto diverso. El predictor de formulaciones debe entrenarse y validarse después de disponer de respuestas reales de fármacos cargados.
