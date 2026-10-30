# Piloto preparado para revisión del laboratorio

Estado: propuesta de caracterización y estimación de variabilidad, **sin aprobación ni resultados experimentales**. Se preparó a partir de las condiciones recibidas (T98G, Fe₃O₄@SiO₂-APTES, MTT, 5–25 ppm, 24/48/72 h y unos tres candidatos). Los supuestos añadidos están en `pilot_design.json`; no sustituyen decisiones del laboratorio ni cambian el estado pendiente de `pilot_protocol.json`.

## 1. Qué está listo

`reports/pilot_preparation/plate_layout.csv` contiene posiciones aleatorizadas y pares de placas. `candidate_review_packet.csv` reúne identidad y ensayos de los tres candidatos; `sample_size_scenarios.csv` cuantifica escenarios de potencia. Cada salida pertenece a una corrida con código, configuración y checksums.

Propuesta de recursos: placas de 96 pozos; tres réplicas técnicas; tres campañas biológicas independientes para el piloto de variabilidad. Por candidato, campaña y tiempo se propone una placa celular y otra sin células para estudiar interferencia. Las mediciones MTT a distintos tiempos se hacen en placas separadas: no se reutiliza un mismo pozo como si el ensayo fuera longitudinal no destructivo.

Esto suma **54 placas** para los tres candidatos, tres campañas y tres tiempos, incluidas las placas sin células. Cada placa tiene 72 posiciones planificadas para lectura y 24 posiciones de medio sin células fuera del análisis. Es una estimación de recursos para discutir, no una orden de laboratorio. Si excede la capacidad disponible, ejecutar primero una etapa de factibilidad para una formulación/tiempo y revisar el presupuesto antes del plan completo; reducir campañas no puede presentarse como un estudio confirmatorio con potencia demostrada.

## 2. Brazos y controles

Cada placa celular propone cinco niveles de dosis para cuatro grupos: libre, NP vacía, NP cargada y mezcla física libre+NP, con tres pozos técnicos por combinación. Son 60 pozos; se añaden seis de vehículo, tres de control positivo del ensayo y tres blancos sin células. El brazo de mezcla física es adicional a los cuatro brazos acordados originalmente y permanece propuesto.

El reactivo del control positivo, condiciones de exposición, lector, vehículo y procedimiento MTT requieren selección/validación del laboratorio. Las placas sin células reproducen posiciones y exposiciones para detectar señal de NP/fármaco/reagentes. Un fondo medible no demuestra que restarlo sea suficiente: la validez de corrección y la linealidad deben comprobarse. Se propone recuento celular por imagen como ensayo ortogonal si está disponible; no se etiqueta una reducción de MTT como muerte celular sin otra evidencia.

La aleatorización tiene semilla fija, conserva los controles y no se vuelve a sortear según los resultados. Los IDs `PLAN-*` son posiciones propuestas, no experimentos realizados. Todas las lecturas, cargas y dosis numéricas de fármaco quedan vacías hasta medirlas.

## 3. Dosis equivalentes y elección de rango

La NP cargada y la mezcla se comparan con el libre a la misma masa nominal de entidad activa incorporada, calculada a partir de la carga real del lote. Se registra por separado cuánto se libera y cuándo: no se asume que una fracción liberada tardíamente equivalga a la misma exposición que fármaco libre añadido al inicio.

La masa de fármaco, del portador total y de hierro no son intercambiables. Para una carga medida `L%` (masa de entidad activa/masa de portador), PM activo y dosis NP en µg/mL:

```text
dosis nominal de fármaco (µM) = dosis NP × L / 100 × 1000 / PM
```

Los cinco niveles 5/10/15/20/25 µg/mL preservan la propuesta recibida. No garantizan un rango informativo para todos los candidatos: docetaxel puede quedar muy por encima de su referencia libre y OSI-027 por debajo, según carga/liberación. El piloto debe determinar si hay respuesta a ambos lados de 50%; una ampliación de rango se documenta antes del estudio comparativo. Si no hay cruce, se informa IC50 no identificable o límite respaldado por la curva, sin fabricar una estimación puntual.

## 4. Resultado principal y análisis propuesto

Se propone **72 h** como tiempo principal para discutir, por la referencia GDSC2 y el antecedente local de oxaliplatino a ese tiempo. Los otros tiempos pueden caracterizar cinética. Esto no certifica comparabilidad de MTT y CellTiter-Glo ni constituye una decisión ya aprobada.

La pregunta farmacológica sigue siendo si la formulación mejora frente al libre a exposición nominal equivalente. Antes del ensayo comparativo deben fijarse el estimando principal, el efecto mínimo relevante y el intervalo de dosis común. Si se elige la razón de IC50, el contraste sería `ln(IC50_cargado / IC50_libre)`, estimado por campaña pareada. Si una o ambas curvas no identifican IC50, no excluirlas silenciosamente ni declarar superioridad a partir de los pares restantes; revisar el rango o usar un método de intervalos/censura justificado y predefinido.

Las réplicas técnicas se resumen dentro de cada campaña. La unidad de inferencia es la campaña independiente, con lotes de formulación y placa registrados. Curvas y supuestos se revisan por campaña: una sigmoide forzada sobre respuesta no monótona no resuelve el problema. Mantener lecturas crudas, controles y normalización; evitar truncar automáticamente respuestas fuera de 0–1.

Para tres comparaciones principales se propone control familiar mediante Holm. La tabla de planificación utiliza Bonferroni bilateral, una aproximación conservadora para presupuestar muestras; no escoge el candidato o tiempo con el mejor p-valor. Intervalos y tamaño del efecto acompañan cualquier contraste.

## 5. Tamaño muestral: escenarios, no una varianza inventada

La tabla calcula potencia con t no central para diferencias pareadas de log-IC50, suponiendo campañas independientes, diferencias normales y ambas IC50 identificables. Considera reducciones de IC50 del 20%, 30% o 50% y SD de log-razón entre 0,2 y 1. Ninguna de esas SD procede de un experimento del proyecto.

Por ejemplo, **si** se eligiera una reducción del 30% y **si** la SD pareada fuese 0,4, se necesitarían 17 campañas para 80% de potencia bajo esa planificación conservadora. Tres campañas tendrían aproximadamente 5,5% de potencia en ese escenario. Es un ejemplo de sensibilidad, no la recomendación definitiva para el estudio. La varianza real debe venir del piloto y la reducción relevante de una decisión científica previa. [Distribución implementada por SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.nct.html).

## 6. Expediente de los candidatos

| Candidato | Evidencia local | Revisión necesaria |
|---|---|---|
| Docetaxel | IDs 1007 y 1819; IC50 libre 5,15–8,23 nM; PubChem CID 148124 | Sitio de modificación/linker, carga y actividad/liberación; rango inicial puede saturarse |
| Gemcitabina | ID 1190; IC50 libre 349 nM; PubChem CID 60750 | Forma comercial y masa activa; comparar rutas de retención/liberación en vez de asumir una química única |
| OSI-027 | ID 1594; IC50 libre 3,57 µM; PubChem CID 135398516, título Cerc 006 | Confirmar alias/forma adquirida; conservación de actividad del derivado; carga y liberación suficientes |

La anotación oficial de GDSC 8.5 confirma esos DRUG_ID y nombres; no aporta por sí misma la identidad del lote comercial. Ninguna revisión química se declara aprobada automáticamente.

Para docetaxel hay evidencia primaria de conjugados con un polímero carboxilado, pero en otro material y otras líneas: apoya revisar la ruta, no copiar su rendimiento a Fe₃O₄@SiO₂-APTES. [Estudio de conjugación HA-docetaxel](https://d-nb.info/110939456X/34).

Para gemcitabina hay un ejemplo de carga no covalente en sílice mesoporosa aminada; difiere de nuestro portador y su capacidad de carga no se extrapola. Justifica considerar retención no covalente como hipótesis adicional, aunque el cribado estructural proponga derivatización. [Estudio original](https://pubmed.ncbi.nlm.nih.gov/27016140/).

Para OSI-027, la caracterización preclínica confirma su farmacología como inhibidor mTORC1/2, pero no establece que conjugar el carboxilo a APTES conserve actividad o libere el compuesto. [Caracterización original](https://pubmed.ncbi.nlm.nih.gov/21673091/).

## 7. Endotelio y Transwell

Antes de interpretar transporte se deben definir analito y método: fármaco libre, total, NP o marcador no son el mismo resultado. Se propone cuantificar compartimentos donante/receptor, fracción asociada a células/membrana cuando sea medible y pérdidas, con tiempos, área y unidades explícitos. El cálculo de permeabilidad aparente sólo tiene sentido si se justifica el flujo, la condición de sumidero y la concentración de referencia correspondiente al analito.

La integridad de barrera debe comprobarse antes y después mediante el método validado por el laboratorio, complementado con permeabilidad de un marcador pertinente. Un aumento de señal en receptor con daño endotelial no demuestra transporte seguro. Los umbrales de viabilidad, integridad y recuperación no se fijan copiando cifras de otros recubrimientos: se validan en controles de esta plataforma. Se conserva su estado pendiente en el protocolo.

## 8. Datos que hay que incorporar

Registrar carga/liberación, revisión de identidad/ruta y fuentes en los contratos ya existentes: `data/candidate_reviews.csv` y `data/formulation_measurements.csv`. Conservar los archivos instrumentales originales para MTT, controles, ensayo ortogonal, carga, liberación y transporte; no completar el registro sólo con una IC50 final.

El antecedente de oxaliplatino requiere sus curvas crudas, carga del lote y controles para reinterpretar los 15 ppm. Hasta entonces no se usa como equivalencia de dosis ni como prueba de superioridad de las nuevas formulaciones.

## Reproducción

```bash
.venv/bin/python scripts/prepare_pilot.py
```

La propuesta es reproducible, revisable y modificable sin borrar versiones. Los campos pendientes del protocolo se cierran con las decisiones y mediciones correspondientes; este documento no firma en nombre del laboratorio.
