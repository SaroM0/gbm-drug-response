# Mapa del protocolo y trabajo pendiente

Actualización: 9 de noviembre de 2026. Fuentes: protocolo completo facilitado en la conversación, fechado el 26 de septiembre de 2026, y mensajes del equipo del 29 de octubre al 4 de noviembre proporcionados por el usuario. **El laboratorio solo autorizó comprar dos fármacos**. El texto del protocolo fue proporcionado en conversación; su versión vigente y aprobación formal deben confirmarse. Este documento registra su contenido, no certifica su aprobación ni sustituye el protocolo firmado.

## Decisiones y alcance

| Tema | Definición o estado | Fuente |
|---|---|---|
| Objetivo | Priorizar fármacos aprobados y evaluar actividad en T98G y transporte aparente con hCMEC/D3 en Transwell, libres y asociados a IONPs. | Protocolo, pregunta y §5 |
| Compra | Dos fármacos; aún no se han seleccionado ni verificado sus requisitos. Los tres nombres históricos son candidatos de análisis, no una lista de compra. | Decisión comunicada por el usuario; protocolo §6.4.1.6 |
| Selección | Dos candidatos disponibles comercialmente, con mecanismos o vías distintos, identidad y aprobación verificadas; el protocolo propone sensibilidad, selectividad, consistencia y factibilidad química, con diversidad estructural. | §6.4.1.6; ponderaciones pendientes |
| Comparación | Libre frente a IONPs-fármaco con dosis nominal equivalente de principio activo; controles de vehículo, IONPs vacías y sin tratamiento. | §6.4.3–6.4.4 |
| Lecturas | MTT como indicador indirecto de viabilidad; controles de interferencia y contraste por azul de tripano si corresponde. Rodamina indica señal asociada al marcador, no concentración de fármaco ni paso intacto demostrado. | §6.4.3–6.4.4 |
| Réplicas | Tres experimentos independientes por condición y tres técnicas, por factibilidad; piloto exploratorio sin potencia confirmatoria establecida. | §6.2.3 |
| Barrera | Microscopía, TEER y FITC-dextrán; faltan criterios numéricos de aceptación y recuperación. | §6.2.2 y §6.3 |
| Modelo | Predice LN_IC50 de fármacos libres en líneas no vistas; el ranking propuesto utiliza respuesta medida en T98G y criterios adicionales, no predicciones del modelo. | §6.4.1.3 |
| Inferencia | Exploratoria con n=3; el protocolo propone contrastes predefinidos y ajuste Holm. Concretar familias, desenlaces y análisis con el asesor metodológico. | §6.5 |
| Otra línea tumoral | El laboratorio ofreció una fase adicional en una línea de otro tipo de cáncer; Sofía mencionó cáncer gástrico o de seno como ejemplos. La línea y el diseño siguen sin seleccionar. | Mensajes de Sofía, 2 de noviembre |

El protocolo propone tres tiempos (24, 48 y 72 h), dosis equivalentes de principio activo y series adaptadas al candidato. Los 5–25 ppm fueron una propuesta anterior para nanopartículas, no una dosis de fármaco universal ni un rango validado. La plataforma Fe₃O₄@SiO₂-APTES figura en la información anterior del laboratorio; el protocolo facilitado habla de IONPs comerciales con recubrimiento potencialmente variable. Confirmar la plataforma definitiva.

## Trabajo computacional: instrucciones y entregables

- [ ] **Auditar §6.4.1 frente al código y salidas**, criterio por criterio, con evidencia de archivo, fuente y estado: implementado, parcial o pendiente. No dar por implementadas características por estar descritas en el protocolo.
- [ ] **Elegibilidad:** obtener o identificar una fuente trazable de aprobación y estado comercial; resolver identificadores y formas químicas, ambigüedades y fecha/jurisdicción. `clinical_phase = Launched` es el criterio propuesto, no prueba automáticamente aprobación vigente en una jurisdicción determinada.
- [ ] **Ranking:** verificar implementación y definir antes de observar nuevos resultados experimentales las ponderaciones, datos faltantes, agregación de ensayos, selectividad, consistencia y factibilidad de formulación. Revisar diversidad Morgan/Tanimoto y diferencia de mecanismo con evidencia. Documentar exclusiones, incertidumbre y sensibilidad a decisiones de ranking. No confundir Z_SCORE entre líneas tumorales con selectividad frente a tejido sano.
- [ ] **Curvas:** revisar datos crudos, rangos y variación por placa. Una IC50 ajustada fuera del rango es extrapolada; no convertirla automáticamente en una cota observada o censura cierta. Documentar los límites del ajuste y la evidencia de cada candidato.
- [ ] **Comparación con otra línea tumoral:** inventariar líneas de cáncer gástrico y de mama disponibles en GDSC2 y su cobertura de los mismos fármacos y ensayos; documentar identificadores, QC, rangos y comparabilidad antes de proponer opciones al laboratorio. Definir comparaciones por fármaco, dosis y tiempo y sus incertidumbres; una sola línea no demuestra especificidad para glioblastoma ni representa todo otro tipo de cáncer. La disponibilidad experimental de la línea aún debe confirmarse.
- [ ] **Modelos:** informar la referencia por fármaco y modelos de expresión ya entrenados, con su desempeño por fármaco y límites. El test ya fue observado; reutilizarlo tras adaptar decisiones es exploratorio. No afirmar validación independiente ni entrenar un predictor de formulación sin datos adecuados.
- [ ] **Lista para decisión:** producir fichas comparables de elegibles y justificar dos finalistas únicamente después de validar aprobación, mecanismo, identidad, disponibilidad y factibilidad con el equipo. Una estructura compatible no prueba carga, liberación ni actividad.
- [ ] **Piloto:** tras selección, actualizar `pilot_design.json`, número de comparaciones y propuesta de placas; recalcular recursos y escenarios. No escoger arbitrariamente dos de la terna anterior. El comando `scripts/prepare_pilot.py` bloquea el diseño histórico incompatible con el límite actual.
- [ ] **Entrega reproducible:** registrar fuentes, hashes, configuraciones, particiones, limitaciones, resultados y versiones; actualizar metodología y documentación del protocolo cuando el equipo confirme su versión. Mantener las ejecuciones y documentos históricos identificados como tales.

## Pendientes para Sofía y el equipo

1. Confirmar versión vigente y aprobación del protocolo, responsables y proceso de selección de los dos candidatos.
2. Precisar jurisdicción y alcance de aprobación farmacológica; revisar relevancia médica, disponibilidad y si se necesita un comparador de referencia, indicando si ya está disponible y cómo encaja en el presupuesto de compra.
3. Acordar efecto mínimo de interés, criterios de continuación y contraste principal; revisar análisis y multiplicidad con el asesor metodológico.
4. Confirmar nanopartícula, recubrimiento, proveedor, marcaje, ruta de asociación y viabilidad química de cada finalista; responsables y métodos para medir carga y liberación por lote.
5. Definir criterios cuantitativos de aceptación de endotelio, TEER, FITC-dextrán y recuperación, así como método para medir fármaco si se pretende cuantificar su transporte. La rodamina por sí sola no lo cuantifica.
6. Confirmar ejecución real y capacidad para las tres campañas y controles, ensayo ortogonal cuando corresponda, selección de dosis y tiempos y quién aprueba el protocolo experimental.
7. Compartir datos previos de oxaliplatino y del laboratorio, si existen: originales por pozo, lotes, controles, dosis, tiempos, imágenes y archivos instrumentales. Confirmar responsables de química, laboratorio y metodología.
8. Sobre la fase adicional ofrecida: ¿qué líneas no GBM están disponibles, con qué recursos y cronograma, y se ensayarían los mismos dos fármacos, formulaciones, controles y tiempos? ¿Qué diferencia entre líneas considerarán relevante y qué comparación buscan hacer?

## Conversación con Sofía: estado y correcciones

- El 29 de octubre Sofía preguntó por los dos compuestos que se habían identificado con carboxilo libre y comunicó que el laboratorio podía comprar dos fármacos. Se le respondieron ácido micofenólico y OSI-027. Ante la lista inicial de docetaxel, gemcitabina y ácido micofenólico, comentó que inicialmente tenía sentido. Esto no constituye selección final, aprobación química ni confirmación de que todos tengan la misma clasificación terapéutica.
- El 2 de noviembre Sofía expresó una valoración favorable del avance y comunicó la oferta de evaluar otra línea tumoral. El objetivo planteado fue explorar si la respuesta es específica para GBM; con una línea por tipo tumoral la interpretación se limita a esos modelos y condiciones. Esa fase no implica por sí sola selectividad frente a tejido sano.
- El 4 de noviembre Sofía informó que ya había entregado la metodología el día anterior y compartió un [documento editable](https://docs.google.com/document/d/1Nr0gf-mzkA_WO8ICSfPVtsN1Hp26GcwamDtIsP_GAHU/edit?usp=drivesdk). Santiago indicó que dejó comentarios sobre la parte computacional, pendientes de revisión por Sofía. No se ha verificado aquí el estado actual del documento ni de los comentarios; entrega no equivale a aprobación final.
- La conversación inicial contenía cifras y afirmaciones después corregidas. Las 66 IC50 son ajustes dentro del rango, no necesariamente mediciones directas precisas; las 214 extrapolaciones no se convierten automáticamente en cotas de resistencia. El cribado actual informa 23 moléculas provisionales bajo su escenario, no 24 finalistas aprobados. Para docetaxel, la revisión usa ambos ensayos y la carga teórica de referencia pasó de aproximadamente 0,017 % a 0,0266 % bajo los supuestos revisados. Son escenarios, no carga medida ni promesa de liberación. Véanse [resultados actuales](resultados_actuales.md) y [correcciones](correcciones_auditoria_2026-10-29.md).
- No está demostrado que el carboxilo sea el único grupo ionizable de ambos compuestos ni que su conjugación afecte necesariamente la interacción con la diana. La derivatización propuesta no demuestra factibilidad, liberación ni conservación de actividad: requiere revisión experimental específica.
- El mensaje de las 18:59 del 29 de octubre está truncado en el material recibido; no se reconstruye su contenido faltante. Las propuestas posteriores sobre modelos de nanoformulación requieren mediciones suficientes y evaluación independiente; no son un modelo ya validado.

## Discrepancias y archivos históricos

- `pilot_design.json`, [piloto para revisión](piloto_para_revision.md) y la corrida de 54 placas representan una propuesta anterior para tres candidatos, no autorizada como plan final. Sus escenarios no se deben extrapolar como potencia actual de dos candidatos. `pilot_protocol.json` refleja el límite de compra, no la aprobación del resto del protocolo.
- La expresión real ya fue preparada y entrenada para tres fármacos seleccionados por cobertura, no para seleccionar los dos fármacos de compra. El texto facilitado que condiciona su ejecución a la disponibilidad de expresión está desactualizado; ver [resultados actuales](resultados_actuales.md).
- La afirmación del protocolo de que el test se usará una sola vez debe reconciliarse con el hecho de que ya se observó. No se obtiene un test independiente cambiando la semilla.
- Corregir el texto truncado «H1bTEER» en §6.1.1 de la fuente y armonizarlo con la hipótesis anterior, sin asumir una interpretación aprobada.
- Los documentos de auditoría y evidencia mantienen sus fechas, cifras y decisiones originales como registro histórico. Los tres fármacos analizados en expresión o curvas no son necesariamente los dos finalistas.
