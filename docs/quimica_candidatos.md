# Identidad química y priorización de candidatos — protocolo corregido

El resultado es una priorización exploratoria de fármacos libres para estudiar formulaciones Fe₃O₄@SiO₂-APTES. Las cifras vigentes se generan en [resultados_actuales.md](resultados_actuales.md); el [informe de auditoría](auditoria_2026-10-29.md) conserva el estado anterior a las correcciones.

## Identidad y grupos funcionales

`scripts/build_compounds.py` conserva la respuesta de PubChem, la estructura original y su identidad de forma (`form_inchikey`). Para compuestos orgánicos se prepara una entidad parental mediante `FragmentParent` y `Uncharger` de RDKit. `compound_id` es el InChIKey completo calculado sobre esa entidad: preserva información estereoquímica. El primer bloque se conserva únicamente como `connectivity_family_id`, no como identidad única de molécula. La masa molecular empleada en el escenario corresponde a la entidad activa; la forma comercial y su conversión de masa deben verificarse por separado.

Los complejos metálicos se conservan completos. Las mezclas con varios fragmentos orgánicos sustanciales se marcan para revisión y no pasan automáticamente como una entidad activa normalizada. Es una regla conservadora para revisión, no una identificación inequívoca de sal/mezcla. No se alteran estructuras metálicas para inventar un grupo conjugable en un ligando desconectado.

Los SMARTS de aminas excluyen amidas, sulfonamidas, sulfamatos y entornos fosforilados. Pevonedistat tiene un sulfamato, que ya no cuenta como amina alifática: la anotación propone revisar sus hidroxilos, sin demostrar una ruta experimental. El ácido borónico no cuenta como asa genérica de derivatización. Los grupos son observaciones estructurales; accesibilidad, reactividad, conservación del farmacóforo y actividad tras conjugación siguen requiriendo revisión.

El balance de grupos ionizables es una heurística (`ionizable_group_balance`), **no una carga a pH 7,4**. No se infiere repulsión inevitable ni imposibilidad de adsorción. Se registran necesidades de pKa, potencial zeta y medio experimental. La existencia de un único carboxilo no demuestra que sea el único grupo ionizable o que medie la unión a la diana.

`reports/chemistry/resolution_review.csv` incluye nombres ambiguos, títulos diferentes, búsquedas modificadas, ausencia de estructura y normalización pendiente. Un título diferente puede ser un sinónimo legítimo; no se aprueba automáticamente. Una consulta sin resultado sólo significa que la consulta no resolvió el nombre: no demuestra inexistencia pública o falta de disponibilidad comercial.

## Potencia: todos los ensayos de la misma entidad

El script conserva cada curva, QC, rango, DRUG_ID, AUC y comparación descriptiva con el panel en `reports/candidates/assays.csv`. No convierte esa comparación con otras líneas tumorales en seguridad endotelial.

El ranking contiene una fila por entidad y usa, como resumen conservador explícito, el **máximo LN_IC50 entre todos sus ensayos**. Guarda todos los DRUG_ID, mínimo/máximo de IC50 y discrepancia relativa. La molécula sólo satisface el criterio de rango cuando todos sus ensayos satisfacen `log(MIN_CONC) ≤ LN_IC50 ≤ log(MAX_CONC)`. Los casos mixtos y las extrapolaciones en cualquier dirección se conservan con su estado y se separan de la lista provisional. No se escoge el ensayo más favorable ni se ocultan los rechazados.

Este resumen es una decisión de cribado, no un estimador de consenso biológico ni una garantía de independencia entre ensayos. Los ajustes extrapolados no se convierten automáticamente en cotas medidas de resistencia.

## Dosis, carga y liberación como escenarios

Las reglas y escenarios están en [screening.json](../screening.json). Su versión, valores y hashes se guardan en cada corrida. Para NP en µg/mL, concentración libre de referencia en µM y PM de entidad activa en g/mol:

```text
carga requerida (% masa fármaco / masa portador)
  = IC50_referencia × PM / (10 × dosis_NP × fracción_liberada)

fármaco equivalente (µM)
  = dosis_NP × carga_porcentual / 100 × fracción_liberada × 1000 / PM
```

Son escenarios de masa que suponen la misma potencia que el fármaco libre, sin cambios de captación ni actividad del portador. No predicen la potencia de la nanoformulación. Se evalúan cargas de 1, 3, 5 y 10%, y liberaciones de 0,1, 0,5 y 1; el escenario de referencia es 5%/1. El corte de cribado de 10% es una elección exploratoria, **no una frontera física ni el máximo publicado**. Cambiarlo requiere identificar una nueva corrida, sin borrar las anteriores.

Las fuentes citadas anteriormente son ejemplos específicos: doxorrubicina sobre ferrita dopada con Gd/sílice/APTES y doxorrubicina sobre sílice carboxilada. No establecen capacidad universal para todos los fármacos, nanopartículas o mecanismos de unión. [Estudio APTES/Gd](https://pubs.rsc.org/en/content/articlehtml/2026/ra/d6ra00215c), [estudio carboxilado](https://researchnow.flinders.edu.au/en/publications/modified-green-synthesis-of-fesub3subosub4subsiosub2sub-nanoparti/).

`dose_scenarios.csv` distingue alcanzar la referencia a la dosis máxima de rodearla con el rango completo. Devuelve `brackets_reference`, `all_doses_above_reference` o `all_doses_below_reference`; ninguna categoría certifica que MTT produzca esa IC50. `scenario_summary.csv` publica la sensibilidad sin escoger a posteriori el escenario más favorable.

## Niveles, conteos y decisiones

Los niveles 1–4 describen rutas estructurales con distinto trabajo químico: COOH, coordinación/adsorción, activación de amina y derivatización. El nivel 5 indica ausencia de asa práctica detectada; 7 supera el escenario elegido; 8 indica rango/base de masa no evaluable; 9 indica estructura pendiente. El orden no expresa probabilidad de éxito, eficacia ni una decisión de compra.

`funnel.csv` aplica filtros acumulativos e informa en columnas distintas filas DRUG_ID y entidades normalizadas. Los conteos marginales de rango se informan aparte. La lista `shortlist.csv` sigue siendo provisional, incluso cuando un compuesto tiene buena potencia y ruta estructural.

Antes de usar un candidato para una comparación formal:

1. Registrar identidad, forma comercial, ruta, revisión del farmacóforo y responsable en [candidate_reviews.csv](../data/candidate_reviews.csv).
2. Completar y aprobar [pilot_protocol.json](../pilot_protocol.json): tiempo y contraste principal, efecto relevante, experimentos independientes y justificación, controles, análisis, seguridad endotelial y definición de transporte.
3. Registrar mediciones y fuente verificable en [formulation_measurements.csv](../data/formulation_measurements.csv): lote, carga sobre entidad activa, liberación, equivalencia de dosis y validación del ensayo/rango. El software verifica el hash del archivo de evidencia, no la veracidad biológica de la medición.

`formulation_readiness.csv` distingue preparación para comparación T98G y soporte para afirmaciones de transporte. Se exigen criterios de barrera íntegra, seguridad endotelial y balance de masa para lo segundo. Los registros iniciales están vacíos: no se fabrican aprobaciones, réplicas, cargas, liberaciones ni resultados de controles. Esto no impide realizar el piloto de caracterización necesario para obtenerlos.

Docetaxel y gemcitabina siguen siendo hipótesis de formulación; OSI-027 requiere resolver su carga/liberación y conservación de actividad. La terna no representa tres categorías de unión distintas: los dos primeros comparten derivatización. Ningún candidato está confirmado por este cribado.

## Reproducción y trazabilidad

```bash
.venv/bin/python scripts/build_compounds.py
.venv/bin/python scripts/rank_candidates.py
.venv/bin/python scripts/update_results_docs.py
.venv/bin/python -m pytest -q
```

Cada corrida guarda salidas separadas en `reports/runs/<etapa>/<run_id>/`, con `run_manifest.json`, hashes de entradas/salidas, snapshot de código/configuración y versiones. Química archiva además la caché de PubChem. Las rutas cortas son enlaces a la corrida actual; resultados anteriores se conservan, y un fallo no recupera silenciosamente métricas de una corrida previa. Para trasladar una investigación hay que conservar originales, corridas, caché y ZIP: Git sigue excluyendo los archivos grandes.

Fuentes de implementación: [RDKit MolStandardize](https://www.rdkit.org/docs/source/rdkit.Chem.MolStandardize.rdMolStandardize.html), [InChI](https://www.inchi-trust.org/technical-faq/), [química de pevonedistat](https://pubs.acs.org/doi/10.1021/bk-2016-1240.ch001), [definiciones de GDSC](https://depmap.sanger.ac.uk/documentation/datasets/drug-sensitivity/).
