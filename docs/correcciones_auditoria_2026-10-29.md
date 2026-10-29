# Correcciones de la auditoría del 29 de octubre de 2026

Se corrigieron los defectos de software y de interpretación que podían resolverse con los archivos disponibles. La ausencia de mediciones experimentales permanece como un requisito explícito: no se sustituyó por supuestos favorables o aprobaciones ficticias. El [informe de auditoría](auditoria_2026-10-29.md) y su evidencia conservan el estado anterior; los [resultados actuales](resultados_actuales.md) se generan desde las corridas corregidas.

## Seguimiento por hallazgo

| Hallazgo | Corrección implementada | Verificación / límite restante |
|---|---|---|
| A01: sulfamato clasificado como amina | SMARTS excluyen sulfamatos, sulfonamidas y otros entornos no amina; el ácido borónico deja de ser asa genérica | Pruebas con pevonedistat real y grupos de control; su ruta pasa a derivatización, pendiente de revisión química |
| A02: capacidades de carga generalizadas | Configuración versionada de escenarios y sensibilidad carga/liberación; fórmula incluye fracción liberada; se eliminan afirmaciones de máximo universal | Prueba de escalamiento y entradas inválidas; capacidad real de cada formulación no medida |
| A03: rango de dosis no rodea IC50 | Cálculo de concentración equivalente mínima/máxima y estado de bracketing por escenario | Prueba donde se alcanza la referencia pero todas las dosis quedan por encima; el rango MTT real requiere piloto |
| A04: conteos y embudo inconsistentes | Conteos marginales separados de embudo acumulativo, con filas DRUG_ID y moléculas en columnas distintas; cifras documentales generadas | 214/280 por encima y 66/280 dentro del rango en T98G; conteos enlazados a run_id |
| A05: mejor ensayo escogido por molécula | Se conservan todos los ensayos, IDs, rangos y discrepancia; referencia conservadora = máximo LN_IC50 de todos los ensayos | Prueba de ensayo discordante; docetaxel conserva ambos y usa 8,23 nM como referencia conservadora |
| A06: rango inferior ignorado | Estado único verifica ambos extremos; no confunde un ajuste bajo MIN_CONC con uno interior | Pruebas de ambos extremos y fracciones por panel; no se borra silenciosamente el dato |
| A07: identidad de sales y estereoquímica | Entidad parental normalizada, forma original separada, ID completo con estereoquímica y PM de entidad activa | Pruebas con sales ionizadas reales, masa equivalente y estereoisómeros; metales y mezclas se preservan para revisión |
| A08: outcome dominado por extrapolación | Desgloses automáticos de métricas por estado de rango sobre los mismos modelos, con interpretación explícita de ajustes exploratorios | No se cambia el test ni se selecciona por subgrupos. Curvas crudas/incertidumbre para un modelo de censura no disponibles; esta limitación biológica no se declara resuelta |
| A09: criterios de comparación/seguridad/transporte incompletos | Contrato del piloto, registro de revisión y registro de mediciones con fuente/checksum; evaluación separada T98G/transporte | Campos pendientes y controles no realizados mantienen estado pendiente. El software comprueba trazabilidad y consistencia, no certifica resultados biológicos |
| A10: resultados sobrescritos y expresión no preparada | Corridas separadas con snapshot de código/configuración/notebooks, versiones, hashes y caché archivada; enlaces actuales sin fallback a métricas/modelos viejos; preparación real de TPM | Pruebas de preservación, fallo posterior, ausencia de expresión y alteración de salidas; matriz de 50 líneas ×37.602 genes y entrenamiento molecular ejecutados |

## Implementación

- `chemistry.py`: normalización, separación de identidades/formas, SMARTS corregidos y carga fisiológica declarada desconocida.
- `screening.py` y `screening.json`: política de potencia, ambos extremos de rango, escenarios y embudo.
- `pilot.py`, `pilot_protocol.json` y registros CSV de `data/`: decisiones y mediciones trazables, sin inventar valores.
- `runs.py`: salidas por corrida, snapshots y verificación; se archivan resultados previos antes de reemplazar rutas cortas.
- `validation.py`: grupos biológicos conocidos y diagnósticos del outcome. Un cruce de donante/familia con el holdout persistido produce error; los grupos relacionados viajan juntos en CV.
- `prepare_expression.py`: lectura completa del miembro TPM histórico, validación de sumas por perfil, transformación una vez y conservación de fuentes/ausentes.
- `update_results_docs.py`: resultados y evidencia versionables desde manifiestos completos verificados. No publica resultados de una ejecución fallida o con archivos alterados.

La selección de genes, imputación, escalado, PCA y ajuste siguen dentro de cada fold. No se añadió una red neuronal ni se cambió la semilla para buscar mejores cifras. La cohorte y las búsquedas del protocolo molecular original se mantienen; su evaluación sigue siendo exploratoria por el holdout previamente observado.

## Qué cambió en la decisión científica

La lista corta pasa de 24 a 23 entidades provisionales con la configuración de referencia. Bortezomib no pasa automáticamente como derivatizable por su ácido borónico; pevonedistat deja de clasificarse como amina. Docetaxel utiliza todos sus ensayos, y su carga teórica de referencia pasa de aproximadamente 0,0166% a 0,0266% por el resumen conservador. No se presenta ninguna molécula como comprada, químicamente aprobada o formulación validada.

La expresión real ya permite ejecutar el modelo molecular para los tres fármacos elegidos por cobertura en desarrollo. El resumen generado muestra sólo la selección por CV para la comparación principal, aunque otro modelo tenga menor error de test. Los intervalos y los desgloses de extrapolación acompañan las cifras; no se interpreta una mejora exploratoria como eficacia antitumoral o transporte confirmado.

## Validación

La suite incluye 48 pruebas del protocolo y de las regresiones detectadas: química, sales/estereoquímica, todos los ensayos, ambos límites de rango, escenarios de liberación/dosis, persistencia sin métricas obsoletas, parentesco/donantes, preparación de expresión y evidencia experimental. Los notebooks 01–03 se ejecutan completos con los datos reales; los resultados y run_id exactos quedan en la documentación generada y `docs/evidencia_correcciones_2026-10-29.json`.

## Requisitos que siguen correspondiendo al laboratorio

`pilot_protocol.json` mantiene nulos el tiempo principal, mejora relevante, número de experimentos independientes y su justificación, modelo de análisis y criterios de seguridad/transporte que no fueron acordados. Los registros de revisión y mediciones están inicialmente vacíos. Estos campos no detienen la caracterización exploratoria; impiden confundir una lista computacional con un ensayo comparativo validado.

Para la comparación formal se necesitan carga y liberación por lote, actividad retenida, equivalencia de dosis y curvas que rodeen la transición en MTT con sus controles. Para atribuir transporte se deben identificar analito, integridad de barrera, seguridad endotelial y balance de masa. Las curvas crudas GDSC necesarias para resolver plenamente el outcome extrapolado tampoco están en los adjuntos. Ningún cambio de software puede producir esas observaciones.
