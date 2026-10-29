# Identidad química y priorización de candidatos

Este documento describe cómo se construye `data/processed/compounds.csv` y cómo se ordenan los candidatos a conjugar en Fe₃O₄@SiO₂-APTES para ensayarlos en T98G. Corresponde a las etapas 1 y 2 del [plan de investigación](plan_investigacion_t98g_nanoparticulas.md), bajo la definición experimental que fijó el laboratorio (sección 0 de ese plan).

**Qué es y qué no es.** Es un embudo de cribado que combina potencia del fármaco libre ya medida en GDSC2, factibilidad de dosis en el rango de nanopartícula que usa el laboratorio y asas químicas presentes en la estructura. No predice la citotoxicidad de la formulación cargada: con aproximadamente tres candidatos no existe conjunto de entrenamiento de formulaciones. Las preguntas C y D del plan se responden en el laboratorio.

## 1. Identidad química

`scripts/build_compounds.py` toma los nombres de fármaco de los adjuntos de respuesta y los resuelve contra PubChem mediante su API REST PUG. Cada respuesta cruda se guarda en `data/external/pubchem_cache/`, de modo que la tabla es reproducible sin volver a consultar la red y queda registro de qué devolvió el servicio y cuándo.

El `compound_id` es el primer bloque del InChIKey, que codifica la conectividad molecular. Así, sales y alias del mismo esqueleto comparten identidad, que es lo que pide el plan antes de contar moléculas.

### Límites de esta resolución

La correspondencia nombre → CID es ambigua por construcción y **no sustituye la verificación de un químico**. El informe `reports/chemistry/resolution_review.csv` recoge las filas que requieren revisión manual:

- nombres que devuelven más de un CID;
- nombres cuyo título en PubChem no coincide con el consultado;
- estructuras que RDKit no puede interpretar o que no tienen estructura asociada.

Los códigos de desarrollo del tipo `AZD7762` o `GSK343` suelen resolverse, pero un nombre comercial o un alias interno puede apuntar a una sal, a un isómero o a otra molécula. Antes de comprar un compuesto hay que confirmar la identidad contra el catálogo del proveedor.

## 2. Anotación de conjugabilidad

`src/gbm_drug_response/chemistry.py` anota cada estructura con patrones SMARTS. La anotación describe **grupos funcionales presentes en la estructura**, no accesibilidad estérica, reactividad real ni eficiencia de carga. Esas tres cosas se determinan experimentalmente.

Cuando la estructura tiene varios fragmentos, como una sal, los grupos se evalúan sobre el fragmento mayor y la fila se marca como `salt_or_mixture`.

### Rutas de unión

La superficie APTES aporta aminas primarias que a pH fisiológico están protonadas, es decir, con carga positiva. De ahí salen dos mecanismos distintos, que el laboratorio ya describió: conjugación covalente y adsorción electrostática. Las rutas se asignan en este orden, que refleja **cuánto trabajo químico adicional exige cada opción, no su eficacia**:

| Ruta | Criterio estructural | Lectura |
|---|---|---|
| `covalente_cooh_edc_nhs` | Carboxilo libre, no esterificado | Caso directo de EDC/NHS con las aminas del APTES |
| `coordinacion_metal` | Contiene Pt, Ru, Pd, Au, Gd o As | Química de coordinación propia; aquí cae el precedente de oxaliplatino |
| `adsorcion_electrostatica` | Otro grupo aniónico: sulfonato, fosfato, tetrazol, acilsulfonamida | Adsorción sobre superficie protonada |
| `covalente_amina_requiere_activacion` | Amina alifática primaria o secundaria | Exige activar la superficie, por ejemplo con glutaraldehído |
| `covalente_requiere_derivatizacion` | Hidroxilo alifático, fenol, tiol o ácido borónico | Requiere un paso previo, como succinilación, que genere el carboxilo |
| `asa_aromatica_poco_reactiva` | Sólo amina unida a anillo aromático | Las aminopirimidinas y similares son poco nucleófilas; no se acoplan en condiciones EDC/NHS habituales |
| `sin_asa_reactiva` | Ninguno de los anteriores | Candidato difícil según el criterio del laboratorio |

### Banderas de revisión

Se registran por separado porque no descalifican un candidato, pero exigen criterio químico:

- `hidrolisis`: éster, lactona, anhídrido, β-lactama, epóxido, haluro de acilo o hidrazona. Son pasivos de estabilidad acuosa durante carga, almacenamiento o liberación.
- `unico_grupo_ionizable`: el carboxilo que se conjugaría es el único grupo ionizable de la molécula. Conjugarlo altera cómo interactúa con su diana. **No es un veredicto de pérdida de actividad**; es exactamente el caso desfavorable que describió el laboratorio y necesita revisión.
- `cationico_repele_superficie`: la molécula queda catiónica a pH fisiológico y no tiene grupos aniónicos. La superficie APTES también está protonada, de modo que se repelen en lugar de adsorberse.
- `logp_alto`, `sal_o_mezcla`, `nombre_no_coincide_con_pubchem`.

La carga se aproxima contando grupos ionizables, **no calculando pKa**. Un grupo límite puede quedar mal clasificado.

## 3. Factibilidad de dosis

Es el criterio que más candidatos descarta y el que usa directamente las condiciones del laboratorio. El laboratorio trabaja entre 5 y 25 ppm de **masa de nanopartícula**. Para que un fármaco alcance su IC50 libre a la dosis más alta de esa serie hace falta una carga en peso de:

```text
carga requerida (% p/p) = IC50(µM) × PM / 250
```

donde 250 sale de 25 µg/mL de nanopartícula convertidos a porcentaje. La fórmula está en `required_loading_pct()` y verificada en las pruebas.

Como referencia, la literatura reporta para esta superficie capacidades de 30 a 52 µg de fármaco por mg de nanopartícula, es decir **3 a 5,2 % p/p**. El cribado usa 5 % como carga alcanzada en la práctica y 10 % como techo optimista, el doble de lo mejor publicado.

**Supuestos, que son fuertes.** Asume liberación completa y que el fármaco cargado actúa con la misma potencia que el libre. Ignora el efecto citotóxico propio de la nanopartícula y cualquier cambio de captación celular. Es una prueba de orden de magnitud para descartar imposibles, no una predicción de respuesta.

## 4. Niveles de prioridad

El ranking es **lexicográfico por niveles**, no una suma ponderada. Una puntuación única mezclaría potencia y conjugabilidad, que no son intercambiables, y sugeriría una precisión que estos datos no sostienen.

Cada nivel exige cumplir los requisitos de los anteriores:

| Nivel | Significado |
|---|---|
| 1 | IC50 medida dentro del rango ensayado, dosis alcanzable y carboxilo libre para EDC/NHS |
| 2 | Igual, con coordinación metálica o adsorción electrostática |
| 3 | Igual, con amina alifática; exige activar la superficie |
| 4 | Igual, con asa que requiere derivatización previa |
| 5 | Dosis alcanzable pero sin asa reactiva práctica |
| 7 | La carga requerida supera el 10 % p/p: no alcanza su IC50 a 25 ppm |
| 8 | IC50 extrapolada por encima del máximo ensayado, o sin peso molecular |
| 9 | Estructura no resuelta o ambigua |

La lista corta son los niveles 1 a 4. Dentro de cada nivel se ordena por carga requerida y luego por potencia.

El primer filtro, y el más importante, es que la IC50 esté **dentro del rango ensayado**. El 76,5 % de las curvas de GBM/glioma tienen la IC50 extrapolada por encima de la concentración máxima: son cotas, no mediciones, y no sirven para elegir un candidato que se va a comprar.

## 5. Resultados de la primera ejecución

### Identidad química

De 295 `DRUG_ID` y 286 nombres en el subconjunto GBM/glioma se resolvieron **243 estructuras**, que corresponden a **233 identidades químicas únicas**. Quedan 59 filas en el informe de revisión.

Los 52 nombres sin estructura son casi todos códigos internos de GDSC sin estructura pública, del tipo `IAP_5620`, `CDK9_5576` o `BDOCA000347a`. No es un fallo de la consulta: esos compuestos no están divulgados y no pueden evaluarse ni comprarse.

### Embudo sobre los 280 fármacos ensayados en T98G

| Filtro | Quedan |
|---|---:|
| Fármacos ensayados en T98G | 280 |
| Con estructura química resuelta | 226 |
| Con IC50 **medida** dentro del rango ensayado | 66 |
| Con dosis alcanzable a 25 ppm (carga ≤ 10 % p/p) | 55 |
| Con asa química práctica: **lista corta** | **24** |

El filtro decisivo es el primero: 174 de los 280 fármacos tienen la IC50 extrapolada por encima de la concentración máxima ensayada. Son cotas inferiores de resistencia, no mediciones, y no sirven para elegir un compuesto que se va a comprar y conjugar.

Otros 18 compuestos tienen potencia y dosis alcanzables pero ninguna asa reactiva práctica. Entre ellos están los dos más potentes de toda la línea, Romidepsina (IC50 4,1 nM) y Daporinad (4,9 nM). Son justamente el caso desfavorable que describió el laboratorio.

### Reparto de rutas en la lista corta

| Ruta | Candidatos |
|---|---:|
| `covalente_cooh_edc_nhs` | 2 |
| `covalente_amina_requiere_activacion` | 5 |
| `covalente_requiere_derivatizacion` | 17 |

**Hallazgo principal, y conviene discutirlo con el laboratorio.** La química que el laboratorio prefiere, acoplamiento directo EDC/NHS sobre un carboxilo libre, sólo encaja con **2 de los 280 fármacos**: Mycophenolic acid y OSI-027. En el conjunto completo hay 15 compuestos con carboxilo libre, pero el resto no supera los filtros de potencia medida o de dosis alcanzable. Además, en ambos el carboxilo es su **único grupo ionizable**, de modo que conjugarlo altera precisamente el grupo que media la interacción con la diana. Es el riesgo que el propio laboratorio anticipó.

La consecuencia práctica es que la vía más probable no es el acoplamiento directo, sino **derivatizar primero**, por ejemplo succinilando un hidroxilo para generar el carboxilo y después aplicar EDC/NHS. Esa ruta cubre 17 de los 24 candidatos.

### Candidatos propuestos para discutir

Los tres cubren rutas químicas distintas, para no apostar todo a un mecanismo de unión:

| Candidato | IC50 libre en T98G | Carga requerida a 25 ppm | Ruta | Por qué |
|---|---:|---:|---|---|
| **Docetaxel** | 5,1 nM | 0,017 % | Derivatización de hidroxilo | Margen de dosis enorme: tres órdenes de magnitud por debajo del límite. La conjugación de taxanos por succinilación está establecida. Vigilar el éster |
| **Gemcitabina** | 349 nM | 0,37 % | Derivatización | Margen amplio y química de conjugación muy documentada. Sin banderas de revisión |
| **OSI-027** | 3,6 µM | 5,8 % | Carboxilo libre, EDC/NHS directo | Única opción, junto con Mycophenolic acid, para la química preferida del laboratorio. Carga por encima de lo reportado y carboxilo como único grupo ionizable: es el candidato de mayor riesgo y mayor valor informativo |

Es una propuesta para decidir con un químico, no una selección cerrada. Tres advertencias sobre la lista:

- **Bortezomib** aparece en el puesto 8 por potencia (5,1 nM) y margen, pero su ácido borónico es el grupo que inhibe el proteasoma. Derivatizarlo anularía la actividad. Está marcado y se desaconseja.
- **Estaurosporina** es el más potente con asa amina (9,3 nM), pero es un inhibidor de quinasas de amplio espectro que se usa como reactivo de laboratorio, no como fármaco. Como citotóxico general, es mal candidato frente al criterio de toxicidad endotelial.
- Once candidatos llevan la bandera `cationico_repele_superficie`: su amina queda protonada a pH fisiológico, igual que la superficie APTES. Para ellos la adsorción no es alternativa si falla la vía covalente.

## 6. Lo que este cribado no decide

- **Disponibilidad comercial y precio.** No se puede determinar desde los datos; la columna queda para el laboratorio.
- **Toxicidad endotelial.** hCMEC/D3 no está en GDSC. Se incluye `t98g_vs_panel_ln`, que compara T98G con la línea mediana del panel completo de GDSC2, como señal de selectividad relativa. Un valor cercano a cero sugiere un citotóxico general, lo que es una precaución a considerar frente al criterio de hCMEC/D3, **no una medida de toxicidad en células sanas**.
- **Comparabilidad con MTT.** Las IC50 de GDSC2 vienen de otro ensayo y otro protocolo. Ordenan candidatos por potencia relativa; no fijan las dosis del laboratorio.
- **Si la conjugación conserva la actividad.** Ninguna anotación estructural lo determina. Es la pregunta que decide el éxito y se responde en el ensayo.

## 7. Reproducción

```bash
.venv/bin/python scripts/build_compounds.py   # identidad química; usa caché local
.venv/bin/python scripts/rank_candidates.py   # potencia, factibilidad y niveles
.venv/bin/python -m pytest tests/test_chemistry.py -q
```

Salidas: `data/processed/compounds.csv`, `reports/chemistry/` y `reports/candidates/`.

## 8. Fuentes

1. PubChem PUG REST, resolución de identidad química: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest
2. RDKit, perfilado de subestructuras con SMARTS: https://www.rdkit.org/docs/RDKit_Book.html
3. Capacidad de carga en Fe₃O₄@SiO₂ funcionalizado con APTES, 51,86 µg/mg con doxorrubicina: https://pubs.rsc.org/ra/article/16/15/13708/1230727/APTES-functionalized-Gd0-18Fe2-82O4-SiO2
4. Capacidad de carga de 34,6 mg/g en Fe₃O₄@SiO₂ con doxorrubicina: https://www.sciencedirect.com/science/article/abs/pii/S0928493119348726
5. Nanotransportadores Fe₃O₄@SiO₂ para entrega de antitumorales: https://pubmed.ncbi.nlm.nih.gov/38797254
6. GDSC, definiciones de respuesta y control de calidad: https://depmap.sanger.ac.uk/documentation/datasets/drug-sensitivity/
