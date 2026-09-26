# Protocolo del primer acercamiento

Este protocolo se escribió antes de ejecutar los modelos. Los hiperparámetros se eligen con validación de desarrollo. El test reservado sólo se usa para evaluación final. Cambiar el protocolo después de observar test convierte ese test en desarrollo: una nueva semilla no restaura independencia científica.

## Pregunta y unidad de análisis

Predecir `LN_IC50` de un fármaco conocido en una línea GBM/glioma no utilizada para ajustar el modelo. La unidad de agrupación es `SANGER_MODEL_ID`; las curvas de una línea viajan juntas. Los 51 SIDM son modelos celulares distintos, no 51 pacientes verificados. Para afirmar independencia por donante hay que incorporar metadatos de origen y agrupar también modelos derivados del mismo donante. Las 12.001 curvas no son 12.001 muestras biológicas independientes.

## Outcome y control de calidad

Outcome principal: `LN_IC50`, logaritmo natural de la IC50 expresada en µM. No aplicar otro logaritmo. La regresión preserva la escala continua. No se construye una etiqueta sensible/resistente usando una mediana global.

Se preespecifica `0 <= RMSE <= 0.3`, además de outcome finito, AUC válida, concentraciones positivas ordenadas e identificadores completos. El RMSE de la tabla mide ajuste de la curva experimental; `prediction_rmse` mide error del predictor y tiene otra escala. GDSC documenta que elimina curvas con RMSE > 0.3 antes de publicar [1]. Por ello se espera que el umbral principal no quite registros de estos adjuntos. Los umbrales 0.15, 0.20 y 0.25 se muestran como análisis de retención, sin entrenar para escoger el que genere mejores métricas.

`LN_IC50` fuera de `log(MIN_CONC)`–`log(MAX_CONC)` se marca como extrapolación. No se trunca ni elimina automáticamente: quitar respuestas altas puede sesgar la cohorte hacia sensibilidad. El modelo inicial trata estos ajustes como valores puntuales; evaluar censura/incertidumbre es trabajo posterior. Un RMSE pequeño no garantiza una IC50 bien determinada fuera del rango.

## Predictores y referencias

Notebook 02: media global, media por DRUG_ID, Elastic Net y Extra Trees usando únicamente DRUG_ID y CANCER_TYPE. Es una referencia con metadatos, no un predictor molecular individualizado. Dentro del mismo tipo tumoral y fármaco asignará la misma predicción a varias líneas. Los nuevos DRUG_ID se registran como no soportados y se excluyen de las métricas principales.

Notebook 03: un modelo por DRUG_ID, con una fila por SIDM. `X` contiene expresión basal real; `y` es LN_IC50. Entrena imputación de medianas, selección de hasta 2.000 genes variables, estandarización, PCA de 5/10 componentes y Elastic Net/Extra Trees. Cada transformación se ajusta dentro de cada fold, sin utilizar test [2]. El baseline por fármaco es la media de desarrollo. No se replica la misma expresión cientos de veces antes de ajustar PCA.

`PATHWAY_NAME` describe la diana del medicamento: no es actividad de vías de una célula. Esta versión implementa PCA. Una representación por vías requerirá colecciones de genes y un método de puntuación documentados; no se simula a partir de PATHWAY_NAME.

No entran en X: AUC, RMSE, Z_SCORE, IDs de curva/resultado, nombres/IDs de línea, concentraciones de ensayo. Excluir concentraciones es una elección conservadora del protocolo para evitar que el diseño del ensayo domine la predicción. DRUG_ID sí es admisible para evaluar fármacos conocidos; no permite extrapolar a una molécula nueva.

## Validación

Holdout fijo de 20% de líneas con semilla 42, persistido y reutilizado entre notebooks. Dentro de desarrollo: GroupKFold de 4 particiones [3], con búsqueda pequeña por MAE. Se selecciona el modelo con menor MAE de CV, incluyendo la referencia simple. La media por fármaco se vuelve a aprender dentro de cada fold. El holdout no participa en la búsqueda.

La primera ejecución molecular procesa hasta tres fármacos por cobertura en desarrollo; empates por DRUG_ID numérico. El límite reduce cómputo y no es selección de fármacos terapéuticos. Se puede preespecificar una lista en config.json antes de mirar resultados. Se requieren 25 líneas de desarrollo y 5 de test por fármaco como mínimos operativos para ejecutar, no como cálculo de potencia ni garantía de validez.

Informar MAE, RMSE predictivo, R² y Spearman, además de errores por fármaco y línea. R²/Spearman agrupados pueden ser altos sólo por diferencias de potencia entre fármacos; las métricas por fármaco y la comparación con la media por DRUG_ID son esenciales. El intervalo bootstrap agrupa por SIDM, se condiciona al modelo ya entrenado y no incluye toda la incertidumbre de selección/entrenamiento. Con aproximadamente 11 líneas de test será impreciso. Las correlaciones constantes o insuficientes se reportan como no definidas.

T98G no es validación externa: es parte de los otros adjuntos. Reservarla antes de entrenar podría constituir un test interno de una línea, nunca una cohorte externa. No se usan las tablas repetidas como datos nuevos.

## Fuentes consultadas

1. Sanger, definiciones de respuestas, QC y múltiples DRUG_ID: https://depmap.sanger.ac.uk/documentation/datasets/drug-sensitivity/
2. Scikit-learn, prevención de fuga por preprocesamiento: https://scikit-learn.org/stable/common_pitfalls.html
3. Scikit-learn, GroupKFold: https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.GroupKFold.html
4. Scikit-learn, Elastic Net: https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.ElasticNet.html
5. Sanger, expresión y escalas de RNA-seq: https://depmap.sanger.ac.uk/documentation/datasets/expression/
6. GDSC, propósito de integrar respuesta y características moleculares: https://www.sanger.ac.uk/tool/gdsc-genomics-drug-sensitivity-cancer/
7. Sanger, relaciones entre modelos: https://depmap.sanger.ac.uk/documentation/cell-models/model-relationships/

Consultadas durante la preparación del repositorio. El sitio histórico cancerrxgene.org devuelve HTTP 410 para la descarga RMA probada. La documentación actual y los adjuntos pueden representar versiones distintas; los conteos del proyecto se calculan sobre los adjuntos, no sobre cifras de una página web.

## Siguiente ampliación

Incorporar expresión con correspondencia SIDM verificada; registrar donante, versión, plataforma y lote. Ampliar el entrenamiento a líneas pan-cáncer, manteniendo GBM/glioma de evaluación completamente fuera de selección/preprocesamiento. Para transferencia a pacientes hacen falta datos y validación externos; el comportamiento en cultivo no demuestra eficacia clínica.
