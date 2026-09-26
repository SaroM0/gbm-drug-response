"""Modelos por fármaco: una fila por línea, PCA y ajuste sólo dentro de cada fold."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.decomposition import PCA
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import ElasticNet
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .data import load_config, sha256, quality_flags, read_responses, FILES
from .modeling import (attach_split, get_split, grouped_folds, fit_candidates,
                       metrics, runtime_info, cluster_bootstrap)


class TopVarianceGenes(TransformerMixin, BaseEstimator):
    """Filtro aprendido en train; evita PCA sobre miles de columnas constantes."""
    def __init__(self, max_features=2000):
        self.max_features = max_features

    def fit(self, X, y=None):
        variance = np.var(np.asarray(X, dtype=float), axis=0)
        eligible = np.flatnonzero(variance > 1e-12)
        if not len(eligible):
            raise ValueError("No hay genes variables en el fold de entrenamiento.")
        self.selected_ = eligible[np.argsort(variance[eligible], kind="stable")[-self.max_features:]]
        return self

    def transform(self, X):
        return np.asarray(X)[:, self.selected_]


class BoundedPCA(TransformerMixin, BaseEstimator):
    """Limita componentes usando únicamente dimensiones del fold de entrenamiento."""
    def __init__(self, n_components=5, random_state=42):
        self.n_components = n_components
        self.random_state = random_state

    def fit(self, X, y=None):
        n = min(self.n_components, X.shape[0] - 1, X.shape[1])
        if n < 1:
            raise ValueError("No hay muestras suficientes para PCA.")
        self.pca_ = PCA(n_components=n, svd_solver="randomized", random_state=self.random_state).fit(X)
        return self

    def transform(self, X):
        return self.pca_.transform(X)


def load_expression(root):
    root = Path(root)
    path = root / "data/external/expression.csv"
    meta_path = root / "data/external/expression_metadata.json"
    if not path.exists() or not meta_path.exists():
        return None, {"status": "missing_expression", "required_files": [str(path), str(meta_path)]}
    meta = json.loads(meta_path.read_text())
    required = {"source_url", "release", "scale", "gene_id_type", "preprocessing", "sha256"}
    if not required <= meta.keys() or any(not meta[k] for k in required):
        raise ValueError(f"Metadatos incompletos: se requieren {sorted(required)}.")
    if meta["scale"] not in {"log2_tpm_plus1", "rma_log2"}:
        raise ValueError("Proporcione expresión ya normalizada en log2(TPM+1) o RMA log2; no counts crudos.")
    if meta["sha256"] != sha256(path):
        raise ValueError("El checksum de expression.csv no coincide con los metadatos.")
    with path.open(newline="") as f:
        header = next(csv.reader(f))
    if len(header) != len(set(header)):
        raise ValueError("Genes duplicados en cabecera. Resolverlos con una regla documentada.")
    expr = pd.read_csv(path, dtype={"SANGER_MODEL_ID": "string"})
    if "SANGER_MODEL_ID" not in expr or expr.SANGER_MODEL_ID.isna().any() or expr.SANGER_MODEL_ID.duplicated().any():
        raise ValueError("Se requiere una única fila por SANGER_MODEL_ID, sin IDs vacíos.")
    if not expr.SANGER_MODEL_ID.str.fullmatch(r"SIDM\d+").all():
        raise ValueError("IDs no son SIDM: mapear COSMIC/DepMap de forma explícita, sin fuzzy matching.")
    genes = [c for c in expr if c.startswith("gene_")]
    if len(genes) < 2 or set(expr.columns) != {"SANGER_MODEL_ID", *genes}:
        raise ValueError("Columnas permitidas: SANGER_MODEL_ID y genes con prefijo gene_.")
    for gene in genes:
        expr[gene] = pd.to_numeric(expr[gene], errors="raise")
    values = expr[genes].to_numpy(dtype=float)
    if np.isinf(values).any() or np.isnan(values).all(axis=1).any():
        raise ValueError("Expresión infinita o muestra completamente vacía.")
    if np.nanmin(values) < 0:
        raise ValueError("Expresión negativa: revisar que no esté centrada/normalizada con toda la cohorte.")
    return expr, meta


def molecular_candidates(cfg):
    def pipeline(regressor):
        return Pipeline([
            ("impute", SimpleImputer(strategy="median", keep_empty_features=True)),
            ("genes", TopVarianceGenes(max_features=2000)),
            ("scale", StandardScaler()),
            ("pca", BoundedPCA(random_state=cfg["seed"])),
            ("pc_scale", StandardScaler()),
            ("regressor", regressor),
        ])
    return {
        "drug_mean": (DummyRegressor(strategy="mean"), {}),
        "pca_elastic_net": (pipeline(ElasticNet(max_iter=20000)),
                            {"pca__n_components": [5, 10], "regressor__alpha": [0.1, 1.0], "regressor__l1_ratio": [0.1, 0.5]}),
        "pca_extra_trees": (pipeline(ExtraTreesRegressor(n_estimators=150, random_state=cfg["seed"], n_jobs=2)),
                            {"pca__n_components": [5, 10], "regressor__min_samples_leaf": [3, 5]}),
    }


def train_molecular(root):
    root = Path(root)
    cfg = load_config(root)
    expr, meta = load_expression(root)
    folder = root / "reports/molecular"
    folder.mkdir(parents=True, exist_ok=True)
    if expr is None:
        status = {**meta, "message_es": "No se entrenó un modelo molecular: falta expresión real. El notebook 02 sí funciona con los adjuntos."}
        (folder / "status.json").write_text(json.dumps(status, indent=2, ensure_ascii=False) + "\n")
        return status, pd.DataFrame(), pd.DataFrame()
    df = read_responses(root / "data/raw" / FILES["glioma_csv"])
    df = df.loc[quality_flags(df, cfg["qc_rmse_max"]).keep].copy()
    genes = [c for c in expr if c.startswith("gene_")]
    # Reusar exactamente el holdout de metadatos; una unión many-to-one validada.
    joined = df.merge(expr, on="SANGER_MODEL_ID", how="inner", validate="many_to_one")
    original_lines = set(df.SANGER_MODEL_ID)
    matched_lines = set(joined.SANGER_MODEL_ID)
    pd.DataFrame({"SANGER_MODEL_ID": sorted(original_lines - matched_lines)}).to_csv(folder / "unmatched_lines.csv", index=False)
    if not len(joined):
        raise ValueError("Ningún SIDM coincide con respuestas. Revisar versión e identificadores.")
    dev, test = attach_split(joined, get_split(root, cfg))
    # La selección inicial depende sólo de cobertura en desarrollo, nunca de y o desempeño test.
    coverage = dev.groupby("DRUG_ID").SANGER_MODEL_ID.nunique().rename("n_development_lines").reset_index()
    coverage["id_numeric"] = coverage.DRUG_ID.astype(int)
    coverage = coverage.sort_values(["n_development_lines", "id_numeric"], ascending=[False, True])
    eligible = coverage.loc[coverage.n_development_lines >= cfg["molecular_min_development_lines"], "DRUG_ID"].tolist()
    requested = [str(x) for x in cfg["molecular_drug_ids"]]
    drugs = requested if requested else eligible[:cfg["molecular_max_drugs"]]
    coverage.to_csv(folder / "development_drug_coverage.csv", index=False)
    all_scores, all_predictions, all_cv, skipped, models, runs = [], [], [], [], {}, []
    for drug in drugs:
        tr = dev.loc[dev.DRUG_ID == drug].reset_index(drop=True)
        te = test.loc[test.DRUG_ID == drug].reset_index(drop=True)
        if drug not in eligible or te.SANGER_MODEL_ID.nunique() < cfg["molecular_min_test_lines"]:
            skipped.append({"DRUG_ID": drug, "n_dev": len(tr), "n_test": len(te), "reason": "cobertura insuficiente; no reemplazar buscando mejor desempeño"})
            continue
        assert not tr.SANGER_MODEL_ID.duplicated().any() and not te.SANGER_MODEL_ID.duplicated().any()
        folds = grouped_folds(tr, cfg["cv_folds"])
        fitted, cv, details = fit_candidates(molecular_candidates(cfg), tr[genes], tr.LN_IC50, folds)
        selected = cv.iloc[0].model
        pred = te[["NLME_CURVE_ID", "SANGER_MODEL_ID", "DRUG_ID", "DRUG_NAME", "LN_IC50"]].copy()
        cv.insert(0, "DRUG_ID", drug)
        all_cv.append(cv)
        details.to_csv(folder / f"cv_search_{drug}.csv", index=False)
        assignment = tr[["NLME_CURVE_ID", "SANGER_MODEL_ID"]].copy()
        assignment["validation_fold"] = -1
        for i, (_, valid) in enumerate(folds):
            assignment.loc[valid, "validation_fold"] = i
        assignment.to_csv(folder / f"cv_assignments_{drug}.csv", index=False)
        for name, model in fitted.items():
            pred[name] = model.predict(te[genes])
            all_scores.append({"DRUG_ID": drug, "DRUG_NAME": te.DRUG_NAME.iloc[0], "model": name,
                               "selected_by_cv": name == selected, **metrics(te.LN_IC50, pred[name])})
        all_predictions.append(pred)
        models[drug] = {"model": fitted[selected], "gene_columns": genes, "selected": selected}
        runs.append({"DRUG_ID": drug, "selected": selected, "n_dev": len(tr), "n_test": len(te),
                     "interval": cluster_bootstrap(pred, selected, repetitions=cfg["bootstrap_repetitions"], seed=cfg["seed"])})
    scores = pd.DataFrame(all_scores)
    predictions = pd.concat(all_predictions, ignore_index=True) if all_predictions else pd.DataFrame()
    scores.to_csv(folder / "test_metrics.csv", index=False)
    predictions.to_csv(folder / "test_predictions.csv", index=False)
    if all_cv:
        pd.concat(all_cv, ignore_index=True).to_csv(folder / "cv_summary.csv", index=False)
    status = {"status": "trained" if models else "insufficient_coverage", "config": cfg,
              "expression_metadata": meta, "versions": runtime_info(), "n_matched_lines": len(matched_lines),
              "n_unmatched_lines": len(original_lines - matched_lines), "n_gene_features": len(genes),
              "selected_drugs_from_development": drugs, "skipped": skipped, "runs": runs,
              "training_partition": "development_only", "scope": "Exploratorio in vitro, por fármaco conocido; no validación clínica ni externa"}
    (folder / "status.json").write_text(json.dumps(status, indent=2, ensure_ascii=False, default=int) + "\n")
    if models:
        joblib.dump({"models": models, "run": status}, root / "models/molecular.joblib")
    return status, scores, predictions
