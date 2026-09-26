"""Validación agrupada, referencias y modelos tabulares sin variables de respuesta en X."""
from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import sklearn
from scipy.stats import spearmanr
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.compose import ColumnTransformer
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GridSearchCV, GroupKFold, GroupShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder
from threadpoolctl import threadpool_limits

from .data import load_config, read_responses, FILES, sha256, quality_flags

BLOCKED_FEATURES = {"LN_IC50", "AUC", "RMSE", "Z_SCORE", "NLME_CURVE_ID",
                    "NLME_RESULT_ID", "SANGER_MODEL_ID", "CELL_LINE_NAME",
                    "MIN_CONC", "MAX_CONC"}


class DrugMeanRegressor(RegressorMixin, BaseEstimator):
    """Media por DRUG_ID aprendida exclusivamente en entrenamiento."""
    def fit(self, X, y):
        self.global_mean_ = float(np.mean(y))
        self.drug_means_ = pd.DataFrame({"drug": np.asarray(X["DRUG_ID"]), "y": np.asarray(y)}).groupby("drug").y.mean()
        return self

    def predict(self, X):
        return X["DRUG_ID"].map(self.drug_means_).fillna(self.global_mean_).to_numpy(dtype=float)


def get_split(root, cfg):
    """Partición persistida por SIDM antes de aplicar filtros o incorporar omics."""
    root = Path(root)
    original = read_responses(root / "data/raw" / FILES["glioma_csv"])
    ids = sorted(original.SANGER_MODEL_ID.dropna().unique())
    signature = hashlib.sha256(json.dumps({"ids": ids, "seed": cfg["seed"], "test_size": cfg["test_size"]}, sort_keys=True).encode()).hexdigest()
    folder = root / "data/processed"
    csv_path, sig_path = folder / "line_split.csv", folder / "line_split_signature.txt"
    if csv_path.exists() or sig_path.exists():
        if not csv_path.exists() or not sig_path.exists() or sig_path.read_text().strip() != signature:
            raise ValueError("Partición existente incompatible. Archive los resultados y la partición antes de crear un protocolo nuevo.")
        split = pd.read_csv(csv_path)
        if split.SANGER_MODEL_ID.duplicated().any() or set(split.SANGER_MODEL_ID) != set(ids):
            raise ValueError("Partición alterada o incompleta.")
        if set(split.partition) != {"development", "test"}:
            raise ValueError("Etiquetas de partición inesperadas.")
    else:
        dev, test = next(GroupShuffleSplit(n_splits=1, test_size=cfg["test_size"], random_state=cfg["seed"]).split(ids, groups=ids))
        split = pd.DataFrame({"SANGER_MODEL_ID": ids, "partition": "development"})
        split.loc[test, "partition"] = "test"
        split.to_csv(csv_path, index=False)
        sig_path.write_text(signature + "\n")
    # Recalcular para detectar también cambios manuales de etiquetas con los mismos IDs.
    _, test_idx = next(GroupShuffleSplit(n_splits=1, test_size=cfg["test_size"], random_state=cfg["seed"]).split(ids, groups=ids))
    if set(split.loc[split.partition == "test", "SANGER_MODEL_ID"]) != {ids[i] for i in test_idx}:
        raise ValueError("La asignación de líneas cambió respecto a la semilla registrada.")
    return split


def attach_split(df, split):
    out = df.merge(split, on="SANGER_MODEL_ID", how="left", validate="many_to_one")
    if out.partition.isna().any():
        raise ValueError("Hay líneas sin partición.")
    dev = out.loc[out.partition == "development"].reset_index(drop=True)
    test = out.loc[out.partition == "test"].reset_index(drop=True)
    assert set(dev.SANGER_MODEL_ID).isdisjoint(test.SANGER_MODEL_ID)
    return dev, test


def grouped_folds(df, n_splits):
    if df.SANGER_MODEL_ID.nunique() < n_splits:
        raise ValueError("Número de líneas insuficiente para CV.")
    splits = list(GroupKFold(n_splits=n_splits).split(df, groups=df.SANGER_MODEL_ID))
    for train, valid in splits:
        assert set(df.iloc[train].SANGER_MODEL_ID).isdisjoint(df.iloc[valid].SANGER_MODEL_ID)
    return splits


def metrics(y, pred):
    y, pred = np.asarray(y, dtype=float), np.asarray(pred, dtype=float)
    corr = float(spearmanr(y, pred).statistic) if len(y) >= 3 and np.ptp(y) > 0 and np.ptp(pred) > 0 else np.nan
    return {"n": len(y), "mae": mean_absolute_error(y, pred),
            "prediction_rmse": float(np.sqrt(mean_squared_error(y, pred))),
            "r2": r2_score(y, pred) if len(y) >= 2 and np.ptp(y) > 0 else np.nan,
            "spearman": corr}


def cluster_bootstrap(predictions, model, baseline="drug_mean", repetitions=1000, seed=42):
    """IC exploratorio por línea; delta MAE negativo favorece al modelo."""
    p = predictions.copy()
    p["error"] = np.abs(p.LN_IC50 - p[model])
    p["base_error"] = np.abs(p.LN_IC50 - p[baseline])
    line = p.groupby("SANGER_MODEL_ID")[["error", "base_error"]].mean()
    rng = np.random.default_rng(seed)
    a = line.to_numpy()
    draws = a[rng.integers(0, len(a), size=(repetitions, len(a)))].mean(axis=1)
    delta = draws[:, 0] - draws[:, 1]
    return {"n_test_lines": len(a), "line_macro_mae": float(a[:, 0].mean()),
            "line_macro_mae_ci95": np.quantile(draws[:, 0], [0.025, 0.975]).tolist(),
            "delta_vs_baseline": float((a[:, 0] - a[:, 1]).mean()),
            "delta_ci95": np.quantile(delta, [0.025, 0.975]).tolist(),
            "method": "bootstrap agrupado por SIDM; percentil; condicionado al entrenamiento fijo"}


def fit_candidates(candidates, X, y, folds):
    fitted, cv_tables, summary = {}, [], []
    with threadpool_limits(limits=2):
        for name, (model, grid) in candidates.items():
            search = GridSearchCV(model, grid, scoring="neg_mean_absolute_error", cv=folds,
                                  refit=True, n_jobs=1, error_score="raise", return_train_score=False)
            search.fit(X, y)
            fitted[name] = search.best_estimator_
            table = pd.DataFrame(search.cv_results_)
            table.insert(0, "model", name)
            cv_tables.append(table)
            summary.append({"model": name, "cv_mae": -search.best_score_,
                            "cv_fold_std": search.cv_results_["std_test_score"][search.best_index_],
                            "best_params": json.dumps(search.best_params_, sort_keys=True)})
    return fitted, pd.DataFrame(summary).sort_values("cv_mae"), pd.concat(cv_tables, ignore_index=True)


def metadata_candidates(cfg):
    features = cfg["metadata_features"]
    assert not set(features) & BLOCKED_FEATURES
    encoder = ColumnTransformer([("categories", OneHotEncoder(handle_unknown="ignore", sparse_output=True), features)])
    return {
        "global_mean": (DummyRegressor(strategy="mean"), {}),
        "drug_mean": (DrugMeanRegressor(), {}),
        "elastic_net_metadata": (Pipeline([("encode", encoder), ("regressor", ElasticNet(max_iter=20000, tol=1e-4))]),
                                 {"regressor__alpha": [0.001, 0.01], "regressor__l1_ratio": [0.1, 0.5]}),
        "extra_trees_metadata": (Pipeline([("encode", encoder), ("regressor", ExtraTreesRegressor(n_estimators=120, random_state=cfg["seed"], n_jobs=2))]),
                                 {"regressor__min_samples_leaf": [5, 15], "regressor__max_features": [1.0]}),
    }


def runtime_info():
    return {"python": platform.python_version(), "pandas": pd.__version__, "numpy": np.__version__, "scikit_learn": sklearn.__version__}


def train_metadata(root):
    root = Path(root)
    cfg = load_config(root)
    raw = read_responses(root / "data/raw" / FILES["glioma_csv"])
    df = raw.loc[quality_flags(raw, cfg["qc_rmse_max"]).keep].copy()
    dev, test = attach_split(df, get_split(root, cfg))
    folds = grouped_folds(dev, cfg["cv_folds"])
    features = cfg["metadata_features"]
    X = dev[features].astype(str)
    fitted, cv, details = fit_candidates(metadata_candidates(cfg), X, dev.LN_IC50, folds)
    selected = cv.iloc[0].model  # Congelado antes de consultar resultados test.
    pred = test[["NLME_CURVE_ID", "SANGER_MODEL_ID", "DRUG_ID", "DRUG_NAME", "LN_IC50"]].copy()
    pred["drug_seen_in_training"] = test.DRUG_ID.isin(dev.DRUG_ID)
    for name, model in fitted.items():
        pred[name] = model.predict(test[features].astype(str))
    supported = pred.loc[pred.drug_seen_in_training].copy()
    if supported.empty:
        raise ValueError("No hay fármacos conocidos en test: este protocolo no evalúa fármacos nuevos.")
    scores, perdrug, perline = [], [], []
    for name in fitted:
        scores.append({"model": name, **metrics(supported.LN_IC50, supported[name])})
        for drug, part in supported.groupby("DRUG_ID"):
            perdrug.append({"model": name, "DRUG_ID": drug, "DRUG_NAME": part.DRUG_NAME.iloc[0], **metrics(part.LN_IC50, part[name])})
        for sidm, part in supported.groupby("SANGER_MODEL_ID"):
            perline.append({"model": name, "SANGER_MODEL_ID": sidm, **metrics(part.LN_IC50, part[name])})
    scores, perdrug, perline = pd.DataFrame(scores), pd.DataFrame(perdrug), pd.DataFrame(perline)
    scores["macro_drug_mae"] = scores.model.map(perdrug.groupby("model").mae.mean())
    scores["macro_line_mae"] = scores.model.map(perline.groupby("model").mae.mean())
    folder = root / "reports/metadata_baseline"
    folder.mkdir(parents=True, exist_ok=True)
    for name, table in [("cv_summary", cv), ("cv_search", details), ("test_metrics", scores),
                        ("test_predictions", pred), ("test_by_drug", perdrug), ("test_by_line", perline)]:
        table.to_csv(folder / f"{name}.csv", index=False)
    fold_manifest = dev[["NLME_CURVE_ID", "SANGER_MODEL_ID"]].copy()
    fold_manifest["validation_fold"] = -1
    for i, (_, valid) in enumerate(folds):
        fold_manifest.loc[valid, "validation_fold"] = i
    assert fold_manifest.validation_fold.ge(0).all()
    fold_manifest.to_csv(folder / "cv_assignments.csv", index=False)
    run = {"config": cfg, "versions": runtime_info(), "selected_on_development_cv": selected,
           "n_dev_lines": dev.SANGER_MODEL_ID.nunique(), "n_test_lines": test.SANGER_MODEL_ID.nunique(),
           "unsupported_test_rows": int((~pred.drug_seen_in_training).sum()),
           "response_source_sha256": sha256(root / "data/raw" / FILES["glioma_csv"]),
           "interval": cluster_bootstrap(supported, selected, repetitions=cfg["bootstrap_repetitions"], seed=cfg["seed"]),
           "scope": "Metadatos: fármacos conocidos en líneas reservadas. No caracteriza heterogeneidad molecular.",
           "training_partition": "development_only; test no usado para ajustar ni seleccionar"}
    (folder / "run.json").write_text(json.dumps(run, indent=2, ensure_ascii=False, default=int) + "\n")
    joblib.dump({"model": fitted[selected], "features": features, "run": run}, root / "models/metadata_baseline.joblib")
    return cv, scores, pred, run
