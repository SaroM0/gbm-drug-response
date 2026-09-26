"""Pruebas del protocolo; datos sintéticos sólo aquí, nunca en resultados científicos."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from gbm_drug_response.data import SCHEMA, FILES, quality_flags, read_responses, sha256
from gbm_drug_response.modeling import DrugMeanRegressor, get_split, attach_split, grouped_folds, cluster_bootstrap
from gbm_drug_response.molecular import load_expression, train_molecular, molecular_candidates

PROJECT = Path(__file__).resolve().parents[1]


def responses(n=40):
    rows = []
    for i in range(n):
        rows.append(dict(zip(SCHEMA, [
            "GDSC2", 343, i + 1, f"test_line_{i}", f"SIDM{i:05d}", "Glioblastoma",
            1003, "synthetic_drug", "TEST", "Test pathway", 0.001, 10,
            i / n, 0.5, 0.1, 0.0,
        ])))
    return pd.DataFrame(rows)


@pytest.fixture
def experiment(tmp_path):
    for directory in ["data/raw", "data/processed", "data/external", "reports", "models"]:
        (tmp_path / directory).mkdir(parents=True)
    cfg = json.loads((PROJECT / "config.json").read_text())
    cfg.update(molecular_max_drugs=1, bootstrap_repetitions=40)
    (tmp_path / "config.json").write_text(json.dumps(cfg))
    responses().to_csv(tmp_path / "data/raw" / FILES["glioma_csv"], index=False)
    return tmp_path, cfg


def test_quality_boundaries_and_extrapolation():
    df = responses(6)
    df.loc[0, "RMSE"] = 0.3
    df.loc[1, "RMSE"] = 0.300001
    df.loc[2, "LN_IC50"] = np.nan
    df.loc[3, "MIN_CONC"] = 0
    df.loc[4, "AUC"] = 1.1
    df.loc[5, "LN_IC50"] = 30  # Extrapolación se marca, no se elimina.
    f = quality_flags(df)
    assert f.keep.tolist() == [True, False, False, False, False, True]
    assert f.loc[5, "ic50_above_range"]


def test_missing_target_annotation_is_not_missing_outcome(tmp_path):
    path = tmp_path / "table.csv"
    df = responses(3)
    df["PUTATIVE_TARGET"] = ["", " ", "TOP1"]
    df.to_csv(path, index=False)
    read = read_responses(path)
    assert read.PUTATIVE_TARGET.isna().sum() == 2
    assert read.DRUG_ID.tolist() == ["1003"] * 3
    assert quality_flags(read).keep.all()


def test_group_split_persists_and_blocks_tampering(experiment):
    root, cfg = experiment
    split = get_split(root, cfg)
    pd.testing.assert_frame_equal(split, get_split(root, cfg))
    df = read_responses(root / "data/raw" / FILES["glioma_csv"])
    df = pd.concat([df, df.assign(DRUG_ID="1004")], ignore_index=True)
    dev, test = attach_split(df, split)
    assert set(dev.SANGER_MODEL_ID).isdisjoint(test.SANGER_MODEL_ID)
    for train, valid in grouped_folds(dev, 4):
        assert set(dev.iloc[train].SANGER_MODEL_ID).isdisjoint(dev.iloc[valid].SANGER_MODEL_ID)
    split.loc[0, "partition"] = "test" if split.loc[0, "partition"] == "development" else "development"
    split.to_csv(root / "data/processed/line_split.csv", index=False)
    with pytest.raises(ValueError, match="asignación"):
        get_split(root, cfg)


def test_drug_mean_uses_training_only_and_handles_unknown():
    model = DrugMeanRegressor().fit(pd.DataFrame({"DRUG_ID": ["a", "a", "b"]}), [1, 3, 20])
    pred = model.predict(pd.DataFrame({"DRUG_ID": ["a", "b", "new"]}))
    np.testing.assert_allclose(pred, [2, 20, 8])


def test_molecular_preprocessing_does_not_refit_on_test():
    cfg = json.loads((PROJECT / "config.json").read_text())
    rng = np.random.default_rng(4)
    X = pd.DataFrame(rng.uniform(0, 10, size=(30, 20)))
    X.iloc[0, 1] = np.nan
    model = molecular_candidates(cfg)["pca_elastic_net"][0]
    model.fit(X.iloc[:24], rng.normal(size=24))
    original_imputer = model.named_steps["impute"].statistics_.copy()
    original_components = model.named_steps["pca"].pca_.components_.copy()
    extreme_test = X.iloc[24:].fillna(0) + 100000
    assert np.isfinite(model.predict(extreme_test)).all()
    np.testing.assert_array_equal(model.named_steps["impute"].statistics_, original_imputer)
    np.testing.assert_array_equal(model.named_steps["pca"].pca_.components_, original_components)
    np.testing.assert_allclose(original_imputer, np.nanmedian(X.iloc[:24], axis=0))


def write_synthetic_expression(root, n=40):
    rng = np.random.default_rng(42)
    expr = pd.DataFrame(rng.uniform(0, 10, size=(n, 25)), columns=[f"gene_test{i}" for i in range(25)])
    expr.insert(0, "SANGER_MODEL_ID", [f"SIDM{i:05d}" for i in range(n)])
    path = root / "data/external/expression.csv"
    expr.to_csv(path, index=False)
    meta = {"source_url": "synthetic-unit-test-only", "release": "test", "scale": "log2_tpm_plus1",
            "gene_id_type": "synthetic", "preprocessing": "Sólo para verificar código en directorio temporal",
            "sha256": sha256(path)}
    (root / "data/external/expression_metadata.json").write_text(json.dumps(meta))
    return expr, path, meta


def test_missing_expression_does_not_fabricate_results(experiment):
    root, _ = experiment
    status, scores, pred = train_molecular(root)
    assert status["status"] == "missing_expression"
    assert scores.empty and pred.empty
    assert not (root / "models/molecular.joblib").exists()


def test_expression_checksum_and_duplicates(experiment):
    root, _ = experiment
    expr, path, meta = write_synthetic_expression(root)
    assert load_expression(root)[0].shape == expr.shape
    expr.iloc[[0, 0]].to_csv(path, index=False)
    with pytest.raises(ValueError, match="checksum"):
        load_expression(root)
    meta["sha256"] = sha256(path)
    (root / "data/external/expression_metadata.json").write_text(json.dumps(meta))
    with pytest.raises(ValueError, match="única fila"):
        load_expression(root)


def test_molecular_end_to_end_on_synthetic_fixture(experiment):
    root, cfg = experiment
    write_synthetic_expression(root)
    status, scores, pred = train_molecular(root)
    assert status["status"] == "trained"
    assert set(scores.model) == {"drug_mean", "pca_elastic_net", "pca_extra_trees"}
    split = get_split(root, cfg)
    assert set(pred.SANGER_MODEL_ID) <= set(split.loc[split.partition == "test", "SANGER_MODEL_ID"])
    assert len(pred) == 8
    assert (root / "models/molecular.joblib").exists()
    assert scores.selected_by_cv.sum() == 1


def test_bootstrap_uses_line_means_not_independent_curves():
    # La línea con más curvas no debe dominar el promedio macro de líneas.
    p = pd.DataFrame({"SANGER_MODEL_ID": ["a"] * 100 + ["b"], "LN_IC50": [0.] * 101,
                      "candidate": [1.] * 100 + [9.], "drug_mean": [2.] * 101})
    report = cluster_bootstrap(p, "candidate", repetitions=50)
    assert report["line_macro_mae"] == 5.0
    assert report["delta_vs_baseline"] == 3.0
