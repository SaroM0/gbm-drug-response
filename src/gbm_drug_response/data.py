"""Lectura, auditoría y control de calidad. Nunca modifica los originales."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

SOURCE = "https://depmap.sanger.ac.uk/documentation/datasets/drug-sensitivity/"
FILES = {
    "gdsc2": "GDSC2_fitted_dose_response_27Oct23.xlsx",
    "glioma_csv": "GBMGliomafiltrado.csv",
    "glioma_xlsx": "GBMGliomafiltrado.xlsx",
    "t98g": "Aislados T98G - Hoja 1.csv",
}
# nombre: (definición, rol, unidad)
SCHEMA = {
    "DATASET": ("Colección de origen", "procedencia", "texto"),
    "NLME_RESULT_ID": ("Identificador del conjunto de ajustes", "identificador", "ID"),
    "NLME_CURVE_ID": ("Identificador de curva ajustada", "clave de observación", "ID"),
    "CELL_LINE_NAME": ("Nombre de línea celular", "identificador", "texto"),
    "SANGER_MODEL_ID": ("Identificador Sanger del modelo celular", "grupo de validación / unión", "ID"),
    "CANCER_TYPE": ("Tipo de cáncer anotado", "covariable opcional", "categoría"),
    "DRUG_ID": ("Identificador del fármaco ensayado", "unión / predictor para fármacos conocidos", "ID"),
    "DRUG_NAME": ("Nombre del fármaco", "anotación", "texto"),
    "PUTATIVE_TARGET": ("Diana farmacológica propuesta", "anotación del fármaco", "texto"),
    "PATHWAY_NAME": ("Vía de la diana farmacológica", "anotación; no actividad de vías celulares", "categoría"),
    "MIN_CONC": ("Concentración mínima ensayada", "control del ensayo", "µM"),
    "MAX_CONC": ("Concentración máxima ensayada", "control del ensayo", "µM"),
    "LN_IC50": ("Logaritmo natural de IC50 ajustada", "outcome primario", "ln(IC50 / µM)"),
    "AUC": ("Área normalizada bajo curva de respuesta", "outcome alternativo; excluir de X", "fracción"),
    "RMSE": ("Error del ajuste dosis–respuesta", "control de calidad; excluir de X", "escala de viabilidad normalizada"),
    "Z_SCORE": ("LN_IC50 estandarizado por fármaco", "derivado del outcome; excluir de X", "adimensional"),
}
NUMERIC = ["MIN_CONC", "MAX_CONC", "LN_IC50", "AUC", "RMSE", "Z_SCORE"]
INTEGER_IDS = ["NLME_RESULT_ID", "NLME_CURVE_ID", "DRUG_ID"]
KEY = ["SANGER_MODEL_ID", "DRUG_ID"]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_config(root):
    cfg = json.loads((Path(root) / "config.json").read_text())
    if cfg["target"] != "LN_IC50" or cfg["group"] != "SANGER_MODEL_ID":
        raise ValueError("Este protocolo implementa LN_IC50 y grupos SANGER_MODEL_ID.")
    if not 0 < cfg["qc_rmse_max"] <= 0.3:
        raise ValueError("RMSE máximo debe estar en (0, 0.3].")
    if cfg["metadata_features"] != ["DRUG_ID", "CANCER_TYPE"]:
        raise ValueError("Revisar explícitamente el protocolo antes de cambiar predictores.")
    return cfg


def read_responses(path):
    path = Path(path)
    if path.suffix == ".xlsx":
        df = pd.read_excel(path, engine="openpyxl", keep_default_na=False)
    else:
        df = pd.read_csv(path, keep_default_na=False)
    if list(df.columns) != list(SCHEMA):
        raise ValueError(f"Esquema inesperado en {path.name}: {list(df.columns)}")
    for col in df:
        if col in NUMERIC:
            # errors=raise: texto inesperado nunca desaparece silenciosamente.
            df[col] = pd.to_numeric(df[col].replace(r"^\s*$", np.nan, regex=True), errors="raise").astype(float)
        elif col in INTEGER_IDS:
            df[col] = pd.to_numeric(df[col].replace(r"^\s*$", np.nan, regex=True), errors="raise").astype("Int64").astype("string")
        else:
            df[col] = df[col].astype("string").str.strip().replace("", pd.NA)
    return df


def data_dictionary(df):
    rows = []
    for col, (definition, role, unit) in SCHEMA.items():
        rows.append({"column": col, "definition_es": definition, "role": role,
                     "unit": unit, "dtype": str(df[col].dtype),
                     "n_missing": int(df[col].isna().sum()),
                     "pct_missing": float(df[col].isna().mean() * 100),
                     "n_unique_nonmissing": int(df[col].nunique()),
                     "examples": " | ".join(df[col].dropna().astype(str).unique()[:3]),
                     "definition_source": SOURCE})
    return pd.DataFrame(rows)


def quality_flags(df, rmse_max=0.3):
    f = pd.DataFrame(index=df.index)
    f["invalid_target"] = ~np.isfinite(df["LN_IC50"])
    f["invalid_rmse"] = ~np.isfinite(df["RMSE"]) | ~df["RMSE"].between(0, rmse_max)
    f["invalid_auc"] = ~np.isfinite(df["AUC"]) | ~df["AUC"].between(0, 1)
    f["invalid_concentrations"] = (
        ~np.isfinite(df["MIN_CONC"]) | ~np.isfinite(df["MAX_CONC"])
        | (df["MIN_CONC"] <= 0) | (df["MAX_CONC"] <= df["MIN_CONC"])
    )
    f["missing_key"] = df[["SANGER_MODEL_ID", "DRUG_ID", "NLME_CURVE_ID", "CANCER_TYPE"]].isna().any(axis=1)
    f["keep"] = ~f.any(axis=1)
    # Son banderas de extrapolación, no reglas para eliminar o truncar y.
    f["ic50_below_range"] = df["LN_IC50"] < np.log(df["MIN_CONC"].where(df["MIN_CONC"] > 0))
    f["ic50_above_range"] = df["LN_IC50"] > np.log(df["MAX_CONC"].where(df["MAX_CONC"] > 0))
    return f


def assert_subset(child, parent):
    if child["NLME_CURVE_ID"].duplicated().any() or parent["NLME_CURVE_ID"].duplicated().any():
        raise ValueError("NLME_CURVE_ID repetido: investigar antes de unir.")
    a = child.set_index("NLME_CURVE_ID").sort_index()
    b = parent.set_index("NLME_CURVE_ID").reindex(a.index)
    pd.testing.assert_frame_equal(a, b, check_dtype=False, check_exact=False, rtol=0, atol=1e-12)


def audit(root):
    root = Path(root)
    cfg = load_config(root)
    report = root / "reports" / "audit"
    report.mkdir(parents=True, exist_ok=True)
    frames, inventory, manifests = {}, [], []
    for label, name in FILES.items():
        path = root / "data/raw" / name
        df = read_responses(path)
        frames[label] = df
        inventory.append({"dataset": label, "file": name, "rows": len(df), "columns": len(df.columns),
                          "cell_lines": df.SANGER_MODEL_ID.nunique(), "drug_ids": df.DRUG_ID.nunique(),
                          "drug_names": df.DRUG_NAME.nunique(), "duplicate_rows": int(df.duplicated().sum()),
                          "duplicate_curve_ids": int(df.NLME_CURVE_ID.duplicated().sum()),
                          "duplicate_model_drug": int(df.duplicated(KEY).sum())})
        manifests.append({"file": name, "bytes": path.stat().st_size, "sha256": sha256(path),
                          "origin": "Adjunto aportado por el usuario; identidad con distribución oficial no verificada"})
        data_dictionary(df).to_csv(report / f"dictionary_{label}.csv", index=False)
    assert_subset(frames["glioma_csv"], frames["glioma_xlsx"])
    if len(frames["glioma_csv"]) != len(frames["glioma_xlsx"]):
        raise ValueError("CSV y XLSX difieren en número de filas.")
    assert_subset(frames["glioma_csv"], frames["gdsc2"])
    assert_subset(frames["t98g"], frames["glioma_csv"])
    for label, df in frames.items():
        if df.duplicated(KEY).any() or df.NLME_CURVE_ID.duplicated().any():
            raise ValueError(f"Duplicados en {label}; se requiere resolverlos explícitamente.")
        if df.groupby("SANGER_MODEL_ID").CELL_LINE_NAME.nunique().max() > 1:
            raise ValueError("Un SIDM tiene varios nombres: revisar mapeo.")
    df = frames["glioma_csv"]
    flags = quality_flags(df, cfg["qc_rmse_max"])
    pd.concat([df, flags], axis=1).to_parquet(root / "data/processed/glioma_flagged.parquet", index=False)
    df.loc[flags.keep].to_parquet(root / "data/processed/glioma_qc.parquet", index=False)
    flags.sum().rename("rows").to_csv(report / "quality_counts.csv")
    sensitivity = []
    for threshold in cfg["qc_sensitivity_thresholds"]:
        keep = quality_flags(df, threshold).keep
        sensitivity.append({"rmse_max": threshold, "retained_rows": int(keep.sum()),
                            "removed_rows": int((~keep).sum()),
                            "retained_lines": df.loc[keep, "SANGER_MODEL_ID"].nunique()})
    pd.DataFrame(sensitivity).to_csv(report / "rmse_sensitivity.csv", index=False)
    (report / "manifest.json").write_text(json.dumps(manifests, indent=2, ensure_ascii=False) + "\n")
    (report / "relationships.json").write_text(json.dumps({
        "glioma_csv_equals_xlsx": True, "glioma_subset_gdsc2": True, "t98g_subset_glioma": True,
        "float_tolerance": 1e-12, "originals_changed": False,
    }, indent=2) + "\n")
    inventory = pd.DataFrame(inventory)
    inventory.to_csv(report / "inventory.csv", index=False)
    df.groupby(["CANCER_TYPE", "SANGER_MODEL_ID", "CELL_LINE_NAME"]).size().rename("curves").to_csv(report / "cell_lines.csv")
    coverage = df.groupby(["DRUG_ID", "DRUG_NAME"]).SANGER_MODEL_ID.nunique().rename("n_lines").reset_index()
    coverage.to_csv(report / "drug_coverage.csv", index=False)
    aliases = df[["DRUG_NAME", "DRUG_ID", "MIN_CONC", "MAX_CONC"]].drop_duplicates()
    aliases[aliases.DRUG_NAME.isin(df.groupby("DRUG_NAME").DRUG_ID.nunique().loc[lambda x: x > 1].index)].to_csv(report / "multiple_drug_ids.csv", index=False)
    return frames, inventory, pd.DataFrame(sensitivity)
