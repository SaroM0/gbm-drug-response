"""Construye la tabla de identidad química de los fármacos ensayados en GBM/glioma.

Entrada: los nombres de fármaco de los adjuntos de respuesta (no se modifican).
Salida: `data/processed/compounds.csv` y un informe de resolución con los casos
ambiguos o sin estructura, que requieren revisión manual antes de usarse.

Uso:
    .venv/bin/python scripts/build_compounds.py
"""
from __future__ import annotations

import json
import sys
import tarfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gbm_drug_response.chemistry import build_table  # noqa: E402
from gbm_drug_response.data import FILES, sha256  # noqa: E402
from gbm_drug_response.runs import tracked, publish_link, write_json

@tracked("chemistry")
def main(root=ROOT, *, _run):
    root = Path(root)
    cache = root / "data/external/pubchem_cache"
    report = _run.path
    publish_link(root / 'data/processed/compounds.csv', root / 'reports/chemistry/compounds.csv',
                 root / 'reports/runs/legacy_processed')
    _run.inputs([root / "data/raw" / FILES[k] for k in ["glioma_csv", "t98g"]])
    responses = pd.read_csv(root / "data" / "raw" / FILES["glioma_csv"])
    t98g = pd.read_csv(root / "data" / "raw" / FILES["t98g"])
    pairs = (responses[["DRUG_ID", "DRUG_NAME", "PUTATIVE_TARGET", "PATHWAY_NAME"]]
             .drop_duplicates().sort_values("DRUG_ID"))
    names = sorted(pairs["DRUG_NAME"].unique())
    print(f"{len(pairs)} DRUG_ID y {len(names)} nombres a resolver.", flush=True)

    chemistry = build_table(names, cache)
    table = pairs.merge(chemistry, left_on="DRUG_NAME", right_on="drug_name", how="left")
    table["assayed_in_t98g"] = table["DRUG_ID"].isin(set(t98g["DRUG_ID"]))
    table = table.rename(columns={"DRUG_ID": "gdsc_drug_id", "PUTATIVE_TARGET": "putative_target",
                                  "PATHWAY_NAME": "pathway_name"}).drop(columns=["DRUG_NAME"])
    order = ["compound_id", "gdsc_drug_id", "drug_name", "canonical_smiles", "inchikey",
             "chemical_form", "conjugation_route", "assayed_in_t98g"]
    table = table[order + [c for c in table.columns if c not in order]]

    table.to_csv(report / 'compounds.csv', index=False)
    cached = sorted(cache.glob('*.json'))
    _run.inputs(cached)
    with tarfile.open(report / 'pubchem_cache.tar.gz', 'w:gz') as archive:
        for path in cached:
            archive.add(path, arcname=path.name)

    resolved = table["canonical_smiles"].notna()
    # Un título distinto puede ser sinónimo válido, pero requiere revisión explícita.
    needs_review = table[(~resolved) | table["ambiguous_name"].fillna(False)
                         | table["query_was_modified"].fillna(False)
                         | ~table["name_matches_title"]
                         | (table["normalization_status"] != "standardized")
                         | (table["parse_status"] != "ok")]
    needs_review.to_csv(report / "resolution_review.csv", index=False)

    summary = {
        "source": "PubChem PUG REST, consulta por nombre",
        "input_file": FILES["glioma_csv"],
        "input_sha256": sha256(root / "data" / "raw" / FILES["glioma_csv"]),
        "drug_ids": int(len(table)),
        "drug_names": int(len(names)),
        "resolved_drug_id_rows": int(resolved.sum()),
        "resolved_names": int(table.loc[resolved, "drug_name"].nunique()),
        "unique_compound_ids": int(table["compound_id"].nunique()),
        "rows_needing_review": int(len(needs_review)),
        "route_counts": table["conjugation_route"].fillna('unresolved').value_counts().to_dict(),
        "note_es": ("La correspondencia nombre->CID es ambigua; las filas del informe de "
                    "revisión no deben usarse sin verificación manual de la identidad."),
    }
    write_json(report / "compounds_run.json", summary)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return summary


if __name__ == "__main__":
    main()
