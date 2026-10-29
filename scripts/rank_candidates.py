"""Prioriza candidatos a conjugar en Fe3O4@SiO2-APTES para ensayarlos en T98G.

No es un predictor de citotoxicidad de la formulación: con ~3 candidatos no existe
conjunto de entrenamiento de formulaciones. Es un embudo de cribado que combina
potencia del fármaco libre medida en GDSC2, factibilidad de dosis en el rango de
nanopartícula que usa el laboratorio y asas químicas presentes en la estructura.

Uso:
    .venv/bin/python scripts/rank_candidates.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gbm_drug_response.data import FILES, load_config, quality_flags, read_responses  # noqa: E402

RAW = ROOT / "data" / "raw"
OUT = ROOT / "reports" / "candidates"

# Condiciones declaradas por el laboratorio (ver docs/plan_investigacion_t98g_nanoparticulas.md).
NP_PPM_MAX = 25.0  # µg/mL de nanopartícula, dosis más alta de la serie 5–25 ppm
# Capacidades publicadas para Fe3O4@SiO2 funcionalizado con APTES: 30–52 µg de fármaco
# por mg de nanopartícula, es decir 3–5,2 % p/p (ver docs/quimica_candidatos.md).
# El límite de 10 % es un techo optimista: el doble de lo mejor reportado.
LOADING_REPORTED_PCT = 5.0
LOADING_FEASIBLE_PCT = 10.0
# Rutas de unión ordenadas por trabajo químico adicional requerido.
ROUTE_TIER = {
    "covalente_cooh_edc_nhs": 1,
    "coordinacion_metal": 2,
    "adsorcion_electrostatica": 2,
    "covalente_amina_requiere_activacion": 3,
    "covalente_requiere_derivatizacion": 4,
    "asa_aromatica_poco_reactiva": 5,
    "sin_asa_reactiva": 5,
}


def required_loading_pct(ic50_um, molecular_weight, np_ppm=NP_PPM_MAX):
    """Carga (% peso fármaco/nanopartícula) necesaria para entregar la IC50 libre a `np_ppm`.

    Supone liberación completa y potencia igual a la del fármaco libre; ignora el efecto
    propio de la nanopartícula y cualquier cambio de captación. Es un orden de magnitud.
    """
    drug_ug_per_ml = ic50_um * molecular_weight / 1000.0
    return drug_ug_per_ml / np_ppm * 100.0


def potency_table(root):
    """Potencia por DRUG_ID en T98G, en el panel GBM/glioma y en todo GDSC2."""
    cfg = load_config(root)
    t98g = read_responses(RAW / FILES["t98g"])
    glioma = read_responses(RAW / FILES["glioma_csv"])
    panel = read_responses(RAW / FILES["gdsc2"])

    for frame in (t98g, glioma, panel):
        flags = quality_flags(frame, cfg["qc_rmse_max"])
        frame["keep"] = flags["keep"]
        frame["above_range"] = flags["ic50_above_range"]
    t98g, glioma, panel = (f[f["keep"]].copy() for f in (t98g, glioma, panel))

    t = (t98g.groupby("DRUG_ID")
         .agg(t98g_ln_ic50=("LN_IC50", "mean"), t98g_auc=("AUC", "mean"),
              t98g_max_conc=("MAX_CONC", "mean"), t98g_above_range=("above_range", "any"))
         .reset_index())
    t["t98g_ic50_um"] = np.exp(t["t98g_ln_ic50"])
    t["t98g_in_range"] = ~t["t98g_above_range"]

    g = (glioma.groupby("DRUG_ID")
         .agg(glioma_lines=("SANGER_MODEL_ID", "nunique"),
              glioma_median_ln_ic50=("LN_IC50", "median"),
              glioma_frac_in_range=("above_range", lambda s: float((~s).mean())))
         .reset_index())

    p = (panel.groupby("DRUG_ID")
         .agg(panel_lines=("SANGER_MODEL_ID", "nunique"),
              panel_median_ln_ic50=("LN_IC50", "median"),
              panel_frac_in_range=("above_range", lambda s: float((~s).mean())))
         .reset_index())

    table = t.merge(g, on="DRUG_ID", how="left").merge(p, on="DRUG_ID", how="left")
    # Positivo: T98G más sensible que la línea mediana del panel completo. Es una
    # señal de selectividad relativa, no una medida de toxicidad en células sanas.
    table["t98g_vs_panel_ln"] = table["panel_median_ln_ic50"] - table["t98g_ln_ic50"]
    return table


def assign_tier(row):
    """Nivel de prioridad. Lexicográfico: cada nivel exige cumplir los requisitos previos."""
    if not row["structure_ok"]:
        return 9, "estructura no resuelta o ambigua"
    if not row["t98g_in_range"]:
        return 8, "IC50 en T98G extrapolada sobre el máximo ensayado"
    if not np.isfinite(row["required_loading_pct"]):
        return 8, "sin peso molecular para evaluar factibilidad de dosis"
    if row["required_loading_pct"] > LOADING_FEASIBLE_PCT:
        return 7, f"requiere carga {row['required_loading_pct']:.1f}% p/p a {NP_PPM_MAX:.0f} ppm"
    tier = ROUTE_TIER.get(row["conjugation_route"], 5)
    return tier, "candidato viable según asa química y dosis alcanzable"


def main():
    compounds = pd.read_csv(ROOT / "data" / "processed" / "compounds.csv")
    compounds["gdsc_drug_id"] = compounds["gdsc_drug_id"].astype(str)
    potency = potency_table(ROOT)
    potency["DRUG_ID"] = potency["DRUG_ID"].astype(str)

    table = compounds.merge(potency, left_on="gdsc_drug_id", right_on="DRUG_ID", how="inner")
    table = table[table["assayed_in_t98g"]].copy()

    table["structure_ok"] = (table["canonical_smiles"].notna()
                             & (table["parse_status"] == "ok")
                             & ~table["ambiguous_name"].fillna(True))
    weight = pd.to_numeric(table["molecular_weight"], errors="coerce")
    table["required_loading_pct"] = required_loading_pct(table["t98g_ic50_um"], weight)
    table["loading_within_reported"] = table["required_loading_pct"] <= LOADING_REPORTED_PCT

    tiers = table.apply(assign_tier, axis=1, result_type="expand")
    table["tier"], table["tier_reason"] = tiers[0], tiers[1]
    table["review_flags"] = table.apply(review_flags, axis=1)

    table = table.sort_values(["tier", "required_loading_pct", "t98g_ln_ic50"],
                              kind="stable").reset_index(drop=True)
    table["rank"] = np.arange(1, len(table) + 1)

    columns = ["rank", "tier", "tier_reason", "drug_name", "gdsc_drug_id", "compound_id",
               "putative_target", "pathway_name", "conjugation_route", "review_flags",
               "t98g_ic50_um", "t98g_ln_ic50", "t98g_in_range", "t98g_max_conc", "t98g_auc",
               "required_loading_pct", "loading_within_reported", "molecular_weight",
               "free_carboxylic_acid", "net_charge_ph74", "hydrolysis_groups",
               "glioma_lines", "glioma_median_ln_ic50", "glioma_frac_in_range",
               "panel_median_ln_ic50", "t98g_vs_panel_ln", "canonical_smiles", "inchikey"]
    ranking = table[columns]

    OUT.mkdir(parents=True, exist_ok=True)
    ranking.to_csv(OUT / "ranking.csv", index=False)
    # Un mismo compuesto puede tener varios DRUG_ID (distintos ensayos de GDSC). La lista
    # corta cuenta moléculas, no identificadores: se conserva la fila mejor clasificada.
    shortlist = ranking[ranking["tier"] <= 4].copy()
    all_ids = (shortlist.groupby("compound_id")["gdsc_drug_id"]
               .apply(lambda s: ";".join(sorted(s.astype(str)))).rename("all_gdsc_drug_ids"))
    shortlist = (shortlist.drop_duplicates("compound_id", keep="first")
                 .merge(all_ids, on="compound_id", how="left"))
    shortlist["rank"] = np.arange(1, len(shortlist) + 1)
    shortlist.to_csv(OUT / "shortlist.csv", index=False)

    summary = {
        "purpose_es": ("Embudo de cribado para elegir candidatos a conjugar; no predice "
                       "la citotoxicidad de la formulación cargada."),
        "np_ppm_max": NP_PPM_MAX,
        "loading_feasible_pct": LOADING_FEASIBLE_PCT,
        "loading_reported_pct": LOADING_REPORTED_PCT,
        "assumptions_es": ("Liberación completa y potencia igual a la del fármaco libre; "
                           "sin efecto propio de la nanopartícula ni cambio de captación."),
        "drugs_assayed_in_t98g": int(len(table)),
        "with_resolved_structure": int(table["structure_ok"].sum()),
        "ic50_in_range_in_t98g": int(table["t98g_in_range"].sum()),
        "dose_feasible_at_max_ppm": int((table["required_loading_pct"] <= LOADING_FEASIBLE_PCT).sum()),
        "shortlist_size": int(len(shortlist)),
        "tier_counts": table["tier"].value_counts().sort_index().to_dict(),
        "route_counts_in_shortlist": shortlist["conjugation_route"].value_counts().to_dict(),
    }
    (OUT / "run.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    print(json.dumps(summary, indent=2, ensure_ascii=False))


def review_flags(row):
    """Banderas que exigen revisión de un químico antes de comprar o conjugar."""
    flags = []
    if row.get("hydrolysis_liability"):
        flags.append("hidrolisis:" + str(row.get("hydrolysis_groups") or ""))
    if row.get("single_handle_review"):
        flags.append("unico_grupo_ionizable")
    # Cuando un fármaco lleva un ácido borónico, casi siempre es su grupo reactivo
    # frente a la diana (p. ej. el proteasoma): derivatizarlo anularía la actividad.
    if row.get("boronic_acid", 0) > 0:
        flags.append("acido_boronico_probable_farmacoforo")
    if row.get("cationic_repulsion_risk"):
        flags.append("cationico_repele_superficie")
    if row.get("high_logp"):
        flags.append("logp_alto")
    if row.get("chemical_form") == "salt_or_mixture":
        flags.append("sal_o_mezcla")
    if row.get("query_was_modified"):
        flags.append("nombre_consultado_modificado")
    return ";".join(flags)


if __name__ == "__main__":
    main()
