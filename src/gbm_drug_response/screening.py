"""Cribado exploratorio por entidad activa, sin seleccionar su mejor ensayo."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .data import FILES, load_config, quality_flags, read_responses, range_state

ROUTE_TIER = {
    'covalente_cooh_edc_nhs': 1, 'coordinacion_metal': 2,
    'adsorcion_electrostatica': 2, 'covalente_amina_requiere_activacion': 3,
    'covalente_requiere_derivatizacion': 4,
    'asa_aromatica_poco_reactiva': 5, 'sin_asa_reactiva': 5,
}


def load_screening(root):
    cfg = json.loads((Path(root) / 'screening.json').read_text())
    for key in ['np_doses_ug_ml', 'loading_scenarios_pct']:
        values = np.asarray(cfg[key], float)
        if not len(values) or not np.isfinite(values).all() or (values <= 0).any():
            raise ValueError(f'{key}: se requieren valores positivos finitos.')
        if len(set(values)) != len(values):
            raise ValueError(f'{key}: valores duplicados.')
    fractions = np.asarray(cfg['released_fraction_scenarios'], float)
    if not len(fractions) or not np.isfinite(fractions).all() or ((fractions <= 0) | (fractions > 1)).any():
        raise ValueError('Fracciones de liberación deben estar en (0, 1].')
    if cfg['reference_loading_pct'] not in cfg['loading_scenarios_pct'] or cfg['reference_released_fraction'] not in fractions:
        raise ValueError('El escenario de referencia debe formar parte de la sensibilidad.')
    if not np.isfinite(cfg['screening_loading_limit_pct']) or cfg['screening_loading_limit_pct'] <= 0:
        raise ValueError('Límite de cribado inválido.')
    if cfg['potency_summary'] != 'maximum_ln_ic50_across_all_assays' or cfg['range_policy'] != 'all_assays_within_range':
        raise ValueError('Política de potencia/rango no implementada.')
    if cfg['mass_basis'] != 'active_parent_drug_mass_per_carrier_mass':
        raise ValueError('Base de masa no implementada.')
    return cfg


def required_loading_pct(ic50_um, molecular_weight, np_ppm=25.0, released_fraction=1.0):
    """Escenario de masa, no capacidad de carga medida ni potencia nano predicha."""
    if not np.isfinite(np_ppm) or np_ppm <= 0 or not 0 < released_fraction <= 1:
        raise ValueError('Dosis positiva y fracción de liberación en (0, 1] requeridas.')
    concentration, weight = np.asarray(ic50_um, float), np.asarray(molecular_weight, float)
    valid = np.isfinite(concentration) & (concentration > 0) & np.isfinite(weight) & (weight > 0)
    return np.where(valid, concentration * weight / (10 * np_ppm * released_fraction), np.nan)


def potency_table(root, frames=None):
    """Una fila por ensayo; las fracciones de rango verifican ambos extremos."""
    root = Path(root)
    cfg = load_config(root)
    if frames is None:
        frames = {k: read_responses(root / 'data/raw' / FILES[k]) for k in ['t98g', 'glioma_csv', 'gdsc2']}
    clean = {}
    for key, frame in frames.items():
        f = quality_flags(frame, cfg['qc_rmse_max'])
        df = frame.loc[f.keep].copy()
        df['range_state'] = range_state(df)
        clean[key] = df
    t = clean['t98g'].rename(columns={
        'DRUG_ID': 'gdsc_drug_id', 'LN_IC50': 't98g_ln_ic50', 'AUC': 't98g_auc',
        'MIN_CONC': 't98g_min_conc', 'MAX_CONC': 't98g_max_conc',
        'range_state': 't98g_range_state', 'RMSE': 'curve_rmse'})
    if t.gdsc_drug_id.duplicated().any():
        raise ValueError('Más de una curva por DRUG_ID en T98G; resolver antes del cribado.')
    columns = ['gdsc_drug_id', 'NLME_CURVE_ID', 't98g_ln_ic50', 't98g_auc',
               't98g_min_conc', 't98g_max_conc', 't98g_range_state', 'curve_rmse']
    table = t[columns].copy()
    table['t98g_ic50_um'] = np.exp(table.t98g_ln_ic50)
    table['t98g_in_range'] = table.t98g_range_state.eq('within_range')
    for key, prefix in [('glioma_csv', 'glioma'), ('gdsc2', 'panel')]:
        summary = clean[key].groupby('DRUG_ID').agg(
            lines=('SANGER_MODEL_ID', 'nunique'), median_ln_ic50=('LN_IC50', 'median'),
            frac_in_range=('range_state', lambda s: float(s.eq('within_range').mean())))
        summary.columns = [prefix + '_' + x for x in summary.columns]
        table = table.merge(summary, left_on='gdsc_drug_id', right_index=True, how='left', validate='one_to_one')
    table['t98g_vs_panel_ln'] = table.panel_median_ln_ic50 - table.t98g_ln_ic50
    return table


def review_flags(row):
    flags = ['identity_and_route_require_manual_review', 'release_and_activity_unmeasured']
    if row.get('hydrolysis_liability') is True or row.get('hydrolysis_liability') == True:
        flags.append('hydrolysis_liability')
    for column, label in [('single_carboxyl_review', 'single_carboxyl_requires_review'),
                          ('electrostatic_review', 'pka_zeta_medium_require_review'),
                          ('high_logp', 'high_logp'), ('query_was_modified', 'modified_name_query')]:
        if pd.notna(row.get(column)) and bool(row.get(column)):
            flags.append(label)
    if pd.notna(row.get('boronic_acid')) and row.get('boronic_acid', 0) > 0:
        flags.append('boronic_acid_pharmacophore_requires_review')
    if row.get('chemical_form') == 'salt_or_mixture':
        flags.append('purchased_form_mass_basis_requires_review')
    return ';'.join(flags)


def assign_tier(row, cfg):
    if not row['structure_ok']:
        return 9, 'estructura ausente/ambigua o mezcla pendiente de resolver'
    if not row['t98g_in_range']:
        return 8, 'uno o más ensayos fuera del rango; no resumir eligiendo el favorable'
    if not np.isfinite(row['required_loading_pct']):
        return 8, 'base de masa o potencia no evaluable'
    if row['required_loading_pct'] > cfg['screening_loading_limit_pct']:
        return 7, 'supera el escenario de carga elegido; no demuestra imposibilidad física'
    return ROUTE_TIER.get(row['conjugation_route'], 5), 'prioridad exploratoria; química y dosis experimentales no verificadas'


def summarize_compounds(assays, cfg):
    """Máximo conservador en escala ln entre TODOS los ensayos de una entidad."""
    assays = assays.copy()
    assays['_identity'] = assays.compound_id.fillna('unresolved:' + assays.gdsc_drug_id.astype(str))
    rows = []
    for identity, part in assays.groupby('_identity', sort=True):
        # Referencia elegida por ID, nunca por potencia ni por prioridad favorable.
        part = part.assign(_id=part.gdsc_drug_id.astype(int)).sort_values('_id')
        row = part.iloc[0].drop(labels=['_identity', '_id']).to_dict()
        row['drug_name'] = ';'.join(sorted(part.drug_name.unique()))
        row['all_gdsc_drug_ids'] = ';'.join(part.gdsc_drug_id.astype(str))
        row['n_assays'] = len(part)
        row['t98g_ln_ic50'] = float(part.t98g_ln_ic50.max())
        row['t98g_ic50_um'] = float(np.exp(row['t98g_ln_ic50']))
        row['ic50_min_um'] = float(part.t98g_ic50_um.min())
        row['ic50_max_um'] = float(part.t98g_ic50_um.max())
        row['assay_disagreement_fold'] = float(np.exp(part.t98g_ln_ic50.max() - part.t98g_ln_ic50.min()))
        row['assay_range_states'] = ';'.join(sorted(part.t98g_range_state.unique()))
        row['t98g_in_range'] = bool(part.t98g_in_range.all())
        row['structure_ok'] = bool(part.structure_ok.all()) and part.conjugation_route.nunique() == 1
        row['required_loading_pct'] = float(required_loading_pct(row['t98g_ic50_um'], row['active_molecular_weight'],
            max(cfg['np_doses_ug_ml']), cfg['reference_released_fraction']))
        row['tier'], row['tier_reason'] = assign_tier(row, cfg)
        row['review_flags'] = ';'.join(sorted({f for x in part.review_flags for f in x.split(';')}))
        # No aparentar que AUC/QC/rango de un ensayo describe la molécula completa.
        for field in ['NLME_CURVE_ID', 'curve_rmse', 't98g_auc', 't98g_min_conc',
                      't98g_max_conc', 't98g_range_state', 'gdsc_drug_id', 't98g_vs_panel_ln']:
            row.pop(field, None)
        rows.append(row)
    table = pd.DataFrame(rows).sort_values(['tier', 'required_loading_pct', 'drug_name'], kind='stable').reset_index(drop=True)
    table.insert(0, 'rank', np.arange(1, len(table) + 1))
    return table


def scenario_table(ranking, cfg):
    rows = []
    for r in ranking.itertuples():
        if not r.structure_ok or not np.isfinite(r.active_molecular_weight) or r.active_molecular_weight <= 0:
            continue
        for load in cfg['loading_scenarios_pct']:
            for release in cfg['released_fraction_scenarios']:
                drug_doses = np.asarray(cfg['np_doses_ug_ml']) * load / 100 * release * 1000 / r.active_molecular_weight
                low, high = drug_doses.min(), drug_doses.max()
                status = 'brackets_reference' if low <= r.t98g_ic50_um <= high else (
                    'all_doses_above_reference' if low > r.t98g_ic50_um else 'all_doses_below_reference')
                required = float(required_loading_pct(r.t98g_ic50_um, r.active_molecular_weight,
                                                       max(cfg['np_doses_ug_ml']), release))
                rows.append({'compound_id': r.compound_id, 'drug_name': r.drug_name,
                             'loading_pct': load, 'released_fraction': release,
                             'equivalent_drug_min_um': low, 'equivalent_drug_max_um': high,
                             'required_loading_pct': required,
                             'passes_loading_scenario': required <= load,
                             'bracketing_status': status,
                             'potency_reference_in_range': r.t98g_in_range,
                             'experimental_prediction': False})
    return pd.DataFrame(rows)


def cumulative_funnel(assays):
    selected = assays.copy()
    rows = []
    stages = [('t98g_qc_assays', None), ('structure_accepted', 'structure_ok'),
              ('all_assays_in_range', 'compound_all_in_range'),
              ('within_chosen_loading_scenario', 'compound_loading_pass'),
              ('structural_route_tier_1_to_4', 'compound_route_pass')]
    for stage, column in stages:
        if column:
            selected = selected[selected[column]]
        rows.append({'stage': stage, 'drug_id_rows': len(selected),
                     'distinct_nonmissing_compound_ids': selected.compound_id.nunique()})
    return pd.DataFrame(rows)
