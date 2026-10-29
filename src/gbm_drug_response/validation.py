"""Dependencia biológica y diagnósticos de ajustes fuera del rango experimental."""
from pathlib import Path
import numpy as np
import pandas as pd
from .data import range_state


def annotate_validation_groups(root, df, folder):
    """Agrupa por relaciones oficiales locales; no rehace el holdout observado."""
    path = Path(root) / 'data/external/model_list_20260921.csv'
    out = df.copy()
    if not path.exists():
        out['VALIDATION_GROUP'] = out.SANGER_MODEL_ID
        return out, {'status': 'unavailable', 'donor_independence_verified': False}
    meta = pd.read_csv(path, dtype='string')
    required = {'model_id', 'patient_id', 'parent_id'}
    if not required <= set(meta) or meta.model_id.duplicated().any():
        raise ValueError('Tabla de relaciones incompleta o modelos duplicados.')
    parent = {}
    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a, b):
        left, right = sorted([find(a), find(b)])
        parent[right] = left
    patients = {}
    for row in meta.itertuples():
        find(row.model_id)
        if pd.notna(row.parent_id):
            union(row.model_id, row.parent_id)
        if pd.notna(row.patient_id):
            if row.patient_id in patients:
                union(row.model_id, patients[row.patient_id])
            patients[row.patient_id] = row.model_id
    ids = out[['SANGER_MODEL_ID']].drop_duplicates()
    joined = ids.merge(meta[['model_id', 'patient_id', 'parent_id']], left_on='SANGER_MODEL_ID',
                       right_on='model_id', how='left', validate='one_to_one')
    joined['VALIDATION_GROUP'] = joined.SANGER_MODEL_ID.map(lambda x: find(x))
    joined['VALIDATION_GROUP'] = joined.groupby('VALIDATION_GROUP').SANGER_MODEL_ID.transform('min')
    joined.to_csv(Path(folder) / 'biological_groups.csv', index=False)
    out['VALIDATION_GROUP'] = out.SANGER_MODEL_ID.map(joined.set_index('SANGER_MODEL_ID').VALIDATION_GROUP)
    return out, {'status': 'local_metadata_checked', 'matched_lines': int(joined.model_id.notna().sum()),
                 'distinct_patient_ids': int(joined.patient_id.nunique()),
                 'donor_independence_verified': bool(joined.patient_id.notna().all()
                     and joined.patient_id.nunique() == len(joined)),
                 'independent_groups': int(joined.VALIDATION_GROUP.nunique()),
                 'official_download_authenticity_verified': False}


def outcome_diagnostics(predictions, model_names, metric_function, folder):
    """Sensibilidad descriptiva con los mismos modelos; no selecciona por test."""
    state = range_state(predictions)
    rows = []
    for group in ['all_fitted_estimates', 'within_range', 'below_range', 'above_range']:
        part = predictions if group == 'all_fitted_estimates' else predictions.loc[state.eq(group)]
        for name in model_names:
            row = {'range_group': group, 'model': name, 'n_lines': part.SANGER_MODEL_ID.nunique()}
            row.update(metric_function(part.LN_IC50, part[name]) if len(part) else
                       {'n': 0, 'mae': np.nan, 'prediction_rmse': np.nan, 'r2': np.nan, 'spearman': np.nan})
            rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(Path(folder) / 'outcome_diagnostics.csv', index=False)
    return {'target_interpretation': 'fitted_point_estimate_not_measured_censoring_bound',
            'evaluation_status': 'exploratory_previously_observed_holdout',
            'selection_uses_test_subgroups': False,
            'n_within_range': int(state.eq('within_range').sum()),
            'n_above_range': int(state.eq('above_range').sum()),
            'n_below_range': int(state.eq('below_range').sum()),
            'validated_formulation_or_clinical_prediction': False}
