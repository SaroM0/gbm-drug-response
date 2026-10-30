"""Normalización descriptiva de pozos públicos; no reproduce el ajuste NLME."""
import numpy as np
import pandas as pd

WELL_KEY = ['BARCODE', 'SCAN_ID', 'POSITION']


def collapse_annotations(frame):
    required = {*WELL_KEY, 'SANGER_MODEL_ID', 'DRUG_ID', 'CONC', 'INTENSITY', 'TAG'}
    if not required <= set(frame) or frame[WELL_KEY + ['SANGER_MODEL_ID']].isna().any().any():
        raise ValueError('Esquema o claves físicas de pozo inválidos.')
    grouped = frame.groupby(WELL_KEY, sort=False, dropna=False)
    # El solvente puede llevar DRUG_ID=4000. Su tag DMSO no es otro fármaco
    # ni debe prevalecer si aparece antes del tag de tratamiento del mismo pozo.
    chemical = frame.copy()
    chemical.loc[chemical.TAG.eq('DMSO'), 'DRUG_ID'] = pd.NA
    chemical.loc[chemical.TAG.eq('DMSO'), 'CONC'] = np.nan
    chemical_grouped = chemical.groupby(WELL_KEY, sort=False, dropna=False)
    stable = ['SANGER_MODEL_ID', 'INTENSITY'] + [k for k in
        ['CELL_ID', 'MASTER_CELL_ID', 'DRUGSET_ID', 'ASSAY', 'DURATION'] if k in frame]
    if grouped[stable].nunique(dropna=True).gt(1).any().any():
        raise ValueError('El mismo pozo/scan tiene metadatos o intensidad contradictorios.')
    wells = grouped.first()
    wells[['DRUG_ID', 'CONC']] = chemical_grouped[['DRUG_ID', 'CONC']].first()
    wells['TAG'] = grouped.TAG.agg(lambda s: ';'.join(sorted(set(s.dropna().astype(str)))))
    wells['annotation_rows'] = grouped.size()
    wells['unsupported_combination'] = chemical_grouped[['DRUG_ID', 'CONC']].nunique().gt(1).any(axis=1)
    wells['unsupported_combination'] |= wells.TAG.str.contains(r'-C(?:;|$)', regex=True)
    return wells.reset_index()


def normalize_controls(wells):
    """Sólo NC-1 y B explícitos. PC-1/UN-USED no se inventan como blancos."""
    out = wells.copy()
    tags = out.TAG.str.split(';')
    out['is_vehicle'] = tags.map(lambda x: 'NC-1' in x)
    out['is_blank'] = tags.map(lambda x: 'B' in x)
    out['is_library_single'] = out.TAG.str.contains(r'(?:^|;)L\d+-D\d+-S(?:;|$)', regex=True)
    out['normalized_viability'] = np.nan
    out['normalization_status'] = 'missing_explicit_controls'
    summaries = []
    for (barcode, scan), part in out.groupby(['BARCODE', 'SCAN_ID'], sort=False):
        vehicle = part.loc[part.is_vehicle & ~part.unsupported_combination, 'INTENSITY']
        blank = part.loc[part.is_blank & ~part.unsupported_combination, 'INTENSITY']
        valid = len(vehicle) >= 2 and len(blank) >= 2
        v, b = vehicle.mean(), blank.mean()
        if not np.isfinite([v, b]).all() or v <= b:
            valid = False
        status = 'normalized_explicit_nc1_b' if valid else 'missing_or_invalid_explicit_controls'
        if valid:
            out.loc[part.index, 'normalized_viability'] = (part.INTENSITY-b)/(v-b)
        out.loc[part.index, 'normalization_status'] = status
        summaries.append({'BARCODE': barcode, 'SCAN_ID': scan, 'n_vehicle_wells': len(vehicle),
            'n_blank_wells': len(blank), 'vehicle_mean': v, 'blank_mean': b, 'status': status,
            'z_prime_descriptive': 1-3*(vehicle.std(ddof=1)+blank.std(ddof=1))/(v-b) if valid else np.nan})
    # Nunca truncar ni considerar un Z' descriptivo como aprobación de QC original.
    return out, pd.DataFrame(summaries)


def crossing_diagnostics(wells):
    eligible = wells[wells.is_library_single & ~wells.unsupported_combination & wells.DRUG_ID.notna()
        & np.isfinite(wells.CONC) & wells.CONC.gt(0) & np.isfinite(wells.normalized_viability)]
    rows = []
    keys = ['SANGER_MODEL_ID', 'DRUG_ID', 'BARCODE', 'SCAN_ID']
    for identity, part in eligible.groupby(keys, sort=False):
        doses = part.groupby('CONC').normalized_viability.median().sort_index()
        above, below = (doses > .5).any(), (doses <= .5).any()
        status = 'points_on_both_sides_of_half' if above and below else (
            'all_points_above_half' if above else 'all_points_at_or_below_half')
        if len(doses) < 2:
            status = 'insufficient_dose_levels'
        rows.append(dict(zip(keys, identity)) | {'n_wells': len(part), 'n_doses': len(doses),
            'raw_min_conc': doses.index.min(), 'raw_max_conc': doses.index.max(),
            'median_viability_at_max': float(doses.iloc[-1]),
            'empirical_range_state': status, 'increasing_adjacent_medians': int((np.diff(doses)>0).sum()),
            'nlme_refit_performed': False})
    return pd.DataFrame(rows, columns=keys+['n_wells', 'n_doses', 'raw_min_conc', 'raw_max_conc',
        'median_viability_at_max', 'empirical_range_state', 'increasing_adjacent_medians', 'nlme_refit_performed'])
