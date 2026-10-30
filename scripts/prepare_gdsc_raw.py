"""Extrae GBM/glioma del ZIP crudo y comprueba vínculo con respuestas aportadas."""
from pathlib import Path
import json
import sys
import zipfile
import argparse

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from gbm_drug_response.data import FILES, read_responses
from gbm_drug_response.runs import ResearchRun, checksum, write_json, verify_run
from gbm_drug_response.gdsc_raw import collapse_annotations, normalize_controls, crossing_diagnostics


def compare_fitted(user, official):
    keys = ['NLME_CURVE_ID', 'SANGER_MODEL_ID', 'DRUG_ID']
    fields = ['LN_IC50', 'MIN_CONC', 'MAX_CONC', 'AUC', 'RMSE']
    if not set(keys+fields) <= set(official):
        raise ValueError('Fuente oficial sin claves/outcomes esperados.')
    for df in [user, official]:
        for col in ['NLME_CURVE_ID', 'DRUG_ID']:
            df[col] = pd.to_numeric(df[col], errors='raise').astype('Int64').astype('string')
    if official.duplicated(keys).any() or user.duplicated(keys).any():
        raise ValueError('Claves duplicadas al comparar ajustes.')
    merged = user[keys+fields].merge(official[keys+fields], on=keys, how='left',
        suffixes=('_user', '_official'), indicator=True, validate='one_to_one')
    same = merged._merge.eq('both')
    for field in fields:
        same &= np.isclose(merged[field+'_user'], merged[field+'_official'], rtol=0, atol=1e-12)
    return {'user_curves': len(user), 'matched_curve_keys': int(merged._merge.eq('both').sum()),
            'matching_response_curves': int(same.sum()), 'all_user_response_fields_match': bool(same.all()),
            'tolerance_absolute': 1e-12, 'annotation_identity_claimed': False}, merged.loc[~same]


def main(root=ROOT, reuse_annotations=None):
    root = Path(root)
    external = root / 'data/external/gdsc_release8_5'
    archive_path = external / 'GDSC2_public_raw_data_27Oct23.csv.zip'
    with ResearchRun(root, 'gdsc_raw_preparation') as run:
        verify_run(root / 'reports/gdsc_source_acquisition')
        run.inputs([archive_path, external / 'GDSC2_fitted_dose_response_27Oct23.xlsx',
                    external / 'screened_compounds_rel_8.5.csv', root / 'data/raw' / FILES['glioma_csv']])
        user = read_responses(root / 'data/raw' / FILES['glioma_csv'])
        official = pd.read_excel(external / 'GDSC2_fitted_dose_response_27Oct23.xlsx')
        comparison, conflicts = compare_fitted(user.copy(), official)
        conflicts.to_csv(run.path / 'fitted_conflicts.csv', index=False)
        write_json(run.path / 'fitted_comparison.json', comparison)
        wanted = set(user.SANGER_MODEL_ID)
        if reuse_annotations is not None:
            reuse_annotations = Path(reuse_annotations).resolve()
            prior = verify_run(reuse_annotations.parent)
            key = str(archive_path.relative_to(root))
            if prior['inputs'].get(key, {}).get('sha256') != checksum(archive_path):
                raise ValueError('La extracción anterior pertenece a otro ZIP.')
            prior_summary = json.loads((reuse_annotations.parent / 'run.json').read_text())
            if not prior_summary['raw_zip_crc_checked_by_full_member_read']:
                raise ValueError('La extracción anterior no acredita lectura completa/CRC.')
            run.inputs([reuse_annotations])
            raw = pd.read_parquet(reuse_annotations)
            if set(raw.SANGER_MODEL_ID) != wanted:
                raise ValueError('La extracción anterior no corresponde a las líneas actuales.')
            total = prior_summary['raw_annotation_rows_scanned']
        else:
            selected, total = [], 0
            with zipfile.ZipFile(archive_path) as z:
                names = [x for x in z.namelist() if x.endswith('.csv')]
                if len(names) != 1:
                    raise ValueError('Se requiere un único CSV crudo identificado.')
                with z.open(names[0]) as stream:
                    chunks = pd.read_csv(stream, chunksize=200000,
                        dtype={'SANGER_MODEL_ID': 'string', 'DRUG_ID': 'Int64', 'BARCODE': 'string', 'SCAN_ID': 'string'})
                    for index, chunk in enumerate(chunks, 1):
                        total += len(chunk)
                        if 'SANGER_MODEL_ID' not in chunk:
                            raise ValueError('No hay SIDM; no se admite unión aproximada por nombres.')
                        part = chunk[chunk.SANGER_MODEL_ID.isin(wanted)].copy()
                        if len(part):
                            selected.append(part)
                        if index % 10 == 0:
                            print(f'Leídas {total:,} anotaciones de pozos', flush=True)
            if not selected:
                raise ValueError('No hay líneas GBM/glioma emparejadas en el crudo.')
            raw = pd.concat(selected, ignore_index=True)
        raw.DRUG_ID = raw.DRUG_ID.astype('string')
        raw.to_parquet(run.path / 'glioma_raw_annotations.parquet', index=False)
        wells = collapse_annotations(raw)
        wells, controls = normalize_controls(wells)
        curves = crossing_diagnostics(wells)
        keys = ['SANGER_MODEL_ID', 'DRUG_ID']
        if user.duplicated(keys).any():
            raise ValueError('Respuestas locales ambiguas por línea/fármaco.')
        if comparison['all_user_response_fields_match']:
            curves = curves.merge(user[keys+['NLME_CURVE_ID', 'DRUG_NAME', 'LN_IC50', 'MIN_CONC', 'MAX_CONC']],
                                  on=keys, how='left', validate='many_to_one')
        wells.to_parquet(run.path / 'glioma_normalized_wells.parquet', index=False)
        controls.to_csv(run.path / 'plate_controls.csv', index=False)
        curves.to_csv(run.path / 'empirical_curves.csv', index=False)
        curves[curves.SANGER_MODEL_ID.eq('SIDM01171')].to_csv(run.path / 't98g_empirical_curves.csv', index=False)
        pd.DataFrame({'SANGER_MODEL_ID': sorted(wanted-set(raw.SANGER_MODEL_ID))}).to_csv(run.path / 'unmatched_lines.csv', index=False)
        summary = {'run_id': run.id, 'source_release': '8.5', 'raw_annotation_rows_scanned': total,
            'selected_annotation_rows': len(raw), 'selected_physical_wells': len(wells),
            'matched_lines': raw.SANGER_MODEL_ID.nunique(), 'selected_plates': controls.shape[0],
            'plates_with_explicit_usable_nc1_b': int(controls.status.eq('normalized_explicit_nc1_b').sum()),
            'empirical_plate_drug_series': len(curves), 'fitted_comparison': comparison,
            'raw_zip_crc_checked_by_full_member_read': True, 'nlme_refit_performed': False,
            'verified_annotation_reuse': str(reuse_annotations) if reuse_annotations else None,
            'interpretation': 'Descriptive raw-well diagnostics, preserving plate/campaign dependence; not a replacement for original NLME fitting, a calibrated censoring model, or validation of nanoformulations.'}
        write_json(run.path / 'run.json', summary)
        print(json.dumps(summary, indent=2, ensure_ascii=False, default=int))
        return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--reuse-annotations', type=Path)
    args = parser.parse_args()
    main(reuse_annotations=args.reuse_annotations)
