import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import pytest

from gbm_drug_response.planning import paired_power, power_scenarios, plate_layout

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from acquire_gdsc_sources import catalog, download
from gbm_drug_response.gdsc_raw import collapse_annotations, normalize_controls, crossing_diagnostics
from prepare_gdsc_raw import compare_fitted


def test_power_matches_null_alpha_and_increases_with_independent_n():
    assert paired_power(10, 0., .4, .05) == pytest.approx(.05)
    assert paired_power(20, .3, .4, .05) > paired_power(3, .3, .4, .05)
    assert paired_power(10, -.3, .4, .05) == paired_power(10, .3, .4, .05)
    with pytest.raises(ValueError):
        paired_power(3, .3, 0., .05)


def test_power_scenarios_do_not_claim_measured_sd_or_sufficient_triplicates():
    cfg = json.loads((ROOT / 'pilot_design.json').read_text())['power_scenarios']
    table = power_scenarios(cfg)
    assert len(table) == 15 and not table.measured_variance.any()
    row = table[(table.loaded_to_free_ic50_ratio == .7) & (table.hypothetical_paired_log_ratio_sd == .4)].iloc[0]
    assert row.n_for_power_80pct == 17
    assert row.power_if_n_3 < .1


def test_randomized_plate_pairs_preserve_controls_and_empty_measurements():
    design = json.loads((ROOT / 'pilot_design.json').read_text())
    design.update(candidate_names=['Fixture'], times_h=[72], independent_experiments_for_variance_pilot=1)
    candidates = pd.DataFrame({'drug_name': ['Fixture'], 'compound_id': ['fixture']})
    first, second = plate_layout(design, candidates), plate_layout(design, candidates)
    pd.testing.assert_frame_equal(first, second)
    assert first.plate_id.nunique() == 2
    assert not first.duplicated(['plate_id', 'well']).any()
    assert first.raw_readout.isna().all() and first.active_drug_conc_um.isna().all()
    assert not first.laboratory_approved.any()
    a = first[first.plate_type == 'cellular'].set_index('well')
    b = first[first.plate_type == 'cell_free_interference'].set_index('well')
    pd.testing.assert_series_equal(a.treatment.sort_index(), b.treatment.sort_index())
    assert not b.cells_present.any()
    assert a[a.treatment == 'vehicle'].shape[0] == 6
    assert first.groupby('plate_id').include_in_measurement_plan.sum().eq(72).all()


def test_layout_refuses_overflow_instead_of_dropping_controls():
    design = json.loads((ROOT / 'pilot_design.json').read_text())
    design['technical_replicates'] = 10
    with pytest.raises(ValueError, match='capacidad'):
        plate_layout(design, pd.DataFrame({'drug_name': ['Fixture'], 'compound_id': ['fixture']}))


def test_catalog_requires_complete_listing():
    xml = '<ListBucketResult xmlns="http://s3.amazonaws.com/doc/2006-03-01/"><IsTruncated>true</IsTruncated></ListBucketResult>'
    with pytest.raises(ValueError, match='incompleto'):
        catalog(xml)


def test_download_checks_size_and_reuses_only_verified_cache(tmp_path, monkeypatch):
    class Response:
        url = 'https://cog.sanger.ac.uk/cancerrxgene/fixture'
        headers = {'ETag': 'fixture-etag'}
        def __enter__(self): return self
        def __exit__(self, *args): pass
        def raise_for_status(self): pass
        def iter_content(self, size): yield b'fixture-bytes'
    monkeypatch.setattr('acquire_gdsc_sources.requests.get', lambda *a, **k: Response())
    target = tmp_path / 'file'
    receipt = download(Response.url, target, len(b'fixture-bytes'), 'fixture-etag')
    assert target.read_bytes() == b'fixture-bytes'
    assert download(Response.url, target, len(b'fixture-bytes'), 'fixture-etag')['cache_reused']
    with pytest.raises(ValueError, match='alterada'):
        download(Response.url, target, len(b'fixture-bytes'), 'new-etag')
    with pytest.raises(ValueError, match='incompleta'):
        download(Response.url, tmp_path / 'short', len(b'fixture-bytes')+1)
    assert not (tmp_path / 'short').exists()


def raw_fixture():
    rows = []
    conditions = [('NC-1', None, None, 110.), ('NC-1', None, None, 110.),
        ('B', None, None, 10.), ('B', None, None, 10.),
        ('L1-D3-S', '1007', .1, 100.), ('L1-D2-S', '1007', .3, 70.),
        ('L1-D1-S', '1007', 1., 30.)]
    for position, (tag, drug, dose, intensity) in enumerate(conditions, 1):
        rows.append({'BARCODE': 'plate', 'SCAN_ID': 'scan', 'POSITION': position,
            'SANGER_MODEL_ID': 'SIDM1', 'DRUG_ID': drug, 'CONC': dose,
            'INTENSITY': intensity, 'TAG': tag})
    # Segundo tag, mismo pozo: no es otra réplica ni otra intensidad.
    rows.append({**rows[-1], 'TAG': 'DMSO', 'DRUG_ID': '4000', 'CONC': None})
    return pd.DataFrame(rows)


def test_raw_annotations_do_not_duplicate_wells_and_controls_normalize():
    raw = raw_fixture(); wells = collapse_annotations(raw)
    assert len(wells) == 7
    assert wells.loc[wells.POSITION.eq(7), 'annotation_rows'].iloc[0] == 2
    normalized, plates = normalize_controls(wells)
    last = normalized[normalized.POSITION.eq(7)].iloc[0]
    assert last.normalized_viability == pytest.approx(.2)
    assert plates.n_vehicle_wells.iloc[0] == 2
    curves = crossing_diagnostics(normalized)
    assert curves.n_doses.iloc[0] == 3
    assert curves.empirical_range_state.iloc[0] == 'points_on_both_sides_of_half'
    assert not curves.nlme_refit_performed.iloc[0]
    assert not wells.unsupported_combination.any()


def test_raw_unused_wells_are_not_automatically_blank_controls():
    raw = raw_fixture(); raw.loc[raw.TAG.eq('B'), 'TAG'] = 'UN-USED'
    normalized, plates = normalize_controls(collapse_annotations(raw))
    assert normalized.normalized_viability.isna().all()
    assert plates.status.iloc[0] == 'missing_or_invalid_explicit_controls'


def test_conflicting_physical_well_intensity_is_rejected():
    raw = raw_fixture(); raw.loc[raw.index[-1], 'INTENSITY'] = 999.
    with pytest.raises(ValueError, match='contradictorios'):
        collapse_annotations(raw)


def test_fit_linkage_does_not_silently_accept_changed_outcomes():
    row = dict(NLME_CURVE_ID='1', SANGER_MODEL_ID='SIDM1', DRUG_ID='1007',
        LN_IC50=0., MIN_CONC=.01, MAX_CONC=1., AUC=.5, RMSE=.1)
    user, official = pd.DataFrame([row]), pd.DataFrame([row])
    assert compare_fitted(user.copy(), official.copy())[0]['all_user_response_fields_match']
    official.LN_IC50 = 2.
    summary, conflicts = compare_fitted(user, official)
    assert not summary['all_user_response_fields_match']
    assert len(conflicts) == 1
