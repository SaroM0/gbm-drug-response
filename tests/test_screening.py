"""Regresiones de la auditoría; fixtures aisladas, sin red ni métricas científicas."""
import io
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

from gbm_drug_response.data import range_state, quality_flags, SCHEMA
from gbm_drug_response.screening import (load_screening, potency_table, required_loading_pct,
    summarize_compounds, scenario_table, cumulative_funnel)
from gbm_drug_response.pilot import protocol_readiness, assess_formulations, chemistry_reviews
from gbm_drug_response.runs import ResearchRun, checksum
from gbm_drug_response.runs import verify_run
from gbm_drug_response.validation import annotate_validation_groups, outcome_diagnostics
from gbm_drug_response.modeling import attach_split, grouped_folds, metrics

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / 'scripts'))
from prepare_expression import read_tpm_member


@pytest.fixture
def cfg():
    return load_screening(PROJECT)


def responses(values):
    return pd.DataFrame([dict(zip(SCHEMA, ['GDSC2', 1, i+1, 'Fixture', 'SIDM00001',
        'Glioblastoma', str(i+1), 'AuditOnly', 'TEST', 'TEST', .001, 10,
        v, .5, .1, 0])) for i, v in enumerate(values)])


def assays():
    return pd.DataFrame([dict(compound_id='same-active-entity', gdsc_drug_id=str(i+1),
        drug_name='Fixture', t98g_ln_ic50=ln, t98g_ic50_um=np.exp(ln),
        t98g_in_range=within, t98g_range_state='within_range' if within else 'above_range',
        active_molecular_weight=500., structure_ok=True, conjugation_route='covalente_cooh_edc_nhs',
        review_flags='pending') for i, (ln, within) in enumerate([(0., True), (2., False)])])


def test_range_state_rejects_both_extrapolation_directions():
    df = responses([-10., 0., 4.])
    assert range_state(df).tolist() == ['below_range', 'within_range', 'above_range']
    # Estas banderas no eliminan silenciosamente ajustes del entrenamiento histórico.
    assert quality_flags(df).keep.all()


def test_potency_table_checks_lower_boundary_and_panel_fractions(tmp_path):
    (tmp_path / 'config.json').write_bytes((PROJECT / 'config.json').read_bytes())
    df = responses([-10., 0., 4.])
    table = potency_table(tmp_path, {k: df.copy() for k in ['t98g', 'glioma_csv', 'gdsc2']})
    assert table.t98g_in_range.tolist() == [False, True, False]
    assert table.glioma_frac_in_range.tolist() == [0., 1., 0.]


def test_compound_summary_keeps_worst_and_all_assays_not_best(cfg):
    table = summarize_compounds(assays(), cfg)
    row = table.iloc[0]
    assert row.all_gdsc_drug_ids == '1;2'
    assert row.n_assays == 2 and row.t98g_ln_ic50 == 2.
    assert not row.t98g_in_range and row.tier == 8
    assert row.assay_disagreement_fold == pytest.approx(np.exp(2))
    assert 'curve_rmse' not in table and 't98g_auc' not in table


def test_scenarios_distinguish_bracketing_from_reaching_maximum(cfg):
    a = assays().iloc[[0]].copy()
    a['t98g_ln_ic50'] = np.log(.005)
    a['t98g_ic50_um'] = .005
    table = scenario_table(summarize_compounds(a, cfg), cfg)
    ref = table[(table.loading_pct == 5) & (table.released_fraction == 1)].iloc[0]
    assert ref.passes_loading_scenario
    assert ref.bracketing_status == 'all_doses_above_reference'
    assert not ref.experimental_prediction


def test_loading_accounts_for_partial_release_and_validates_units():
    assert float(required_loading_pct(1., 500., 25., .1)) == pytest.approx(20.)
    assert float(required_loading_pct(1., 500., 25., 1.)) == pytest.approx(2.)
    assert np.isnan(required_loading_pct(1., -500.))
    for release in [0, -1, 1.1]:
        with pytest.raises(ValueError):
            required_loading_pct(1., 500., 25., release)


def test_screening_rejects_unsupported_policies_and_invalid_release(tmp_path):
    d = json.loads((PROJECT / 'screening.json').read_text())
    d['released_fraction_scenarios'] = [0.]
    (tmp_path / 'screening.json').write_text(json.dumps(d))
    with pytest.raises(ValueError):
        load_screening(tmp_path)


def test_funnel_is_cumulative_and_preserves_units():
    df = pd.DataFrame({'compound_id': ['a', 'a', None, 'b'],
        'structure_ok': [True, True, False, True], 'compound_all_in_range': [True, True, True, False],
        'compound_loading_pass': [True, True, True, True], 'compound_route_pass': [True, True, True, True]})
    funnel = cumulative_funnel(df)
    assert funnel.drug_id_rows.tolist() == [4, 3, 2, 2, 2]
    assert funnel.distinct_nonmissing_compound_ids.tolist() == [2, 2, 1, 1, 1]


def test_draft_protocol_and_empty_measurements_are_not_approval(tmp_path):
    (tmp_path / 'data').mkdir()
    for name in ['candidate_reviews.csv', 'formulation_measurements.csv']:
        (tmp_path / 'data' / name).write_bytes((PROJECT / 'data' / name).read_bytes())
    p = json.loads((PROJECT / 'pilot_protocol.json').read_text())
    (tmp_path / 'pilot_protocol.json').write_text(json.dumps(p))
    assert not protocol_readiness(p)['ready']
    assert 'primary_time_h' in protocol_readiness(p)['missing_decisions']
    assert chemistry_reviews(tmp_path) == set()
    assert assess_formulations(tmp_path, set()).empty


def test_run_preserves_legacy_and_does_not_expose_previous_model(tmp_path):
    old = tmp_path / 'reports/molecular'; old.mkdir(parents=True)
    (old / 'test_metrics.csv').write_text('legacy metrics')
    model = tmp_path / 'models/molecular.joblib'; model.parent.mkdir()
    model.write_text('legacy model')
    with ResearchRun(tmp_path, 'molecular') as run:
        run.model_link('molecular.joblib')
        assert not model.exists()
        assert not (tmp_path / 'reports/molecular/test_metrics.csv').exists()
        (run.path / 'status.json').write_text('{"status":"missing_expression"}')
        first = run.path
    assert any(p.read_text() == 'legacy model' for p in (tmp_path / 'reports/runs/legacy_models').iterdir())
    assert any(p.read_text() == 'legacy metrics' for p in (tmp_path / 'reports/runs/legacy').rglob('test_metrics.csv'))
    with pytest.raises(ValueError):
        with ResearchRun(tmp_path, 'molecular') as failed:
            raise ValueError('bad new expression')
    assert json.loads((failed.path / 'run_manifest.json').read_text())['status'] == 'failed'
    assert (first / 'status.json').exists()
    assert not (tmp_path / 'reports/molecular/status.json').exists()


def test_related_models_cannot_cross_holdout_or_cv(tmp_path):
    (tmp_path / 'data/external').mkdir(parents=True)
    meta = pd.DataFrame({'model_id': ['SIDM1', 'SIDM2', 'SIDM3', 'SIDM4'],
        'patient_id': ['P1', 'P1', 'P3', 'P4'], 'parent_id': [None, 'SIDM1', None, None]})
    meta.to_csv(tmp_path / 'data/external/model_list_20260921.csv', index=False)
    df = pd.DataFrame({'SANGER_MODEL_ID': meta.model_id})
    grouped, info = annotate_validation_groups(tmp_path, df, tmp_path)
    split = pd.DataFrame({'SANGER_MODEL_ID': meta.model_id, 'partition': ['development', 'test', 'development', 'test']})
    with pytest.raises(ValueError, match='donante'):
        attach_split(grouped, split)
    for tr, va in grouped_folds(grouped, 2):
        assert set(grouped.iloc[tr].VALIDATION_GROUP).isdisjoint(grouped.iloc[va].VALIDATION_GROUP)


def test_outcome_diagnostics_do_not_drop_extrapolated_rows(tmp_path):
    df = responses([-10., 0., 4.]); df['model'] = 0.
    status = outcome_diagnostics(df, ['model'], metrics, tmp_path)
    report = pd.read_csv(tmp_path / 'outcome_diagnostics.csv')
    assert report.loc[report.range_group == 'all_fitted_estimates', 'n'].iloc[0] == 3
    assert status['n_within_range'] == status['n_above_range'] == status['n_below_range'] == 1
    assert not status['selection_uses_test_subgroups']


def tpm_bytes(a=400000, b=600000, duplicate=False):
    return io.BytesIO(('model_id,,SIDM1\nmodel_name,,Fixture\ndataset_name,,Merged\ndata_source,,Sanger\n'
        f'gene_id,symbol,\nSIDG1,G1,{a}\n' + f'{"SIDG1" if duplicate else "SIDG2"},G2,{b}\n').encode())


def test_expression_transform_once_and_keep_source_profiles():
    expr, profiles, genes, missing = read_tpm_member(tpm_bytes(), ['SIDM1', 'SIDM2'])
    assert missing == ['SIDM2']
    assert expr.gene_SIDG1.iloc[0] == pytest.approx(np.log2(400001))
    assert profiles.data_source.iloc[0] == 'Sanger'
    assert genes.gene_symbol.tolist() == ['G1', 'G2']


def test_expression_rejects_already_logged_values_and_duplicate_genes():
    with pytest.raises(ValueError, match='escala'):
        read_tpm_member(tpm_bytes(18., 19.), ['SIDM1'])
    with pytest.raises(ValueError, match='duplicados'):
        read_tpm_member(tpm_bytes(duplicate=True), ['SIDM1'])


def test_verified_measurements_can_pass_but_bad_source_cannot(tmp_path):
    (tmp_path / 'data').mkdir()
    p = json.loads((PROJECT / 'pilot_protocol.json').read_text())
    from gbm_drug_response.pilot import REQUIRED_DECISIONS
    p.update({k: 'Fixture decision, not laboratory data' for k in REQUIRED_DECISIONS})
    p.update(status='approved', primary_time_h=72, independent_experiments=3,
             minimum_relevant_improvement=.1)
    (tmp_path / 'pilot_protocol.json').write_text(json.dumps(p))
    source = tmp_path / 'fixture.txt'; source.write_text('Unit test evidence only')
    table = pd.read_csv(PROJECT / 'data/formulation_measurements.csv')
    row = {k: 'true' for k in table.columns}
    row.update(formulation_id='fixture', compound_id='fixture-compound', batch_id='fixture-batch',
        loading_pct=5, loading_basis='active_parent_drug_mass_per_carrier_mass', released_fraction=.5,
        release_time_h=72, independent_experiments=3, source_file='fixture.txt', source_sha256=checksum(source))
    pd.DataFrame([row]).to_csv(tmp_path / 'data/formulation_measurements.csv', index=False)
    result = assess_formulations(tmp_path, {'fixture-compound'})
    assert result.ready_for_t98g_comparison.iloc[0]
    assert result.ready_for_transport_claim.iloc[0]
    source.write_text('Evidence changed')
    result = assess_formulations(tmp_path, {'fixture-compound'})
    assert not result.ready_for_t98g_comparison.iloc[0]
    assert 'source_not_verified' in result.pending_or_failed.iloc[0]


def test_manifest_detects_changed_outputs(tmp_path):
    with ResearchRun(tmp_path, 'test') as run:
        (run.path / 'result.csv').write_text('original')
    assert verify_run(run.path)['status'] == 'complete'
    (run.path / 'result.csv').write_text('altered')
    with pytest.raises(ValueError, match='alterada'):
        verify_run(run.path)


def test_notebook_socket_restriction_fails_before_kernel_start(monkeypatch):
    import run_notebooks
    def denied():
        raise PermissionError('fixture sandbox')
    monkeypatch.setattr(run_notebooks.socket, 'socket', denied)
    with pytest.raises(RuntimeError, match='sockets locales'):
        run_notebooks.check_kernel_transport()
