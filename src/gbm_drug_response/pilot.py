"""Decisiones experimentales explícitas; no inventa respuestas del laboratorio."""
import json
from pathlib import Path
import pandas as pd
import numpy as np
from .runs import checksum

REQUIRED_DECISIONS = [
    'primary_time_h', 'minimum_relevant_improvement', 'independent_experiments',
    'sample_size_rationale', 'normalization_and_curve_model', 'multiple_comparisons_plan',
    'randomization_and_plate_layout', 'orthogonal_assay', 'physical_mixture_control_decision',
    'endothelial_acceptance_criterion', 'transwell_analyte', 'barrier_integrity_criterion',
    'transport_mass_balance_criterion', 'approved_by', 'approved_at',
]


def protocol_readiness(protocol):
    missing = [key for key in REQUIRED_DECISIONS if protocol.get(key) is None or protocol.get(key) == '']
    errors = []
    if protocol.get('status') != 'approved':
        errors.append('protocol_not_approved')
    if protocol.get('primary_time_h') is not None and protocol['primary_time_h'] not in protocol.get('planned_times_h', []):
        errors.append('primary_time_not_planned')
    if protocol.get('independent_experiments') is not None:
        n = protocol['independent_experiments']
        if not isinstance(n, int) or isinstance(n, bool) or n < 2:
            errors.append('independent_experiments_requires_integer_at_least_2_and_rationale')
    value = protocol.get('minimum_relevant_improvement')
    if value is not None and (not isinstance(value, (int, float)) or not 0 < value <= 1):
        errors.append('improvement_requires_fraction_in_0_1')
    for key in ['vehicle_control', 'empty_nanoparticle_control', 'cell_free_interference_controls']:
        if protocol.get(key) is not True:
            errors.append(key + '_required')
    return {'ready': not missing and not errors, 'missing_decisions': missing, 'errors': errors,
            'meaning': 'Decisions only; not evidence that controls or measurements were performed.'}


def chemistry_reviews(root):
    table = pd.read_csv(Path(root) / 'data/candidate_reviews.csv', dtype='string')
    required = ['compound_id', 'identity_verified', 'route_verified', 'pharmacophore_review',
                'supplier_form_id', 'reviewer', 'reviewed_at', 'notes']
    if list(table.columns) != required or table.compound_id.isna().any() or table.compound_id.duplicated().any():
        raise ValueError('Registro químico inválido o con identidades duplicadas.')
    accepted = set()
    for row in table.itertuples():
        flags = all(str(getattr(row, k)).lower() == 'true' for k in
                    ['identity_verified', 'route_verified', 'pharmacophore_review'])
        signed = all(pd.notna(getattr(row, k)) and str(getattr(row, k)).strip() for k in
                     ['supplier_form_id', 'reviewer', 'reviewed_at'])
        if flags and signed:
            accepted.add(row.compound_id)
    return accepted


def readiness(root):
    return protocol_readiness(json.loads((Path(root) / 'pilot_protocol.json').read_text()))


def assess_formulations(root, reviewed):
    """Exige mediciones trazables; distingue comparación celular de transporte."""
    root = Path(root).resolve()
    protocol = json.loads((root / 'pilot_protocol.json').read_text())
    ready = protocol_readiness(protocol)
    data = pd.read_csv(root / 'data/formulation_measurements.csv', dtype='string')
    required = ['formulation_id', 'compound_id', 'batch_id', 'loading_pct', 'loading_basis',
        'released_fraction', 'release_time_h', 'release_medium', 'active_drug_equivalence_verified',
        'interference_controls_passed', 'orthogonal_assay_passed', 'dose_range_brackets_half_response',
        'independent_experiments', 'endothelial_criterion_passed', 'barrier_integrity_passed',
        'transport_mass_balance_passed', 'source_file', 'source_sha256', 'reviewer', 'reviewed_at']
    if list(data.columns) != required or data.formulation_id.isna().any() or data.formulation_id.duplicated().any():
        raise ValueError('Registro de formulaciones inválido o IDs repetidos.')
    rows = []
    for _, row in data.iterrows():
        reasons = []
        if not ready['ready']:
            reasons.append('protocol_pending')
        if row.compound_id not in reviewed:
            reasons.append('chemistry_review_pending')
        if row[required].isna().any():
            reasons.append('missing_measurements_or_provenance')
        for key in ['loading_pct', 'released_fraction', 'release_time_h', 'independent_experiments']:
            value = pd.to_numeric(row[key], errors='coerce')
            if not np.isfinite(value) or value <= 0:
                reasons.append(key + '_invalid')
        release = pd.to_numeric(row.released_fraction, errors='coerce')
        if np.isfinite(release) and release > 1:
            reasons.append('released_fraction_above_1')
        count = pd.to_numeric(row.independent_experiments, errors='coerce')
        if np.isfinite(count) and (count != int(count) or count < (protocol.get('independent_experiments') or 2)):
            reasons.append('insufficient_independent_experiments')
        if pd.to_numeric(row.release_time_h, errors='coerce') != protocol.get('primary_time_h'):
            reasons.append('release_time_does_not_match_primary_time')
        if str(row.loading_basis) != 'active_parent_drug_mass_per_carrier_mass':
            reasons.append('wrong_mass_basis')
        for key in ['active_drug_equivalence_verified', 'interference_controls_passed',
                    'orthogonal_assay_passed', 'dose_range_brackets_half_response']:
            if str(row[key]).lower() != 'true':
                reasons.append(key + '_pending_or_failed')
        source = (root / str(row.source_file)).resolve()
        verified = source.is_relative_to(root) and source.is_file() and checksum(source) == row.source_sha256
        if not verified:
            reasons.append('measurement_source_not_verified')
        cellular_ready = not reasons
        transport_ready = cellular_ready and all(str(row[k]).lower() == 'true' for k in
            ['endothelial_criterion_passed', 'barrier_integrity_passed', 'transport_mass_balance_passed'])
        rows.append({'formulation_id': row.formulation_id, 'compound_id': row.compound_id,
                     'ready_for_t98g_comparison': cellular_ready, 'ready_for_transport_claim': transport_ready,
                     'pending_or_failed': ';'.join(reasons)})
    return pd.DataFrame(rows, columns=['formulation_id', 'compound_id', 'ready_for_t98g_comparison',
                                      'ready_for_transport_claim', 'pending_or_failed'])
