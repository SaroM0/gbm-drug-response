"""Cribado por molécula con ensayos completos y sensibilidad de dosis/liberación."""
from pathlib import Path
import json
import sys
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from gbm_drug_response.data import FILES
from gbm_drug_response.runs import ResearchRun, write_json
from gbm_drug_response.screening import (required_loading_pct, potency_table, assign_tier,
    load_screening, summarize_compounds, scenario_table, cumulative_funnel, review_flags)
from gbm_drug_response.pilot import readiness, chemistry_reviews, assess_formulations

NP_PPM_MAX = 25.0  # Sólo default de la función; la corrida lee screening.json.


def main(root=ROOT):
    root = Path(root)
    with ResearchRun(root, 'candidates') as run:
        run.inputs([root / 'data/raw' / FILES[k] for k in ['t98g', 'glioma_csv', 'gdsc2']] +
                   [root / 'data/processed/compounds.csv', root / 'screening.json',
                    root / 'pilot_protocol.json', root / 'data/candidate_reviews.csv',
                    root / 'data/formulation_measurements.csv'])
        cfg = load_screening(root)
        compounds = pd.read_csv(root / 'data/processed/compounds.csv', dtype={'gdsc_drug_id': 'string'})
        potency = potency_table(root)
        assays = potency.merge(compounds, on='gdsc_drug_id', how='left', validate='one_to_one', indicator=True)
        if not assays._merge.eq('both').all():
            raise ValueError('Ensayos de T98G sin fila química, incluso para estructuras no resueltas.')
        assays = assays.drop(columns='_merge')
        assays['structure_ok'] = (assays.parse_status.eq('ok') & assays.compound_id.notna()
            & ~assays.ambiguous_name.fillna(True).astype(bool)
            & assays.normalization_status.eq('standardized'))
        assays['review_flags'] = assays.apply(review_flags, axis=1)
        ranking = summarize_compounds(assays, cfg)
        accepted = chemistry_reviews(root)
        formulations = assess_formulations(root, accepted)
        ranking['chemistry_review_approved'] = ranking.compound_id.isin(accepted)
        ready_compounds = set(formulations.loc[formulations.ready_for_t98g_comparison, 'compound_id'])
        ranking['comparative_assay_ready'] = ranking.compound_id.isin(ready_compounds)
        ranking['decision_status'] = np.where(ranking.chemistry_review_approved,
            'chemistry_reviewed_formulation_and_assay_validation_pending', 'manual_chemistry_review_pending')
        short = ranking[ranking.tier.le(4)].copy()
        scenarios = scenario_table(ranking, cfg)
        identity = ranking.dropna(subset=['compound_id']).set_index('compound_id')
        assays['compound_all_in_range'] = assays.compound_id.map(identity.t98g_in_range).astype('boolean').fillna(False).astype(bool)
        assays['compound_loading_pass'] = assays.compound_id.map(identity.required_loading_pct).le(cfg['screening_loading_limit_pct'])
        assays['compound_route_pass'] = assays.compound_id.map(identity.tier).le(4)
        funnel = cumulative_funnel(assays)
        sensitivity = scenarios.groupby(['loading_pct', 'released_fraction'], as_index=False).agg(
            compound_scenarios=('compound_id', 'size'),
            passing_mass_scenario=('passes_loading_scenario', 'sum'),
            bracketing_reference=('bracketing_status', lambda s: int(s.eq('brackets_reference').sum())))
        for name, table in [('assays', assays), ('ranking', ranking), ('shortlist', short),
                            ('dose_scenarios', scenarios), ('scenario_summary', sensitivity), ('funnel', funnel),
                            ('formulation_readiness', formulations)]:
            table.to_csv(run.path / f'{name}.csv', index=False)
        summary = {'run_id': run.id, 'protocol': cfg, 'purpose': 'Exploratory prioritization, not formulation efficacy prediction',
            'qc_assay_rows': len(assays), 't98g_above_range': int(assays.t98g_range_state.eq('above_range').sum()),
            't98g_below_range': int(assays.t98g_range_state.eq('below_range').sum()),
            't98g_within_range': int(assays.t98g_in_range.sum()),
            'ranked_entities_including_unresolved': len(ranking), 'shortlist_molecules': len(short),
            'chemistry_reviews_approved': int(ranking.chemistry_review_approved.sum()),
            'pilot_protocol': readiness(root), 'comparative_assay_ready': bool(ranking.comparative_assay_ready.any()),
            'unmeasured_requirements': ['loading per formulation batch', 'release over time in assay conditions',
                'retained activity', 'MTT interference and orthogonal validation', 'dose equivalence',
                'endothelial safety and intact-barrier transport'],
            'funnel': funnel.to_dict('records')}
        write_json(run.path / 'run.json', summary)
        print(json.dumps(summary, ensure_ascii=False, indent=2))
        return summary


if __name__ == '__main__':
    main()
