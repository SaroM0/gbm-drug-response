"""Materializa un piloto propuesto con aleatorización y sensibilidad de tamaño muestral."""
import json
from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from gbm_drug_response.planning import plate_layout, power_scenarios
from gbm_drug_response.runs import ResearchRun, verify_run, write_json


def main(root=ROOT):
    root = Path(root)
    with ResearchRun(root, 'pilot_preparation') as run:
        verify_run(root / 'reports/candidates')
        run.inputs([root / 'pilot_design.json', root / 'pilot_protocol.json', root / 'reports/candidates/ranking.csv'])
        design = json.loads((root / 'pilot_design.json').read_text())
        ranking = pd.read_csv(root / 'reports/candidates/ranking.csv')
        candidates = ranking[ranking.drug_name.isin(design['candidate_names'])].copy()
        if len(candidates) != len(design['candidate_names']) or not candidates.compound_id.is_unique:
            raise ValueError('Candidatos ausentes o ambiguos; revisar la propuesta antes del diseño.')
        layout = plate_layout(design, candidates)
        power = power_scenarios(design['power_scenarios'])
        layout.to_csv(run.path / 'plate_layout.csv', index=False)
        power.to_csv(run.path / 'sample_size_scenarios.csv', index=False)
        candidates[['compound_id', 'drug_name', 'pubchem_cid', 'pubchem_title', 'form_inchikey',
                    'active_molecular_weight', 'active_smiles', 'all_gdsc_drug_ids', 'conjugation_route', 'review_flags',
                    'ic50_min_um', 'ic50_max_um', 'required_loading_pct']].to_csv(run.path / 'candidate_review_packet.csv', index=False)
        summary = {'run_id': run.id, 'status': 'planning_only_not_laboratory_approved',
            'design': design, 'n_candidates': len(candidates), 'n_plates': layout.plate_id.nunique(),
            'planned_well_slots': len(layout), 'planned_measurement_wells': int(layout.include_in_measurement_plan.sum()),
            'independent_campaigns': design['independent_experiments_for_variance_pilot'],
            'formal_power_claimed': False, 'observations_generated': 0,
            'notes': ['Numeric active-drug doses require measured loading, supplied separately.',
                      'A positive-control reagent and readout protocol require laboratory selection.',
                      'Matched cell-free controls diagnose interference, not guarantee additivity or valid subtraction.',
                      'This proposal does not approve candidates or modify the holdout or trained models.']}
        write_json(run.path / 'run.json', summary)
        print(f'Propuesta preparada: {len(candidates)} candidatos, {summary["n_plates"]} placas, cero mediciones inventadas.')
        return summary


if __name__ == '__main__':
    main()
