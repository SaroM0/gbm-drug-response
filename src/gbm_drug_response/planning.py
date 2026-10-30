"""Diseño de placas y potencia por escenarios; no produce mediciones ficticias."""
import numpy as np
import pandas as pd
from scipy.stats import t, nct


def paired_power(n, effect, paired_sd, alpha):
    if not isinstance(n, int) or isinstance(n, bool) or n < 2:
        raise ValueError('n cuenta experimentos independientes y debe ser entero >=2.')
    if not all(np.isfinite(x) for x in [effect, paired_sd, alpha]) or paired_sd <= 0 or not 0 < alpha < 1:
        raise ValueError('Efecto finito, SD positiva y alpha en (0,1) requeridos.')
    critical = t.ppf(1-alpha/2, n-1)
    noncentrality = abs(effect) * np.sqrt(n) / paired_sd
    return float(nct.cdf(-critical, n-1, noncentrality) + nct.sf(critical, n-1, noncentrality))


def power_scenarios(cfg):
    alpha = cfg['familywise_alpha'] / cfg['candidate_comparisons']
    rows = []
    for ratio in cfg['loaded_to_free_ic50_ratios']:
        if not 0 < ratio < 1:
            raise ValueError('Las razones propuestas deben estar en (0,1).')
        for sd in cfg['sd_of_paired_log_ratio']:
            powers = [(n, paired_power(n, np.log(ratio), sd, alpha))
                      for n in range(2, cfg['max_independent_experiments']+1)]
            row = {'loaded_to_free_ic50_ratio': ratio, 'hypothetical_paired_log_ratio_sd': sd,
                   'per_comparison_alpha': alpha, 'power_if_n_3': paired_power(3, np.log(ratio), sd, alpha),
                   'measured_variance': False}
            for target in cfg['power_targets']:
                row[f'n_for_power_{int(target*100)}pct'] = next((n for n, p in powers if p >= target), None)
            rows.append(row)
    return pd.DataFrame(rows)


def plate_layout(design, candidates):
    rows, cols = design['plate_rows'], design['plate_columns']
    if rows != 8 or cols != 12:
        raise ValueError('Este diseño implementa placas de 96 pozos.')
    rng = np.random.default_rng(design['randomization_seed'])
    wells = [f'{chr(65+r)}{c:02d}' for r in range(rows) for c in range(1, cols+1)]
    output = []
    groups = ['free_drug', 'empty_nanoparticle', 'loaded_nanoparticle']
    if design['include_physical_mixture']:
        groups.append('physical_mixture')
    for candidate in candidates.itertuples():
        for experiment in range(1, design['independent_experiments_for_variance_pilot']+1):
            for time in design['times_h']:
                treatments = []
                for group in groups:
                    for dose in design['np_doses_ug_ml']:
                        for replicate in range(1, design['technical_replicates']+1):
                            treatments.append((group, dose, replicate))
                for group, key in [('vehicle', 'vehicle_wells'), ('positive_assay_control', 'positive_control_wells'),
                                   ('medium_blank', 'cell_free_blank_wells')]:
                    treatments += [(group, None, replicate+1) for replicate in range(design[key])]
                if len(treatments) > len(wells):
                    raise ValueError('El diseño excede la capacidad; no eliminar controles para hacerlo caber.')
                treatments += [('unused_medium_only', None, None)] * (len(wells)-len(treatments))
                assignment = rng.permutation(len(wells))
                pair = f'PLAN-{candidate.drug_name}-B{experiment}-T{time}'
                for index, (group, dose, replicate) in enumerate(treatments):
                    for plate_type in ['cellular', 'cell_free_interference']:
                        has_cells = plate_type == 'cellular' and group not in ['medium_blank', 'unused_medium_only']
                        output.append({'plate_id': pair + '-' + plate_type, 'plate_pair_id': pair,
                            'well': wells[assignment[index]], 'plate_type': plate_type,
                            'compound_id': candidate.compound_id, 'drug_name': candidate.drug_name,
                            'planned_independent_experiment': experiment, 'cell_model': design['cell_model'],
                            'cells_present': has_cells, 'treatment': group,
                            'technical_replicate': replicate, 'time_h': time,
                            'matched_loaded_np_dose_ug_ml': dose,
                            'np_dose_ug_ml': dose if group in ['empty_nanoparticle', 'loaded_nanoparticle', 'physical_mixture'] else None,
                            'active_drug_conc_um': None, 'loading_pct': None,
                            'raw_readout': None, 'measurement_status': 'not_collected',
                            'include_in_measurement_plan': group != 'unused_medium_only',
                            'laboratory_approved': False})
    return pd.DataFrame(output)
