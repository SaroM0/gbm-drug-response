"""Genera cifras documentales desde corridas completas y verificadas."""
from pathlib import Path
import json
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from gbm_drug_response.runs import verify_run, write_json


def markdown_table(frame):
    def cell(x):
        return '' if pd.isna(x) else str(x).replace('|', '\\|').replace('\n', ' ')
    rows = ['| ' + ' | '.join(map(str, frame.columns)) + ' |',
            '| ' + ' | '.join(['---'] * len(frame.columns)) + ' |']
    rows += ['| ' + ' | '.join(cell(x) for x in row) + ' |' for row in frame.itertuples(index=False, name=None)]
    return '\n'.join(rows)


def main(root=ROOT):
    root = Path(root)
    kinds = ['chemistry', 'candidates', 'metadata_baseline', 'molecular', 'expression_preparation', 'audit']
    manifests = {k: verify_run(root / 'reports' / k) for k in kinds}
    candidates = json.loads((root / 'reports/candidates/run.json').read_text())
    chemistry = json.loads((root / 'reports/chemistry/compounds_run.json').read_text())
    molecular = json.loads((root / 'reports/molecular/status.json').read_text())
    expression = json.loads((root / 'reports/expression_preparation/expression_metadata.json').read_text())
    base = pd.read_csv(root / 'reports/metadata_baseline/test_metrics.csv')
    shortlist = pd.read_csv(root / 'reports/candidates/shortlist.csv')
    body = ['# Resultados actuales generados', '',
        'Generado por `scripts/update_results_docs.py` desde salidas con checksums verificados. '
        'Estas ejecuciones son exploratorias; el holdout ya se había observado. '
        'Ninguna cifra acredita eficacia de formulaciones o validación clínica.', '',
        '## Química y priorización', '',
        f"Se resolvieron {chemistry['resolved_drug_id_rows']} filas DRUG_ID y {chemistry['resolved_names']} nombres, "
        f"con {chemistry['unique_compound_ids']} entidades normalizadas distintas. "
        f"Hay {chemistry['rows_needing_review']} filas en el informe de revisión; un sinónimo puede ser correcto y seguir requiriendo verificación.", '',
        f"T98G: {candidates['t98g_above_range']} ajustes por encima del máximo, "
        f"{candidates['t98g_below_range']} por debajo del mínimo y {candidates['t98g_within_range']} dentro del rango.", '',
        markdown_table(pd.DataFrame(candidates['funnel'])), '',
        f"Lista corta: {candidates['shortlist_molecules']} moléculas provisionales; "
        f"revisiones químicas aprobadas: {candidates['chemistry_reviews_approved']}. "
        f"Protocolo comparativo listo: {candidates['pilot_protocol']['ready']}. "
        'Los umbrales de carga son escenarios elegidos, no capacidades experimentales universales.', '',
        markdown_table(shortlist[['rank', 'drug_name', 'tier', 'n_assays', 'ic50_min_um', 'ic50_max_um',
                                   'required_loading_pct']].round(6)), '',
        'La potencia representativa usa el máximo LN_IC50 entre todos los ensayos, conservando IDs, '
        'discrepancias y ambos extremos de rango. Ver `reports/candidates/assays.csv`, '
        '`dose_scenarios.csv`, `scenario_summary.csv` y `formulation_readiness.csv`.', '',
        '## Metadatos', '', markdown_table(base.round(6)), '',
        '## Expresión real', '',
        f"Matriz preparada: {expression['n_matched_lines']} líneas × {expression['n_genes']} genes. "
        f"Líneas sin perfil: {', '.join(expression['unmatched_lines'])}. "
        'Se verificó escala TPM cruda mediante sus sumas por perfil y se aplicó log2(TPM+1) una vez. '
        'Se conservaron los ausentes; no se imputó antes de CV.', '',
        f"Estado molecular: `{molecular['status']}`. "
        'Los fármacos se eligieron por cobertura en desarrollo, no por error de test.', '']
    if molecular['status'] == 'trained':
        scores = pd.read_csv(root / 'reports/molecular/test_metrics.csv', dtype={'DRUG_ID': str})
        selected = scores[scores.selected_by_cv].copy()
        baseline = scores[scores.model == 'drug_mean'][['DRUG_ID', 'mae']].rename(columns={'mae': 'baseline_mae'})
        selected = selected.merge(baseline, on='DRUG_ID', validate='one_to_one')
        intervals = pd.read_csv(root / 'reports/molecular/paired_comparisons.csv', dtype={'DRUG_ID': str})
        selected = selected.merge(intervals[['DRUG_ID', 'model', 'delta_ci95_low', 'delta_ci95_high']],
                                  on=['DRUG_ID', 'model'], validate='one_to_one')
        body += [markdown_table(selected[['DRUG_ID', 'DRUG_NAME', 'model', 'n', 'mae', 'baseline_mae',
                                         'delta_ci95_low', 'delta_ci95_high']].round(6)), '',
                 'Se muestran los modelos elegidos por CV, aunque otro tenga menor error de test. '
                 'Los intervalos bootstrap son condicionales al entrenamiento fijo; tres fármacos '
                 'y once líneas de test por fármaco no constituyen validación general.', '']
    body += ['Los desgloses `outcome_diagnostics.csv` informan métricas con el mismo modelo '
             'en ajustes dentro/por encima/por debajo de rango. No se selecciona modelo con esos subgrupos '
             'ni se convierten extrapolaciones en cotas de censura medidas.', '', '## Ejecuciones', '',
             markdown_table(pd.DataFrame([{'kind': k, 'run_id': v['run_id'], 'status': v['status']}
                                         for k, v in manifests.items()])), '']
    (root / 'docs/resultados_actuales.md').write_text('\n'.join(body))
    write_json(root / 'docs/evidencia_correcciones_2026-10-29.json', {
        'scope': 'Software correction and exploratory rerun; laboratory validation remains pending',
        'runs': manifests, 'candidates': candidates, 'chemistry': chemistry,
        'expression': expression, 'molecular_status': molecular})
    print('Documentación y evidencia generadas desde ejecuciones verificadas.')


if __name__ == '__main__':
    main()
