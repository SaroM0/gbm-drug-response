"""Informe y figuras del avance: datos crudos disponibles y piloto aún propuesto."""
from pathlib import Path
import json
import shutil
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from gbm_drug_response.runs import ResearchRun, verify_run, write_json
from gbm_drug_response.pilot import readiness


def table(frame):
    def cell(x):
        return '' if pd.isna(x) else str(x).replace('|', '\\|')
    return '\n'.join(['| '+' | '.join(frame.columns)+' |', '| '+' | '.join(['---']*len(frame.columns))+' |'] +
        ['| '+' | '.join(cell(x) for x in row)+' |' for row in frame.itertuples(index=False, name=None)])


def main(root=ROOT):
    root = Path(root)
    kinds = ['gdsc_source_acquisition', 'gdsc_raw_preparation', 'pilot_preparation', 'candidates']
    verified = {kind: verify_run(root / 'reports' / kind) for kind in kinds}
    with ResearchRun(root, 'pending_progress') as run:
        paths = [root / 'reports/gdsc_raw_preparation' / name for name in
            ['run.json', 't98g_empirical_curves.csv', 'glioma_normalized_wells.parquet']]
        paths += [root / 'reports/pilot_preparation/run.json', root / 'reports/pilot_preparation/sample_size_scenarios.csv',
                  root / 'reports/candidates/ranking.csv', root / 'pilot_protocol.json', root / 'data/raw/Aislados T98G - Hoja 1.csv']
        run.inputs(paths)
        raw = json.loads(paths[0].read_text())
        series = pd.read_csv(paths[1], dtype={'DRUG_ID': str})
        wells = pd.read_parquet(paths[2])
        ranking = pd.read_csv(root / 'reports/candidates/ranking.csv')
        fitted = pd.read_csv(root / 'data/raw/Aislados T98G - Hoja 1.csv', dtype={'DRUG_ID': str})
        pilot = json.loads((root / 'reports/pilot_preparation/run.json').read_text())
        protocol = readiness(root)
        missing = fitted[~fitted.DRUG_ID.isin(series.DRUG_ID)][['DRUG_ID', 'DRUG_NAME']]
        missing.to_csv(run.path / 't98g_fitted_drugs_without_usable_raw.csv', index=False)
        rows = []
        names = ['Docetaxel', 'Gemcitabine', 'OSI-027']
        fig, axes = plt.subplots(1, 3, figsize=(12, 4.4), sharey=True)
        for ax, name in zip(axes, names):
            candidate = ranking[ranking.drug_name.eq(name)].iloc[0]
            ids = candidate.all_gdsc_drug_ids.split(';')
            part = series[series.DRUG_ID.isin(ids)]
            for drug in ids:
                d = part[part.DRUG_ID.eq(drug)]
                rows.append({'candidate': name, 'DRUG_ID': drug, 'plate_drug_series': len(d),
                    'series_with_points_on_both_sides_of_half': int(d.empirical_range_state.eq('points_on_both_sides_of_half').sum()),
                    'series_all_above_half': int(d.empirical_range_state.eq('all_points_above_half').sum()),
                    'raw_max_conc_um_by_series': ';'.join(map(str, sorted(d.raw_max_conc.unique())))})
            selected = wells[wells.SANGER_MODEL_ID.eq('SIDM01171') & wells.DRUG_ID.isin(ids)
                & wells.is_library_single & ~wells.unsupported_combination & np.isfinite(wells.normalized_viability)
                & np.isfinite(wells.CONC) & wells.CONC.gt(0)]
            for (drug, barcode), points in selected.groupby(['DRUG_ID', 'BARCODE'], sort=True):
                median = points.groupby('CONC').normalized_viability.median().sort_index()
                ax.plot(median.index, median.values, '.-', linewidth=1, markersize=4, alpha=.8,
                        label=f'{drug} · placa {barcode}')
            ax.axhline(.5, color='#333333', linestyle='--', linewidth=.8)
            ax.set_xscale('log'); ax.set_title(name); ax.set_xlabel('Concentración libre (µM)')
            ax.set_ylim(-.15, 1.45); ax.grid(alpha=.15)
            ax.legend(fontsize=6.5, frameon=False, loc='lower left')
        axes[0].set_ylabel('Readout relativo a controles NC-1 y B')
        fig.suptitle('T98G · medianas por dosis y placa del crudo público', fontsize=12)
        fig.text(.5, .01, 'Diagnóstico descriptivo; placas no equivalen a campañas biológicas independientes. No es MTT ni reajuste NLME.',
                 ha='center', fontsize=8)
        fig.tight_layout(rect=[0, .04, 1, .95])
        fig.savefig(run.path / 'candidate_raw_curves.png', dpi=180)
        fig.savefig(run.path / 'candidate_raw_curves.svg')
        plt.close(fig)
        evidence = pd.DataFrame(rows)
        evidence.to_csv(run.path / 'candidate_raw_support.csv', index=False)
        (root / 'docs/figures').mkdir(exist_ok=True)
        shutil.copyfile(run.path / 'candidate_raw_curves.png', root / 'docs/figures/candidate_raw_curves.png')
        power = pd.read_csv(root / 'reports/pilot_preparation/sample_size_scenarios.csv')
        report = ['# Avance de los pendientes', '',
            'Se adquirieron fuentes oficiales y se preparó un piloto revisable. Las mediciones de las nuevas formulaciones siguen pendientes; no se fabricaron ni se aprobaron candidatos automáticamente.', '',
            '## Datos crudos y procedencia', '',
            f"Se descargó el ZIP público de GDSC2 8.5 y se leyeron {raw['raw_annotation_rows_scanned']:,} anotaciones sin extraer permanentemente el CSV de 2,11 GB. "
            f"Se conservaron {raw['selected_annotation_rows']:,} anotaciones, que corresponden a {raw['selected_physical_wells']:,} pozos físicos de {raw['matched_lines']} líneas y {raw['selected_plates']} placas. "
            'El CRC del miembro ZIP se verificó por lectura completa; cada archivo tiene recibo HTTP y SHA-256.', '',
            f"Los {raw['fitted_comparison']['matching_response_curves']:,} ajustes locales coinciden con la distribución oficial en claves y cinco campos de respuesta, con tolerancia absoluta de 1e-12. "
            'Esto verifica la correspondencia de respuestas, no igualdad de bytes del Excel ni autenticidad de todas las anotaciones tumorales.', '',
            'No se sumaron intensidades repetidas por etiquetas. El tag DMSO puede llevar DRUG_ID=4000: se conserva como solvente y no como un segundo fármaco. '
            'La normalización descriptiva usa sólo controles NC-1 y B explícitos; PC-1 y UN-USED no se reinterpretan como blancos. '
            f"Las {raw['plates_with_explicit_usable_nc1_b']} placas tuvieron controles suficientes bajo esta regla operativa. "
            'No se sustituyó el QC original ni se reprodujo el modelo NLME. [Documentación primaria del formato](https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC_Raw_Data_Description.pdf).', '',
            '## Qué muestran las curvas de los candidatos', '', table(evidence), '',
            '![Curvas crudas de T98G](figures/candidate_raw_curves.png)', '',
            'Hay variación entre placas y rangos de concentración. Para gemcitabina aparecen máximos de 0,1, 1 y 10 µM; el rango agregado del ajuste no debe atribuirse a cada experimento. '
            'Tener puntos a ambos lados de 50% es un diagnóstico empírico; no garantiza monotonicidad, precisión de IC50 ni una cota de censura calibrada. '
            'Las diferencias entre placas no se cuentan como réplicas biológicas independientes sin verificar expansión/lote.', '',
            f"En T98G se obtuvieron {len(series)} series placa–fármaco normalizables. No todos los ajustes tienen una serie pública utilizable bajo estas reglas; "
            'se conserva la lista de ausencias sin inventar puntos ni sustituirlas por resistencia:', '', table(missing), '',
            '## Piloto y tamaño muestral', '',
            f"Se preparó una propuesta de {pilot['n_plates']} placas, con posiciones aleatorizadas, blancos, controles y una mezcla física adicional para discutir. "
            'Usa tres campañas para estimación inicial de variabilidad y tres pozos técnicos por condición; no afirma potencia confirmatoria. '
            'Las lecturas y dosis de fármaco dependientes de carga siguen vacías. Ver [propuesta completa para el laboratorio](piloto_para_revision.md).', '',
            'Escenarios ilustrativos para una reducción de IC50 del 30%, con alpha familiar 0,05/3 y contraste pareado bilateral:', '',
            table(power[power.loaded_to_free_ic50_ratio.eq(.7)].round(5)), '',
            'La SD de cada fila es hipotética. El tamaño confirmatorio no está fijado: requiere efecto relevante acordado y variabilidad medida en campañas independientes.', '',
            '## Lo que sigue abierto', '',
            '- Revisar identidad, forma comercial, ruta, sitio de unión y actividad/liberación de cada formulación.',
            '- Medir carga/liberación por lote, controles MTT y readout ortogonal; incorporar el antecedente de oxaliplatino con sus datos instrumentales.',
            '- Fijar efecto relevante, número confirmatorio y criterios de endotelio, barrera y recuperación en Transwell.',
            '- Reajustar o modelar incertidumbre/censura con un protocolo justificado a partir de estos crudos; no convertir directamente las categorías descriptivas en etiquetas ciertas.', '',
            f"Listo para comparación formal: `{protocol['ready']}`; sigue sin aprobación. Decisiones aún sin valor: " + ', '.join(protocol['missing_decisions']) + '.', '',
            '## Reproducir', '', '```bash', '.venv/bin/python scripts/acquire_gdsc_sources.py',
            '.venv/bin/python scripts/prepare_gdsc_raw.py', '.venv/bin/python scripts/prepare_pilot.py',
            '.venv/bin/python scripts/report_pending_progress.py', '```', '',
            'Todos los análisis y cambios de protocolo posteriores al holdout observado siguen siendo exploratorios.', '']
        (root / 'docs/avance_pendientes.md').write_text('\n'.join(report))
        summary = {'run_id': run.id, 'verified_source_runs': verified, 'raw': raw,
            'candidate_support': rows, 'pilot': pilot, 'protocol_readiness': protocol,
            'laboratory_observations_invented': 0, 'nlme_refit_performed': False}
        write_json(run.path / 'run.json', summary)
        write_json(root / 'docs/evidencia_avance_pendientes.json', summary)
        print('Informe, figura y evidencia del avance generados.')
        return summary


if __name__ == '__main__':
    main()
