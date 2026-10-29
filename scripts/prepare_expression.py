"""Prepara TPM del ZIP local; valida formato/escala y conserva perfiles y símbolos."""
from pathlib import Path
import csv
import io
import sys
import zipfile

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from gbm_drug_response.data import FILES, read_responses
from gbm_drug_response.runs import ResearchRun, publish_link, write_json, checksum


def read_tpm_member(stream, wanted):
    reader = csv.reader(io.TextIOWrapper(stream, encoding='utf-8'))
    header = next(reader)
    if header[:2] != ['model_id', '']:
        raise ValueError('Cabecera TPM histórica inesperada.')
    if len(header[2:]) != len(set(header[2:])):
        raise ValueError('Perfiles SIDM repetidos; elegir perfil requiere protocolo explícito.')
    metadata = [next(reader) for _ in range(4)]
    if [r[0] for r in metadata] != ['model_name', 'dataset_name', 'data_source', 'gene_id']:
        raise ValueError('Filas de metadatos históricas inesperadas.')
    matched = sorted(set(wanted) & set(header[2:]))
    if not matched:
        raise ValueError('Ningún SIDM solicitado tiene expresión.')
    positions = [header.index(x) for x in matched]
    genes, symbols, values, seen = [], [], [], set()
    for row in reader:
        if len(row) != len(header) or not row[0].startswith('SIDG') or row[0] in seen:
            raise ValueError('Genes duplicados, ID inválido o fila de longitud inesperada.')
        seen.add(row[0])
        v = np.array([float(row[p]) if row[p] else np.nan for p in positions])
        if np.isinf(v).any() or (v < 0).any():
            raise ValueError('TPM infinito o negativo; los ausentes se conservan, no se rellenan.')
        genes.append('gene_' + row[0]); symbols.append(row[1]); values.append(v)
    matrix = np.array(values, float).T
    # Una distribución ya log-transformada no suma ~1e6. No inferir escala sólo por nombre.
    sums = np.nansum(matrix, axis=1)
    if len(genes) < 2 or np.isnan(matrix).all(axis=1).any() or not np.allclose(sums, 1e6, rtol=.02, atol=0):
        raise ValueError('No se valida escala TPM cruda: sumas por perfil fuera de 1e6 ±2%. No aplicar otro log.')
    expr = pd.DataFrame(np.log2(matrix + 1), columns=genes)
    expr.insert(0, 'SANGER_MODEL_ID', matched)
    profiles = pd.DataFrame({'SANGER_MODEL_ID': matched,
        'model_name': [metadata[0][p] for p in positions],
        'dataset_name': [metadata[1][p] for p in positions],
        'data_source': [metadata[2][p] for p in positions], 'raw_tpm_sum': sums,
        'n_missing_genes': np.isnan(matrix).sum(axis=1)})
    annotation = pd.DataFrame({'feature': genes, 'gene_symbol': symbols})
    return expr, profiles, annotation, sorted(set(wanted) - set(matched))


def main(root=ROOT):
    root = Path(root)
    source = root / 'data/external/rnaseq_all_20220624.zip'
    with ResearchRun(root, 'expression_preparation') as run:
        for name in ['expression.csv', 'expression_metadata.json']:
            publish_link(root / 'data/external' / name, root / 'reports/expression_preparation' / name,
                         root / 'reports/runs/legacy_expression')
        run.inputs([source, root / 'data/raw' / FILES['glioma_csv']])
        wanted = read_responses(root / 'data/raw' / FILES['glioma_csv']).SANGER_MODEL_ID.unique()
        with zipfile.ZipFile(source) as z:
            with z.open('rnaseq_tpm_20220624.csv') as stream:
                expr, profiles, annotation, unmatched = read_tpm_member(stream, wanted)
        # Lectura completa del miembro TPM verifica también su CRC, sin extraer 5 GB largos.
        expr.to_csv(run.path / 'expression.csv', index=False)
        profiles.to_csv(run.path / 'profiles.csv', index=False)
        annotation.to_csv(run.path / 'gene_annotation.csv', index=False)
        pd.DataFrame({'SANGER_MODEL_ID': unmatched}).to_csv(run.path / 'unmatched_lines.csv', index=False)
        meta = {'source_url': 'https://cog.sanger.ac.uk/cmp/download/rnaseq_all_20220624.zip',
            'source_authenticity': 'local_file_official_byte_identity_not_verified',
            'release': '20220624', 'scale': 'log2_tpm_plus1',
            'gene_id_type': 'Sanger SIDG gene identifiers; symbols retained separately',
            'preprocessing': 'Merged wide TPM member; unique SIDM profiles as supplied; select response SIDMs; validate nonnegative TPM, no infinities and observed profile sums within 2% of 1e6; log2(TPM+1) once; transpose. Preserve partial missing values; no outcome-based filtering, cross-sample scaling or imputation.',
            'input_scale': 'raw_tpm', 'input_scale_evidence': 'Historical TPM member with unlogged abundances and profile sums approximately 1e6; current documentation concerns newer releases.',
            'sha256': checksum(run.path / 'expression.csv'), 'zip_sha256': checksum(source),
            'member': 'rnaseq_tpm_20220624.csv', 'run_id': run.id,
            'n_matched_lines': len(expr), 'n_genes': expr.shape[1]-1, 'unmatched_lines': unmatched,
            'profile_data_sources': profiles.data_source.value_counts().to_dict()}
        write_json(run.path / 'expression_metadata.json', meta)
        print(f'Expresión preparada: {len(expr)} líneas, {expr.shape[1]-1} genes; no emparejadas: {unmatched}')
        return meta


if __name__ == '__main__':
    main()
