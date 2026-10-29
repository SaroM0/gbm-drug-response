"""Adquiere fuentes oficiales de GDSC 8.5 sin sobrescribir los adjuntos."""
from pathlib import Path
import json
import sys
import xml.etree.ElementTree as ET

import requests

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from gbm_drug_response.runs import ResearchRun, checksum, write_json

BASE = 'https://cog.sanger.ac.uk/cancerrxgene/'
PREFIX = 'GDSC_release8.5/'
FILES = ['GDSC2_public_raw_data_27Oct23.csv.zip', 'GDSC2_fitted_dose_response_27Oct23.xlsx',
         'screened_compounds_rel_8.5.csv', 'GDSC_Raw_Data_Description.pdf']
MAX_BYTES = 256 * 1024 * 1024


def catalog(xml):
    document = ET.fromstring(xml)
    ns = {'s': 'http://s3.amazonaws.com/doc/2006-03-01/'}
    if document.findtext('s:IsTruncated', namespaces=ns) != 'false':
        raise ValueError('Índice incompleto; no declarar ausencia de una fuente.')
    return {e.findtext('s:Key', namespaces=ns): {
                'bytes': int(e.findtext('s:Size', namespaces=ns)),
                'etag': e.findtext('s:ETag', namespaces=ns)}
            for e in document.findall('s:Contents', ns)}


def download(url, target, expected_size):
    target = Path(target)
    if expected_size <= 0 or expected_size > MAX_BYTES:
        raise ValueError('Tamaño fuera del límite de adquisición.')
    receipt = target.with_suffix(target.suffix + '.receipt.json')
    if target.exists() and receipt.exists():
        meta = json.loads(receipt.read_text())
        if (meta['url'] == url and target.stat().st_size == expected_size
                and checksum(target) == meta['sha256']):
            return {**meta, 'cache_reused': True}
        raise ValueError(f'Caché alterada: {target}; conservar e investigar antes de reemplazar.')
    if target.exists():
        raise ValueError(f'Archivo sin recibo verificable: {target}')
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + '.partial')
    received = 0
    with requests.get(url, stream=True, timeout=(15, 60)) as response:
        response.raise_for_status()
        if not response.url.startswith(BASE):
            raise ValueError('Redirección fuera de la fuente oficial esperada.')
        with partial.open('wb') as stream:
            for block in response.iter_content(1024 * 1024):
                received += len(block)
                if received > expected_size or received > MAX_BYTES:
                    raise ValueError('Descarga excede el tamaño publicado.')
                stream.write(block)
        if received != expected_size:
            raise ValueError('Descarga incompleta respecto al índice oficial.')
        meta = {'url': url, 'final_url': response.url, 'bytes': received,
                'sha256': checksum(partial), 'etag': response.headers.get('ETag'),
                'last_modified': response.headers.get('Last-Modified'), 'cache_reused': False}
    partial.replace(target)
    write_json(receipt, meta)
    return meta


def main(root=ROOT):
    root = Path(root)
    with ResearchRun(root, 'gdsc_source_acquisition') as run:
        response = requests.get(BASE, params={'prefix': PREFIX}, timeout=(15, 30))
        response.raise_for_status()
        (run.path / 'official_catalog.xml').write_text(response.text)
        entries = catalog(response.text)
        records = []
        for name in FILES:
            key = PREFIX + name
            if key not in entries:
                raise ValueError(f'Fuente esperada ausente del índice: {key}')
            target = root / 'data/external/gdsc_release8_5' / name
            print(f'Adquiriendo {name} ({entries[key]["bytes"]} bytes)', flush=True)
            meta = download(BASE + key, target, entries[key]['bytes'])
            meta['local_path'] = str(target.relative_to(root))
            records.append(meta)
            run.inputs([target, target.with_suffix(target.suffix + '.receipt.json')])
        write_json(run.path / 'sources.json', {'release': '8.5', 'files': records,
            'raw_vs_fitted_linkage': 'not_yet_verified', 'user_originals_overwritten': False})
        print('Fuentes oficiales adquiridas y checksums registrados.', flush=True)
        return records


if __name__ == '__main__':
    main()
