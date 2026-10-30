"""Adquiere fuentes oficiales de GDSC 8.5 sin sobrescribir los adjuntos."""
from pathlib import Path
import json
import sys
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor

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


def download(url, target, expected_size, expected_etag=None):
    target = Path(target)
    if expected_size <= 0 or expected_size > MAX_BYTES:
        raise ValueError('Tamaño fuera del límite de adquisición.')
    receipt = target.with_suffix(target.suffix + '.receipt.json')
    if target.exists() and receipt.exists():
        meta = json.loads(receipt.read_text())
        if (meta['url'] == url and target.stat().st_size == expected_size
                and (expected_etag is None or meta.get('etag') == expected_etag)
                and checksum(target) == meta['sha256']):
            return {**meta, 'cache_reused': True}
        raise ValueError(f'Caché alterada: {target}; conservar e investigar antes de reemplazar.')
    if target.exists():
        raise ValueError(f'Archivo sin recibo verificable: {target}')
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + '.parallel.partial')
    if expected_size > 8 * 1024 * 1024 and expected_etag:
        block_size = 4 * 1024 * 1024
        ranges = [(start, min(start+block_size-1, expected_size-1))
                  for start in range(0, expected_size, block_size)]
        with partial.open('wb') as stream:
            stream.truncate(expected_size)
        def fetch_byte_range(bounds):
            start, end = bounds
            headers = {'Range': f'bytes={start}-{end}', 'If-Match': expected_etag, 'Accept-Encoding': 'identity'}
            with requests.get(url, headers=headers, stream=True, timeout=(15, 60)) as response:
                response.raise_for_status()
                wanted_range = f'bytes {start}-{end}/{expected_size}'
                if (response.status_code != 206 or response.headers.get('Content-Range') != wanted_range
                        or response.headers.get('ETag') != expected_etag or not response.url.startswith(BASE)):
                    raise ValueError('El servidor no preservó rango, versión o procedencia.')
                count = 0
                with partial.open('r+b') as stream:
                    stream.seek(start)
                    for block in response.iter_content(256 * 1024):
                        count += len(block)
                        if count > end-start+1:
                            raise ValueError('Rango descargado excedido.')
                        stream.write(block)
                if count != end-start+1:
                    raise ValueError('Rango descargado incompleto.')
        with ThreadPoolExecutor(max_workers=4) as pool:
            for index, _ in enumerate(pool.map(fetch_byte_range, ranges), 1):
                print(f'  {target.name}: {index}/{len(ranges)} bloques verificados', flush=True)
        meta = {'url': url, 'final_url': url, 'bytes': expected_size, 'sha256': checksum(partial),
                'etag': expected_etag, 'cache_reused': False,
                'transport': '4 bounded range workers, matching ETag and Content-Range'}
        partial.replace(target)
        write_json(receipt, meta)
        return meta
    received = 0
    with requests.get(url, stream=True, timeout=(15, 60), headers={'Accept-Encoding': 'identity'}) as response:
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
            meta = download(BASE + key, target, entries[key]['bytes'], entries[key]['etag'])
            meta['local_path'] = str(target.relative_to(root))
            records.append(meta)
            run.inputs([target, target.with_suffix(target.suffix + '.receipt.json')])
        write_json(run.path / 'sources.json', {'release': '8.5', 'files': records,
            'raw_vs_fitted_linkage': 'not_yet_verified', 'user_originals_overwritten': False})
        print('Fuentes oficiales adquiridas y checksums registrados.', flush=True)
        return records


if __name__ == '__main__':
    main()
