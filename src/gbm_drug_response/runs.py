"""Ejecuciones separadas, originales preservados y punteros a la ejecución actual."""
from __future__ import annotations

import functools
import inspect
import hashlib
import json
import os
import platform
import subprocess
import tarfile
import uuid
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path


def checksum(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    """JSON interoperable: los valores no definidos se representan como null."""
    import numpy as np
    def clean(x):
        if isinstance(x, dict):
            return {str(k): clean(v) for k, v in x.items()}
        if isinstance(x, (list, tuple, np.ndarray)):
            return [clean(v) for v in x]
        if isinstance(x, np.generic):
            return clean(x.item())
        if isinstance(x, float) and not np.isfinite(x):
            return None
        return x
    path = Path(path)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(clean(value), indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    temporary.replace(path)


def publish_link(path, target, archive):
    """Preserva archivos/carpetas anteriores; cambia el puntero de forma atómica."""
    path, target, archive = Path(path), Path(target), Path(archive)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not path.is_symlink():
        archive.mkdir(parents=True, exist_ok=True)
        path.rename(archive / (uuid.uuid4().hex + '-' + path.name))
    temporary = path.with_name('.' + path.name + '-' + uuid.uuid4().hex)
    temporary.symlink_to(os.path.relpath(target, path.parent), target_is_directory=target.is_dir())
    temporary.replace(path)


class ResearchRun:
    """Cada salida pertenece a un run; un fallo nunca recupera métricas viejas."""
    def __init__(self, root, kind):
        self.root = Path(root).resolve()
        self.kind = kind
        self.id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ') + '-' + uuid.uuid4().hex[:8]
        self.path = self.root / 'reports/runs' / kind / self.id
        self.manifest = {'run_id': self.id, 'kind': kind, 'status': 'running',
                         'started_utc': datetime.now(timezone.utc).isoformat(),
                         'python': platform.python_version(), 'inputs': {}, 'versions': {}}

    def __enter__(self):
        self.path.mkdir(parents=True)
        for package in ['numpy', 'pandas', 'scipy', 'scikit-learn', 'rdkit', 'requests', 'pyarrow']:
            try:
                self.manifest['versions'][package] = version(package)
            except PackageNotFoundError:
                pass
        try:
            self.manifest['git_commit'] = subprocess.check_output(
                ['git', '-C', str(self.root), 'rev-parse', 'HEAD'], stderr=subprocess.DEVNULL, text=True).strip()
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.manifest['git_commit'] = None
        sources = [p for pattern in ['src/**/*.py', 'scripts/*.py', 'tests/*.py',
                                     '*.json', 'pyproject.toml', 'uv.lock', 'notebooks/*.ipynb',
                                     'data/candidate_reviews.csv', 'data/formulation_measurements.csv',
                                     'docs/metodologia.md', 'docs/quimica_candidatos.md', 'docs/expresion.md']
                   for p in self.root.glob(pattern) if p.is_file()]
        self.manifest['source_sha256'] = {str(p.relative_to(self.root)): checksum(p) for p in sorted(sources)}
        with tarfile.open(self.path / 'source.tar.gz', 'w:gz') as archive:
            for p in sorted(sources):
                archive.add(p, arcname=str(p.relative_to(self.root)))
        self._save()
        publish_link(self.root / 'reports' / self.kind, self.path,
                     self.root / 'reports/runs/legacy')
        return self

    def inputs(self, paths):
        for path in paths:
            path = Path(path)
            if not path.is_file():
                raise FileNotFoundError(path)
            try:
                label = str(path.relative_to(self.root))
            except ValueError:
                label = str(path)
            self.manifest['inputs'][label] = {'sha256': checksum(path), 'bytes': path.stat().st_size}
        self._save()

    def model_link(self, filename):
        # El destino pertenece al run actual incluso cuando no se entrena un modelo.
        publish_link(self.root / 'models' / filename, self.root / 'reports' / self.kind / 'model.joblib',
                     self.root / 'reports/runs/legacy_models')

    def _save(self):
        write_json(self.path / 'run_manifest.json', self.manifest)

    def __exit__(self, exc_type, exc, tb):
        self.manifest.update(status='failed' if exc else 'complete',
                             finished_utc=datetime.now(timezone.utc).isoformat())
        if exc:
            self.manifest['error'] = {'type': exc_type.__name__, 'message': str(exc)}
        self.manifest['outputs'] = {str(p.relative_to(self.path)): checksum(p)
                                    for p in sorted(self.path.rglob('*'))
                                    if p.is_file() and p.name != 'run_manifest.json'}
        self._save()
        return False


def tracked(kind):
    def decorate(function):
        @functools.wraps(function)
        def wrapped(root=None, *args, **kwargs):
            if root is None:
                root = inspect.signature(function).parameters['root'].default
                if root is inspect.Parameter.empty:
                    raise TypeError('root is required')
            with ResearchRun(root, kind) as run:
                return function(Path(root), *args, _run=run, **kwargs)
        return wrapped
    return decorate


def verify_run(folder):
    """No publicar resultados de una corrida incompleta o alterada."""
    folder = Path(folder)
    manifest = json.loads((folder / 'run_manifest.json').read_text())
    if manifest['status'] != 'complete':
        raise ValueError(f'Ejecución no completa: {folder}')
    for name, expected in manifest['outputs'].items():
        path = folder / name
        if not path.is_file() or checksum(path) != expected:
            raise ValueError(f'Salida alterada o ausente: {path}')
    return manifest
