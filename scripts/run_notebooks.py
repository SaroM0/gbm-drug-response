"""Ejecuta notebooks con el Python activo y conserva outputs incluso en caso de error."""
from pathlib import Path
import sys
import tempfile

import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager
from jupyter_client.kernelspec import KernelSpecManager
from ipykernel.kernelspec import write_kernel_spec

ROOT = Path(__file__).resolve().parents[1]


def main():
    names = sys.argv[1:] or [p.name for p in sorted((ROOT / "notebooks").glob("*.ipynb"))]
    # Kernel privado y temporal: no instala ni altera kernels del usuario.
    with tempfile.TemporaryDirectory(prefix="gbm-kernel-") as kernel_dir:
        write_kernel_spec(path=str(Path(kernel_dir) / "python3"))
        specs = KernelSpecManager(kernel_dirs=[kernel_dir])
        for name in names:
            path = ROOT / "notebooks" / name
            notebook = nbformat.read(path, as_version=4)
            nbformat.validate(notebook)
            manager = KernelManager(kernel_name="python3", kernel_spec_manager=specs)
            client = NotebookClient(notebook, km=manager, timeout=900,
                                    resources={"metadata": {"path": str(ROOT)}})
            print(f"Ejecutando {path.name}...", flush=True)
            try:
                client.execute()
            finally:
                if manager.has_kernel:
                    manager.shutdown_kernel(now=True)
                if client.kc is not None:
                    client.kc.stop_channels()
                nbformat.write(notebook, path)
            print(f"Completado: {path.name}", flush=True)


if __name__ == "__main__":
    main()
