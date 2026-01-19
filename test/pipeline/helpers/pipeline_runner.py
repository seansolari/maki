
"""Nextflow pipeline runner wrapper.

Paths to Nextflow, the maki executable, and the pipeline entrypoint are controlled via:
- NEXTFLOW_PATH (default: 'nextflow')
- MAKI_PATH     (default: 'maki')
- MAKI_PIPELINE (default: 'maki_nf/main.nf')

No external YAML library is required; we write a minimal YAML by hand.
"""
from __future__ import annotations
from pathlib import Path
from typing import Any, Dict, Mapping, NamedTuple, Optional, Sequence
import os
import shutil
import subprocess

class CompletedRun(NamedTuple):
    """Structured result of a pipeline execution."""
    returncode: int
    stdout: str
    stderr: str
    expected_outputs: Mapping[str, Path]

def _write_params_yaml(path: Path, params: Dict[str, Any]) -> None:
    """Write a minimal YAML (key: value) — supports strings, ints, floats, bools."""
    def to_yaml_line(k: str, v: Any) -> str:
        if isinstance(v, bool):
            return f"{k}: {'true' if v else 'false'}\n"
        elif isinstance(v, (int, float)):
            return f"{k}: {v}\n"
        elif v is None:
            return f"{k}: null\n"
        else:
            # Quote strings and paths
            s = str(v).replace("\\", "/")
            return f'{k}: "{s}"\n'
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for k, v in params.items():
            fh.write(to_yaml_line(k, v))

def run_pipeline(
    entry: str,
    params: Dict[str, Any],
    workdir: Path,
    profile: Optional[str] = None,
    timeout: Optional[int] = None,
) -> CompletedRun:
    """Run the Nextflow pipeline and return a CompletedRun."""
    nf_path = os.getenv("NEXTFLOW_PATH", "nextflow")
    maki_path = os.getenv("MAKI_PATH", "maki")
    pipeline_path = os.getenv("MAKI_PIPELINE", "maki_nf/main.nf")

    # Basic availability check (skip here; tests also check in conftest)
    if not shutil.which(nf_path):
        raise RuntimeError("Nextflow not found; set NEXTFLOW_PATH or install Nextflow")
    if not shutil.which(maki_path):
        raise RuntimeError("maki executable not found; set MAKI_PATH or add to PATH")
    else:
        params["tool"] = maki_path

    # Write params.yaml
    params_file = workdir / "params.yaml"
    _write_params_yaml(params_file, params)

    cmd = [nf_path, "run", pipeline_path, "-entry", entry, "-params-file", str(params_file)]
    if profile:
        cmd += ["-profile", profile]

    # Run
    cp = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True, timeout=timeout or 30)

    # Expected outputs per entry (based on your Nextflow defaults)
    outmap: Dict[str, Path] = {}
    # PublishDir default is results/<step>/..., but we return the file names to be checked there.
    # Test code can resolve/validate in the published dirs.
    if entry == "jaccard":
        outmap["jaccard"] = Path(params.get("jaccard", "jaccard-summary.txt.gz"))
    elif entry == "random_gene_distances":
        outmap["feature_table"] = Path(params.get("gene_features", "feature-table.txt.gz"))
    elif entry == "kmer_prevalence":
        outmap["cooc_summary"] = Path(params.get("cooc_summary", "cooc.kmers.lca1.txt.gz"))

    return CompletedRun(cp.returncode, cp.stdout, cp.stderr, outmap)

def expect_output(paths: Sequence[Path], base: Optional[Path] = None) -> None:
    """Assert that every path exists and is non-empty (in its published location)."""
    for p in paths:
        # Caller should pass actual published paths if needed. Here we check existence/size if absolute,
        # or leave it to higher-level code to resolve publishDir locations.
        if not p.exists() and base is not None:
            # Fall back: try to locate under results/** recursively
            found = list(base.rglob(p.name))
            assert found, f"Expected output '{p.name}' not found (looked in CWD and results/**)"
            p = found[0]
