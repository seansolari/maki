
"""
Pytest configuration and shared fixtures.
"""
import os
import shutil
import subprocess
import pytest

# Global timeout for pipeline runs (seconds) — change via env if needed
GLOBAL_TIMEOUT = int(os.getenv("MAKI_TEST_TIMEOUT", "30"))
DEFAULT_PROFILE = os.getenv("MAKI_TEST_PROFILE", "conda")  # Preferred Nextflow profile

@pytest.fixture(scope="session")
def rng():
    """Return a reproducible random generator (NumPy if available, else stdlib).

    Seed is controlled by env var ``MAKI_TEST_SEED`` (default: 42).
    """
    seed = int(os.getenv("MAKI_TEST_SEED", "42"))
    try:
        import numpy as np
        return np.random.default_rng(seed)
    except Exception:
        import random
        return random.Random(seed)

@pytest.fixture(scope="session")
def nf_runner():
    """Check Nextflow and maki availability; provide a runner callable.

    Uses environment variables:
      - NEXTFLOW_PATH (default: 'nextflow')
      - MAKI_PATH     (default: 'maki')
      - MAKI_PIPELINE (default: 'maki_nf/main.nf')
    """
    nf_path = os.getenv("NEXTFLOW_PATH", "nextflow")
    maki_path = os.getenv("MAKI_PATH", "maki")
    pipeline_path = os.getenv("MAKI_PIPELINE", "maki_nf/main.nf")

    if not shutil.which(nf_path):
        pytest.skip("Nextflow not found; set NEXTFLOW_PATH or install Nextflow")
    if not shutil.which(maki_path):
        pytest.skip("maki executable not found; set MAKI_PATH or add to PATH")

    def _run(entry: str, params_file: str, workdir: str, profile: str = DEFAULT_PROFILE):
        cmd = [nf_path, "run", pipeline_path, "-entry", entry, "-params-file", params_file]
        if profile:
            cmd += ["-profile", profile]
        return subprocess.run(
            cmd, cwd=workdir, capture_output=True, text=True, timeout=GLOBAL_TIMEOUT
        )

    return _run
