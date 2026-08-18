from __future__ import annotations

import os
from pathlib import Path
import platform
import socket
import subprocess
from datetime import datetime

from .schema import BenchmarkMetadata


def _git_value(args: list[str]) -> str | None:
    try:
        return subprocess.check_output(
            args,
            text=True,
            cwd=Path(__file__).parent
        ).strip()
    except Exception:
        return None


def collect_metadata(
    threads: int,
) -> BenchmarkMetadata:

    return BenchmarkMetadata(
        git_commit=_git_value(
            ["git", "rev-parse", "HEAD"]
        ),
        git_branch=_git_value(
            ["git", "branch", "--show-current"]
        ),
        python_version=platform.python_version(),
        hostname=socket.gethostname(),
        cpu_count=os.cpu_count() or 1,
        thread_count=threads,
        timestamp_utc=datetime.utcnow().isoformat()
    )