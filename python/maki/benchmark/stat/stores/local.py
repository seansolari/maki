from __future__ import annotations

import json
from pathlib import Path

from ..schema import WorkflowBenchmarkResult
from .base import ResultStore


class LocalResultStore(ResultStore):
    def __init__(self, output_dir: str | Path):
        self.output_dir = Path(output_dir)

    def save(self, result: WorkflowBenchmarkResult) -> None:
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        file_path = self.output_dir / f"{result.run_id}.json"
        
        with file_path.open("w", encoding="utf-8") as handle:
            json.dump(result.to_json(), handle, indent=2)
