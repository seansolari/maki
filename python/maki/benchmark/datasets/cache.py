from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import json
import shutil

from .base import BenchmarkDataset
from .schema import CACHE_MANIFEST_SCHEMA_VERSION, MaterializedDataset


class DatasetCache:
    def __init__(self, root: str | Path) -> None:
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)

    def path_for(self, dataset: BenchmarkDataset) -> Path:
        return self.root / dataset.cache_key()

    def data_path_for(self, dataset: BenchmarkDataset) -> Path:
        return self.path_for(dataset) / "data"

    def manifest_path_for(self, dataset: BenchmarkDataset) -> Path:
        return self.path_for(dataset) / "manifest.json"

    def exists(self, dataset: BenchmarkDataset) -> bool:
        return self.data_path_for(dataset).exists()

    def prepare(self, dataset: BenchmarkDataset) -> Path:
        path = self.path_for(dataset)
        path.mkdir(parents=True, exist_ok=True)
        return path

    def remove(self, dataset: BenchmarkDataset) -> None:
        path = self.path_for(dataset)

        if path.exists():
            shutil.rmtree(path)

    def write_manifest(
        self,
        dataset: BenchmarkDataset,
        materialized: MaterializedDataset,
    ) -> None:
        manifest = {
            "schema_version": CACHE_MANIFEST_SCHEMA_VERSION,
            "dataset_id": dataset.dataset_id,
            "dataset_version": dataset.version,
            "cache_key": dataset.cache_key(),
            "materialized": asdict(materialized),
        }

        path = self.manifest_path_for(dataset)
        path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True),
            encoding="utf-8",
        )

    def read_manifest(self, dataset: BenchmarkDataset) -> dict:
        return json.loads(
            self.manifest_path_for(dataset).read_text(
                encoding="utf-8"
            )
        )