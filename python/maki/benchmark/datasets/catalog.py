from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import json

from .base import BenchmarkDataset
from .registry import DatasetRegistry


CATALOG_SCHEMA_VERSION = 1


class DatasetCatalog:
    def __init__(self, registry: DatasetRegistry):
        self.registry = registry

    def to_dict(self) -> dict:
        return {
            "schema_version":
                CATALOG_SCHEMA_VERSION,
            "datasets": [
                asdict(dataset)
                for dataset
                in self.registry.datasets.values()
            ],
        }

    def save(self, path: Path) -> None:
        with path.open("w", encoding="utf-8") as handle:
            json.dump(self.to_dict(), handle, indent=2)

    @classmethod
    def load(cls, path: Path) -> DatasetRegistry:
        with path.open(encoding="utf-8") as handle:
            payload = json.load(handle)

        registry = DatasetRegistry()

        for item in payload["datasets"]:
            dataset = BenchmarkDataset(**item)
            registry.register(dataset)

        return registry