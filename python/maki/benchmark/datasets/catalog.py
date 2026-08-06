from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Any
import json

from .base import (
    BenchmarkDataset,
    DerivedDataset,
    LocalDataset,
    RealDataset,
    SyntheticDataset,
)
from .registry import DatasetRegistry
from .schema import CATALOG_SCHEMA_VERSION, DatasetFormat, DatasetMetadata, DatasetType


_DATASET_CLASSES: dict[str, type[BenchmarkDataset]] = {
    DatasetType.REAL.value: RealDataset,
    DatasetType.LOCAL.value: LocalDataset,
    DatasetType.SYNTHETIC.value: SyntheticDataset,
    DatasetType.DERIVED.value: DerivedDataset,
}


class DatasetCatalog:
    def __init__(
        self,
        registry: DatasetRegistry | None = None,
        schema_version: int = CATALOG_SCHEMA_VERSION,
        catalog_version: str = "1",
    ) -> None:
        self.registry = registry or DatasetRegistry()
        self.schema_version = schema_version
        self.catalog_version = catalog_version

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "catalog_version": self.catalog_version,
            "datasets": [
                _dataset_to_dict(dataset)
                for dataset in self.registry.values()
            ],
        }

    def save(self, path: str | Path) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)

        path.write_text(
            json.dumps(
                self.to_dict(),
                indent=2,
                sort_keys=True,
            ),
            encoding="utf-8",
        )

    @classmethod
    def load(cls, path: str | Path) -> DatasetCatalog:
        path = Path(path)
        payload = json.loads(path.read_text(encoding="utf-8"))

        schema_version = payload["schema_version"]

        if schema_version != CATALOG_SCHEMA_VERSION:
            raise ValueError(
                f"Unsupported catalog schema version: {schema_version}"
            )

        registry = DatasetRegistry()

        for item in payload["datasets"]:
            registry.register(_dataset_from_dict(item))

        return cls(
            registry=registry,
            schema_version=schema_version,
            catalog_version=payload.get("catalog_version", "1"),
        )


def _dataset_to_dict(dataset: BenchmarkDataset) -> dict[str, Any]:
    payload = asdict(dataset)

    payload["dataset_type"] = dataset.dataset_type.value
    payload["metadata"]["format"] = dataset.metadata.format.value

    return payload


def _dataset_from_dict(payload: dict[str, Any]) -> BenchmarkDataset:
    payload = dict(payload)

    dataset_type = DatasetType(payload["dataset_type"])
    dataset_class = _DATASET_CLASSES[dataset_type.value]

    metadata_payload = dict(payload["metadata"])
    metadata_payload["format"] = DatasetFormat(metadata_payload["format"])

    payload["dataset_type"] = dataset_type
    payload["metadata"] = DatasetMetadata(**metadata_payload)

    return dataset_class(**payload)
