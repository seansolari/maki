from __future__ import annotations

from dataclasses import dataclass, field

from .schema import DatasetMetadata, DatasetType


@dataclass(slots=True)
class BenchmarkDataset:
    dataset_id: str
    name: str
    version: str
    dataset_type: DatasetType

    tags: list[str] = field(default_factory=list)
    metadata: DatasetMetadata = field(
        default_factory=DatasetMetadata
    )


@dataclass(slots=True)
class DerivedDataset(BenchmarkDataset):
    parent_dataset_id: str = ""
    transformations: list[str] = field(
        default_factory=list
    )
