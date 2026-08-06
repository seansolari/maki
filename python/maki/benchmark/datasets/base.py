from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .schema import DatasetMetadata, DatasetType


@dataclass(slots=True)
class BenchmarkDataset:
    dataset_id: str
    name: str
    version: str
    dataset_type: DatasetType
    metadata: DatasetMetadata = field(default_factory=DatasetMetadata)

    def cache_key(self) -> str:
        return f"{self.dataset_id}@{self.version}"


@dataclass(slots=True)
class RealDataset(BenchmarkDataset):
    urls: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.dataset_type = DatasetType.REAL


@dataclass(slots=True)
class LocalDataset(BenchmarkDataset):
    path: str = ""

    def __post_init__(self) -> None:
        self.dataset_type = DatasetType.LOCAL


@dataclass(slots=True)
class SyntheticDataset(BenchmarkDataset):
    generator: str = ""
    parameters: dict[str, Any] = field(default_factory=dict)
    seed: int | None = None

    def __post_init__(self) -> None:
        self.dataset_type = DatasetType.SYNTHETIC


@dataclass(slots=True)
class DerivedDataset(BenchmarkDataset):
    parent_dataset_id: str = ""
    transformations: list[dict[str, Any]] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.dataset_type = DatasetType.DERIVED
