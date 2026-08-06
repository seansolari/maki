from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any


CATALOG_SCHEMA_VERSION = 1
CACHE_MANIFEST_SCHEMA_VERSION = 1


class DatasetType(str, Enum):
    SYNTHETIC = "synthetic"
    REAL = "real"
    LOCAL = "local"
    DERIVED = "derived"


class DatasetFormat(str, Enum):
    FASTQ = "fastq"
    FASTQ_GZ = "fastq.gz"
    FASTA = "fasta"
    FASTA_GZ = "fasta.gz"
    GRAPH = "graph"
    DIRECTORY = "directory"
    OTHER = "other"


@dataclass(slots=True)
class DatasetMetadata:
    description: str = ""
    source: str | None = None
    checksum: str | None = None
    checksum_algorithm: str = "sha256"
    format: DatasetFormat = DatasetFormat.OTHER
    tags: list[str] = field(default_factory=list)
    extra: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class MaterializedDataset:
    dataset_id: str
    version: str
    path: str
    metadata: DatasetMetadata
    