
from dataclasses import dataclass, field
from typing import Any


@dataclass
class Checksum:
    algorithm: str
    value: str


@dataclass
class DatasetFile:
    role: str
    url: str
    filename: str
    size_bytes: int | None = None
    checksum: Checksum | None = None


@dataclass
class Dataset:
    id: str
    name: str
    description: str | None
    domain: str
    data_type: str
    organism: str | None
    assay: str | None
    source: str | None
    license: str | None
    metadata: dict[str, Any] = field(default_factory=dict)
    files: list[DatasetFile] = field(default_factory=list)