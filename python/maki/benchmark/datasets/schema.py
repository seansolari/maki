from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any


class DatasetType(str, Enum):
    SYNTHETIC = "synthetic"
    ISOLATE_GENOME = "isolate_genome"
    SIMULATED_METAGENOME = "simulated_metagenome"
    REAL_METAGENOME = "real_metagenome"


@dataclass(slots=True)
class DatasetMetadata:
    description: str = ""
    source: str | None = None
    checksum: str | None = None
    extra: dict[str, Any] = field(default_factory=dict)
