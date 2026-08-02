from __future__ import annotations

from dataclasses import asdict
from dataclasses import dataclass
from dataclasses import field
from datetime import datetime
from uuid import uuid4


@dataclass(slots=True)
class BenchmarkMetrics:
    wall_time_seconds: float
    cpu_time_seconds: float
    peak_rss_bytes: int


@dataclass(slots=True)
class BenchmarkMetadata:
    git_commit: str | None
    git_branch: str | None
    python_version: str
    hostname: str
    cpu_count: int
    thread_count: int
    timestamp_utc: str


@dataclass(slots=True)
class BenchmarkResult:
    benchmark_name: str
    metrics: BenchmarkMetrics
    metadata: BenchmarkMetadata
    schema_version: int = 1
    run_id: str = field(default_factory=lambda: str(uuid4()))

    def to_dict(self) -> dict:
        return asdict(self)

    def to_json(self) -> dict:
        return self.to_dict()

