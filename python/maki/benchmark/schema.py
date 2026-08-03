from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(slots=True)
class DatasetReference:
    id: str


@dataclass(slots=True)
class Dataset:
    id: DatasetReference


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
class PhaseBenchmarkResult:
    """
    Metrics for an individual workflow stage.
    """
    phase_name: str
    metrics: BenchmarkMetrics

    def to_dict(self) -> dict:
        return {
            "phase_name": self.phase_name,
            "metrics": asdict(self.metrics),
        }
        
    def to_json(self) -> dict:
        return self.to_dict()
    
    
@dataclass(slots=True)
class WorkflowBenchmarkResult:
    """
    Benchmark results for an entire workflow.
    """
    workflow_name: str
    dataset_id: DatasetReference
    phase_results: list[PhaseBenchmarkResult]
    total_metrics: BenchmarkMetrics
    metadata: BenchmarkMetadata
    schema_version: int = 2
    run_id: str = ""

    def to_dict(self) -> dict:
        return {
            "workflow_name": self.workflow_name,
            "dataset_id": asdict(self.dataset_id),
            "phase_results": [
                p.to_dict()
                for p in self.phase_results
            ],
            "total_metrics": asdict(
                self.total_metrics
            ),
            "metadata": asdict(
                self.metadata
            ),
            "schema_version": self.schema_version,
            "run_id": self.run_id,
        }
        
    def to_json(self) -> dict:
        return self.to_dict()
