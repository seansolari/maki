from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Iterable


@dataclass(slots=True)
class BenchmarkMetrics:
    """
    Statistics returned by GNU /usr/bin/time.
    """

    elapsed_seconds: float
    user_seconds: float
    system_seconds: float
    cpu_percent: float

    max_rss_kb: int

    major_page_faults: int
    minor_page_faults: int

    voluntary_context_switches: int
    involuntary_context_switches: int

    filesystem_inputs: int
    filesystem_outputs: int

    exit_code: int

    samples: int = 1

    def __add__(self, other: "BenchmarkMetrics") -> "BenchmarkMetrics":
        if not isinstance(other, BenchmarkMetrics):
            return NotImplemented

        total_samples = self.samples + other.samples

        return BenchmarkMetrics(
            elapsed_seconds=self.elapsed_seconds + other.elapsed_seconds,
            user_seconds=self.user_seconds + other.user_seconds,
            system_seconds=self.system_seconds + other.system_seconds,
            cpu_percent=(
                (self.cpu_percent * self.samples)
                + (other.cpu_percent * other.samples)
            )
            / total_samples,
            max_rss_kb=max(self.max_rss_kb, other.max_rss_kb),
            major_page_faults=self.major_page_faults + other.major_page_faults,
            minor_page_faults=self.minor_page_faults + other.minor_page_faults,
            voluntary_context_switches=(
                self.voluntary_context_switches
                + other.voluntary_context_switches
            ),
            involuntary_context_switches=(
                self.involuntary_context_switches
                + other.involuntary_context_switches
            ),
            filesystem_inputs=(
                self.filesystem_inputs + other.filesystem_inputs
            ),
            filesystem_outputs=(
                self.filesystem_outputs + other.filesystem_outputs
            ),
            exit_code=other.exit_code,
            samples=total_samples,
        )

    def __iadd__(self, other: "BenchmarkMetrics") -> "BenchmarkMetrics":
        combined = self + other
        self.__dict__.update(combined.__dict__)
        return self

    @classmethod
    def aggregate(cls, stats: Iterable["BenchmarkMetrics"]) -> "BenchmarkMetrics":
        iterator = iter(stats)

        try:
            result = next(iterator)
        except StopIteration:
            raise ValueError("No statistics supplied")

        for item in iterator:
            result = result + item

        return result

    @property
    def average_elapsed_seconds(self) -> float:
        return self.elapsed_seconds / self.samples

    @property
    def average_user_seconds(self) -> float:
        return self.user_seconds / self.samples

    @property
    def average_system_seconds(self) -> float:
        return self.system_seconds / self.samples


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
    dataset_id: str | None
    phase_results: list[PhaseBenchmarkResult]
    total_metrics: BenchmarkMetrics
    metadata: BenchmarkMetadata
    schema_version: int = 2
    run_id: str = ""

    def to_dict(self) -> dict:
        return {
            "workflow_name": self.workflow_name,
            "dataset_id": self.dataset_id,
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
