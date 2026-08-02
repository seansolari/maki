from __future__ import annotations

from collections.abc import Callable
from typing import Any

from .metadata import collect_metadata
from .metrics import MetricsCollector
from .schema import BenchmarkMetrics
from .schema import BenchmarkResult
from .stores.base import ResultStore


class BenchmarkRunner:
    def run(
        self,
        name: str,
        workload: Callable[[], Any],
        *,
        threads: int = 1,
        store: ResultStore | None = None,
    ) -> BenchmarkResult:

        metadata = collect_metadata(threads=threads)
        metrics = MetricsCollector()
        metrics.start()

        workload()

        (
            wall_time,
            cpu_time,
            peak_rss
        ) = metrics.stop()

        result = BenchmarkResult(
            benchmark_name=name,
            metadata=metadata,
            metrics=BenchmarkMetrics(
                wall_time_seconds=wall_time,
                cpu_time_seconds=cpu_time,
                peak_rss_bytes=peak_rss
            )
        )

        if store is not None:
            store.save(result)

        return result