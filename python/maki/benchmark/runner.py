from __future__ import annotations

import uuid

from .metadata import collect_metadata
from .metrics import MetricsCollector
from .schema import BenchmarkMetrics, PhaseBenchmarkResult, WorkflowBenchmarkResult


class BenchmarkRunner:

    def run_workflow(
        self,
        workflow,
        dataset=None,
        *,
        threads=1,
        store=None,
        materializer_registry=None,
    ) -> WorkflowBenchmarkResult:
        """
        Execute all workflow phases.

        Dataset is passed to the first phase.
        Output from each phase becomes input
        to the next phase.
        """

        metadata = collect_metadata(
            threads=threads,
        )
        
        current = dataset
        dataset_id = None
    
        if dataset is not None:
            dataset_id = getattr(dataset, "dataset_id", None)
    
            if materializer_registry is not None:
                materializer = materializer_registry.get(dataset)
                materialized = materializer.materialize(dataset)
                current = materialized

        phase_results: list[
            PhaseBenchmarkResult
        ] = []

        workflow_wall = 0.0
        workflow_cpu = 0.0
        workflow_peak_rss = 0

        for phase in workflow.phases():

            collector = MetricsCollector()
            
            collector.start()

            current = phase.execute(current)

            (
                wall_time,
                cpu_time,
                peak_rss,
            ) = collector.stop()

            metrics = BenchmarkMetrics(
                wall_time_seconds=wall_time,
                cpu_time_seconds=cpu_time,
                peak_rss_bytes=peak_rss,
            )

            phase_results.append(
                PhaseBenchmarkResult(
                    phase_name=phase.name,
                    metrics=metrics,
                )
            )

            workflow_wall += wall_time
            workflow_cpu += cpu_time

            workflow_peak_rss = max(
                workflow_peak_rss,
                peak_rss,
            )

        total_metrics = BenchmarkMetrics(
            wall_time_seconds=workflow_wall,
            cpu_time_seconds=workflow_cpu,
            peak_rss_bytes=workflow_peak_rss,
        )

        result = WorkflowBenchmarkResult(
            workflow_name=workflow.name,
            dataset_id=dataset_id,
            phase_results=phase_results,
            total_metrics=total_metrics,
            metadata=metadata,
            run_id=str(uuid.uuid4()),
        )

        if store is not None:
            store.save(result)

        return result
