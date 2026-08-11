from __future__ import annotations

from typing import Optional
import uuid

from . import workflow
from .data import materialise
from .data.models import DatasetReference
from .stat import *


def run_workflow(
    workflow: "workflow.BenchmarkWorkflow",
    dataset: Optional["DatasetReference"] = None,
    *,
    threads: int = 1,
    store = None,
    **kwargs
) -> WorkflowBenchmarkResult:
    """
    Execute all workflow phases.

    Dataset is passed to the first phase.
    Output from each phase becomes input
    to the next phase.
    """
    
    current = dataset
    dataset_id = None
    
    if (None if not dataset else dataset.data_type) not in workflow.compatible_data_types():
        raise RuntimeError(
            f"Workflow [{getattr(workflow, "name", "")}] is incompatible with dataset type "
            f"{None if not dataset else dataset.data_type}, requires one "
            f"of [{", ".join(workflow.compatible_data_types())}]."
        )

    if dataset is not None:
        dataset_id = dataset.id
        
        mgr = materialise.DatasetMaterialiser.from_dataset(dataset)(**kwargs)
        materialised = mgr.materialise(dataset)
        current = materialised

    metadata = collect_metadata(
        threads=threads,
    )
    
    phase_results: list[PhaseBenchmarkResult] = []
    
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
        workflow_name=getattr(workflow, "name", ""),
        dataset_id=dataset_id,
        phase_results=phase_results,
        total_metrics=total_metrics,
        metadata=metadata,
        run_id=str(uuid.uuid4()),
    )

    if store is not None:
        store.save(result)

    return result
