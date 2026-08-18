from __future__ import annotations

from contextlib import AbstractContextManager
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
    
    if (None if not dataset else dataset.data_type) not in workflow.compatible_data_types:
        raise RuntimeError(
            f"Workflow [{getattr(workflow, "name", "")}] is incompatible with dataset type "
            f"{None if not dataset else dataset.data_type}, requires one "
            f"of [{", ".join(workflow.compatible_data_types)}]."
        )

    if dataset is not None:
        dataset_id = dataset.id
        
        mgr = materialise.DatasetMaterialiser.from_dataset(dataset)(**kwargs)
        current = mgr.materialise(dataset)

    metadata = collect_metadata(
        threads=getattr(workflow, "threads", 1),
    )
    
    time_cmd = TimeWrapper()
    phase_results: list[PhaseBenchmarkResult] = []
    
    for phase in workflow.phases():
        cmd, next_data = phase.transform(current)
        
        if isinstance(current, AbstractContextManager):
            with current:
                metrics = time_cmd.run(cmd)
        else:
            metrics = time_cmd.run(cmd)

        phase_results.append(
            PhaseBenchmarkResult(
                phase_name=phase.name,
                metrics=metrics,
            )
        )
        
        current = next_data

    result = WorkflowBenchmarkResult(
        workflow_name=getattr(workflow, "name", ""),
        dataset_id=dataset_id,
        phase_results=phase_results,
        total_metrics=BenchmarkMetrics.aggregate(r.metrics for r in phase_results),
        metadata=metadata,
        run_id=str(uuid.uuid4()),
    )

    if store is not None:
        store.save(result)

    return result
