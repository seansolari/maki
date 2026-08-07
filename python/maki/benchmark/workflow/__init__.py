
from collections.abc import Mapping
from typing import List

from .base import BenchmarkPhase, BenchmarkWorkflow
from .graph.build import GraphBuildWorkflow
from .graph.unitig import UnitigWorkflow


workflow_instances: List[BenchmarkWorkflow] = [
    GraphBuildWorkflow(),
    UnitigWorkflow()
]


all_workflows: Mapping[str, BenchmarkWorkflow] = {
    workflow.name: workflow
    for workflow in workflow_instances
}


__all__ = [
    "BenchmarkPhase",
    "BenchmarkWorkflow",
    "all_workflows"
]
