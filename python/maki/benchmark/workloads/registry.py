from __future__ import annotations
from typing import List

from .base import BenchmarkWorkflow


class WorkflowRegistry:
    """
    Registry of available benchmark workflows.
    """

    def __init__(self) -> None:
        self._workflows: dict[str, BenchmarkWorkflow] = {}

    def register(
        self,
        workflow: BenchmarkWorkflow,
    ) -> None:
        if workflow.name in self._workflows:
            raise ValueError(
                f"Workflow '{workflow.name}' already registered."
            )

        self._workflows[workflow.name] = workflow

    def get(
        self,
        name: str,
    ) -> BenchmarkWorkflow:
        try:
            return self._workflows[name]
        except KeyError as exc:
            raise KeyError(
                f"Unknown workflow '{name}'."
            ) from exc

    def list(self) -> List[str]:
        return sorted(self._workflows.keys())