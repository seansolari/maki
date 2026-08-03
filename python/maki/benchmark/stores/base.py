from __future__ import annotations

from abc import ABC
from abc import abstractmethod

from ..schema import WorkflowBenchmarkResult


class ResultStore(ABC):

    @abstractmethod
    def save(
        self,
        result: WorkflowBenchmarkResult,
    ) -> None:
        ...
        