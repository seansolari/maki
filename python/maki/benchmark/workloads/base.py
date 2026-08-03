from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, Sequence


class BenchmarkPhase(ABC):
    """
    A single measurable stage within a benchmark workflow.

    The framework measures execution time and memory.
    Implementations only perform computation.
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """Human-readable phase identifier."""

    @abstractmethod
    def execute(
        self,
        data: Any,
    ) -> Any:
        """
        Execute the phase.

        Parameters
        ----------
        data
            Input object produced by the previous phase.

        Returns
        -------
        Any
            Output passed to the next phase.
        """
        

class BenchmarkWorkflow(ABC):
    """
    Ordered sequence of benchmark phases.

    The framework executes phases sequentially
    and automatically propagates outputs.
    """

    @property
    @abstractmethod
    def name(self) -> str:
        ...

    @abstractmethod
    def phases(
        self,
    ) -> Sequence[BenchmarkPhase]:
        """
        Return ordered workflow phases.
        """