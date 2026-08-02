from __future__ import annotations

from abc import ABC
from abc import abstractmethod

from ..schema import BenchmarkResult


class ResultStore(ABC):
    @abstractmethod
    def save(self, result: BenchmarkResult) -> None:
        ...
        