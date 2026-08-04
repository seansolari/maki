from __future__ import annotations

from abc import ABC, abstractmethod

from .base import BenchmarkDataset


class DatasetMaterializer(ABC):

    @abstractmethod
    def materialize(
        self,
        dataset: BenchmarkDataset,
    ) -> BenchmarkDataset:
        ...
