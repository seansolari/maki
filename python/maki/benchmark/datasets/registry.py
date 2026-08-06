from __future__ import annotations

from .base import BenchmarkDataset
from .errors import DatasetAlreadyRegisteredError, DatasetNotFoundError


class DatasetRegistry:
    def __init__(self) -> None:
        self._datasets: dict[str, BenchmarkDataset] = {}

    def register(self, dataset: BenchmarkDataset) -> None:
        if dataset.dataset_id in self._datasets:
            raise DatasetAlreadyRegisteredError(
                f"Dataset already registered: {dataset.dataset_id}"
            )

        self._datasets[dataset.dataset_id] = dataset

    def get(self, dataset_id: str) -> BenchmarkDataset:
        try:
            return self._datasets[dataset_id]
        except KeyError as error:
            raise DatasetNotFoundError(
                f"Dataset not found: {dataset_id}"
            ) from error

    def list(self) -> list:
        return sorted(self._datasets.keys())

    def values(self) -> list:
        return [
            self._datasets[dataset_id]
            for dataset_id in self.list()
        ]

    def as_dict(self) -> dict[str, BenchmarkDataset]:
        return dict(self._datasets)
