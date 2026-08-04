from __future__ import annotations

from .base import BenchmarkDataset


class DatasetRegistry:
    def __init__(self):
        self._datasets: dict[str, BenchmarkDataset] = {}

    def register(self, dataset: BenchmarkDataset) -> None:
        if dataset.dataset_id in self._datasets:
            raise ValueError(
                f"Dataset already registered: "
                f"{dataset.dataset_id}"
            )

        self._datasets[dataset.dataset_id] = dataset

    def get(self, dataset_id: str) -> BenchmarkDataset:
        return self._datasets[dataset_id]

    def list(self) -> list:
      return sorted(self._datasets.keys())

    @property
    def datasets(self) -> dict[str, BenchmarkDataset]:
        return dict(self._datasets)
