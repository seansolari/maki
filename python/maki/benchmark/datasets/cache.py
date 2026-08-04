from __future__ import annotations

from pathlib import Path


class DatasetCache:
    def __init__(self, root: Path):
        self.root = root
        self.root.mkdir(parents=True, exist_ok=True)

    def path_for(self, dataset_id: str) -> Path:
        return self.root / dataset_id

    def exists(self, dataset_id: str) -> bool:
        return self.path_for(dataset_id).exists()

    def remove(self, dataset_id: str) -> None:
        path = self.path_for(dataset_id)

        if path.exists():
            path.unlink()
