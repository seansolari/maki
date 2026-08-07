
from pathlib import Path
import os
import shutil
import time


class DatasetCache:
    def __init__(self, cache_dir: str | Path, max_size_bytes: int):
        self.cache_dir = Path(cache_dir).expanduser().resolve()
        self.max_size_bytes = max_size_bytes
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def dataset_dir(self, dataset_id: str) -> Path:
        path = self.cache_dir / dataset_id
        path.mkdir(parents=True, exist_ok=True)
        return path

    def path_for(self, dataset_id: str, filename: str) -> Path:
        return self.dataset_dir(dataset_id) / filename

    def touch(self, path: str | Path) -> None:
        now = time.time()
        os.utime(path, (now, now))

    def total_size(self) -> int:
        total = 0
        for path in self.cache_dir.rglob("*"):
            if path.is_file():
                total += path.stat().st_size
        return total

    def prune(self) -> None:
        total = self.total_size()

        if total <= self.max_size_bytes:
            return

        files = [
            path
            for path in self.cache_dir.rglob("*")
            if path.is_file()
        ]

        files.sort(key=lambda p: p.stat().st_mtime)

        for path in files:
            if total <= self.max_size_bytes:
                break

            size = path.stat().st_size
            path.unlink()
            total -= size

        self._remove_empty_dirs()

    def clear(self) -> None:
        if self.cache_dir.exists():
            shutil.rmtree(self.cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def _remove_empty_dirs(self) -> None:
        dirs = sorted(
            [p for p in self.cache_dir.rglob("*") if p.is_dir()],
            key=lambda p: len(p.parts),
            reverse=True,
        )

        for d in dirs:
            try:
                d.rmdir()
            except OSError:
                pass