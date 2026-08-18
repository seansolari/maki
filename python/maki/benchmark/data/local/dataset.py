
from dataclasses import dataclass
from pathlib import Path
import shutil

from maki.benchmark.data.models import Dataset


@dataclass
class LocalBenchmarkDataset(Dataset):
    temporary_directory: Path | None = None
    
    def cleanup(self) -> None:
        if self.temporary_directory is not None and self.temporary_directory.exists():
            shutil.rmtree(self.temporary_directory)

    def __enter__(self) -> "LocalBenchmarkDataset":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.cleanup()
