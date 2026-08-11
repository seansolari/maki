"""Dataset object and convenience functions for synthetic benchmarks."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import shutil

from maki.benchmark.data.models import Dataset

from .topologies import SyntheticSequence


@dataclass
class SyntheticBenchmarkDataset(Dataset):
    """Generated benchmark dataset suitable for workflow registration layers.

    `files` maps logical roles such as `forward`, `reverse`, `unpaired`, or
    `long_reads` to generated FASTQ paths. When `cleanup_on_exit` is true and
    the simulator created a temporary directory, use this object as a context
    manager or call `cleanup()` when downstream workflows have finished.
    """

    sequences: list[SyntheticSequence] = field(default_factory=list)
    cleanup_on_exit: bool = True
    temporary_directory: Path | None = None

    def cleanup(self) -> None:
        """Remove generated temporary files if they are owned by this dataset."""

        if self.temporary_directory is not None and self.temporary_directory.exists():
            shutil.rmtree(self.temporary_directory)

    def __enter__(self) -> "SyntheticBenchmarkDataset":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        if self.cleanup_on_exit:
            self.cleanup()

