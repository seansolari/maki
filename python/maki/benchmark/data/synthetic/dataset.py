"""Dataset object and convenience functions for synthetic benchmarks."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import shutil

from .sequencing import SequencingResult, SequencingSimulator
from .topologies import SyntheticSequence


@dataclass
class SyntheticBenchmarkDataset:
    """Generated benchmark dataset suitable for workflow registration layers.

    `files` maps logical roles such as `forward`, `reverse`, `unpaired`, or
    `long_reads` to generated FASTQ paths. When `cleanup_on_exit` is true and
    the simulator created a temporary directory, use this object as a context
    manager or call `cleanup()` when downstream workflows have finished.
    """

    sequences: list[SyntheticSequence]
    files: dict[str, Path]
    data_type: str
    cleanup_on_exit: bool = True
    temporary_directory: Path | None = None
    metadata: dict = field(default_factory=dict)

    def cleanup(self) -> None:
        """Remove generated temporary files if they are owned by this dataset."""

        if self.temporary_directory is not None and self.temporary_directory.exists():
            shutil.rmtree(self.temporary_directory)

    def __enter__(self) -> "SyntheticBenchmarkDataset":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        if self.cleanup_on_exit:
            self.cleanup()


def generate_synthetic_benchmark(
    *,
    sequences: list[SyntheticSequence],
    simulator: SequencingSimulator,
    cleanup_on_exit: bool = True,
) -> SyntheticBenchmarkDataset:
    """Generate canonical synthetic topology sequences and simulated reads."""

    result: SequencingResult = simulator.generate(sequences)
    return SyntheticBenchmarkDataset(
        sequences=sequences,
        files=result.files,
        data_type=result.data_type,
        cleanup_on_exit=cleanup_on_exit,
        temporary_directory=result.temporary_directory,
        metadata={
            "topologies": sorted({s.topology for s in sequences}),
            "sequence_count": len(sequences),
        },
    )
