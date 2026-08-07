"""Very small synthetic sequencing layer for benchmark data generation."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
import random
import tempfile
from typing import Sequence

from .topologies import SyntheticSequence

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1].upper()


@dataclass(frozen=True)
class SequencingResult:
    """Locations produced by a sequencing simulator."""

    files: dict[str, Path]
    data_type: str
    temporary_directory: Path | None = None


class SequencingSimulator(ABC):
    """Base class for minimal sequencing simulators.

    Subclasses expose `data_type` and return a dictionary mapping file roles to
    generated file locations. No attempt is made to simulate real instrument
    quality distributions; this layer is only intended to pressure-test graph
    construction and workflow plumbing.
    """

    @property
    @abstractmethod
    def data_type(self) -> str:
        ...

    def __init__(
        self,
        *,
        read_length: int,
        read_count: int,
        seed: int = 1,
        output_dir: str | Path | None = None,
        prefix: str = "synthetic",
    ) -> None:
        if read_length < 1:
            raise ValueError("read_length must be positive")
        if read_count < 1:
            raise ValueError("read_count must be positive")
        self.read_length = read_length
        self.read_count = read_count
        self.seed = seed
        self.output_dir = Path(output_dir) if output_dir is not None else None
        self.prefix = prefix

    @abstractmethod
    def generate(self, sequences: Sequence[SyntheticSequence]) -> SequencingResult:
        ...

    def _prepare_output_dir(self) -> tuple[Path, Path | None]:
        if self.output_dir is not None:
            self.output_dir.mkdir(parents=True, exist_ok=True)
            return self.output_dir, None
        tmp = Path(tempfile.mkdtemp(prefix=f"{self.prefix}_"))
        return tmp, tmp

    def _sample_read(self, rng: random.Random, sequences: Sequence[SyntheticSequence], length: int | None = None) -> tuple[str, str]:
        if not sequences:
            raise ValueError("at least one sequence is required")
        seq = rng.choice(sequences)
        n = min(length or self.read_length, len(seq.sequence))
        start = rng.randint(0, len(seq.sequence) - n)
        return seq.name, seq.sequence[start : start + n]

    @staticmethod
    def _write_fastq_record(handle, name: str, bases: str) -> None:
        quality = "I" * len(bases)
        handle.write(f"@{name}\n{bases}\n+\n{quality}\n")


class ShortUnpairedSimulator(SequencingSimulator):
    """Generate single-end short-read FASTQ data."""

    @property
    def data_type(self) -> str:
        return "unpaired_metagenome"

    def __init__(self, *, read_length: int = 150, read_count: int = 10_000, **kwargs) -> None:
        super().__init__(read_length=read_length, read_count=read_count, **kwargs)

    def generate(self, sequences: Sequence[SyntheticSequence]) -> SequencingResult:
        rng = random.Random(self.seed)
        out_dir, tmp = self._prepare_output_dir()
        path = out_dir / f"{self.prefix}.single.fastq"
        with path.open("w", encoding="utf-8") as handle:
            for i in range(self.read_count):
                source, read = self._sample_read(rng, sequences)
                self._write_fastq_record(handle, f"{source}_{i}/SE", read)
        return SequencingResult(files={"unpaired": path}, data_type=self.data_type, temporary_directory=tmp)


class ShortPairedEndSimulator(SequencingSimulator):
    """Generate paired-end short-read FASTQ data."""

    @property
    def data_type(self) -> str:
        return "shotgun_metagenome"

    def __init__(
        self,
        *,
        read_length: int = 150,
        read_count: int = 10_000,
        insert_size: int = 350,
        **kwargs,
    ) -> None:
        super().__init__(read_length=read_length, read_count=read_count, **kwargs)
        if insert_size < read_length:
            raise ValueError("insert_size must be at least read_length")
        self.insert_size = insert_size

    def generate(self, sequences: Sequence[SyntheticSequence]) -> SequencingResult:
        rng = random.Random(self.seed)
        out_dir, tmp = self._prepare_output_dir()
        fwd = out_dir / f"{self.prefix}.R1.fastq"
        rev = out_dir / f"{self.prefix}.R2.fastq"
        eligible = [s for s in sequences if len(s.sequence) >= self.insert_size]
        if not eligible:
            raise ValueError("no sequence is long enough for the requested insert_size")
        with fwd.open("w", encoding="utf-8") as h1, rev.open("w", encoding="utf-8") as h2:
            for i in range(self.read_count):
                seq = rng.choice(eligible)
                start = rng.randint(0, len(seq.sequence) - self.insert_size)
                fragment = seq.sequence[start : start + self.insert_size]
                r1 = fragment[: self.read_length]
                r2 = _revcomp(fragment[-self.read_length :])
                self._write_fastq_record(h1, f"{seq.name}_{i}/1", r1)
                self._write_fastq_record(h2, f"{seq.name}_{i}/2", r2)
        return SequencingResult(files={"forward": fwd, "reverse": rev}, data_type=self.data_type, temporary_directory=tmp)


class LongReadSimulator(SequencingSimulator):
    """Generate long-read-like unpaired FASTQ data."""

    @property
    def data_type(self) -> str:
        return "long_read_metagenome"

    def __init__(
        self,
        *,
        read_length: int = 5_000,
        read_count: int = 2_000,
        min_read_length: int = 500,
        **kwargs,
    ) -> None:
        super().__init__(read_length=read_length, read_count=read_count, **kwargs)
        if min_read_length < 1:
            raise ValueError("min_read_length must be positive")
        if min_read_length > read_length:
            raise ValueError("min_read_length cannot exceed read_length")
        self.min_read_length = min_read_length

    def generate(self, sequences: Sequence[SyntheticSequence]) -> SequencingResult:
        rng = random.Random(self.seed)
        out_dir, tmp = self._prepare_output_dir()
        path = out_dir / f"{self.prefix}.long.fastq"
        with path.open("w", encoding="utf-8") as handle:
            for i in range(self.read_count):
                length = rng.randint(self.min_read_length, self.read_length)
                source, read = self._sample_read(rng, sequences, length=length)
                self._write_fastq_record(handle, f"{source}_{i}/LR", read)
        return SequencingResult(files={"long_reads": path}, data_type=self.data_type, temporary_directory=tmp)
