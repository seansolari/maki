

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

from maki.benchmark.data.enums import SyntheticSequencingDatasetType

from .sequencing import *


@dataclass(frozen=True)
class SequencingSimulatorOptions:
    sequencing: SyntheticSequencingDatasetType
    read_length: int
    read_count: int
    sequencing_seed: int
    output_dir: Optional[Path]
    prefix: str
    insert_size: int
    min_read_length: int


def build_simulator(opts: SequencingSimulatorOptions) -> SequencingSimulator:
    """Construct the requested sequencing simulator from CLI parameters."""

    kwargs: dict[str, Any] = {
        "read_length": opts.read_length,
        "read_count": opts.read_count,
        "seed": opts.sequencing_seed,
        "output_dir": opts.output_dir,
        "prefix": opts.prefix,
    }

    if opts.sequencing == SyntheticSequencingDatasetType.PairedEnd:
        return ShortPairedEndSimulator(insert_size=opts.insert_size, **kwargs)
    elif opts.sequencing == SyntheticSequencingDatasetType.Unpaired:
        return ShortUnpairedSimulator(**kwargs)
    elif opts.sequencing == SyntheticSequencingDatasetType.LongRead:
        return LongReadSimulator(min_read_length=opts.min_read_length, **kwargs)
    else:
        raise RuntimeError(f"Sequencing must be one of: paired-end, unpaired, long-read. Received {opts.sequencing}.")
