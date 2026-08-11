

from enum import Enum


class SyntheticSequencingDatasetType(str, Enum):
    PairedEnd = "paired-end"
    Unpaired = "unpaired"
    LongRead = "long-read"
