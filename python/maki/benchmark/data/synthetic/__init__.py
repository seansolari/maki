"""Synthetic benchmark generation for graph benchmarking.

This submodule intentionally has no non-stdlib runtime dependencies. It creates
synthetic linear reference sequences with graph-stressing topology motifs, then
simulates lightweight FASTQ-style outputs for benchmark workflows.
"""

from .topologies import (
    SyntheticSequence,
    generate_linear_chains,
    generate_bubbles,
    generate_deep_branching,
    generate_high_degree_repeats,
    generate_random_sequences,
)


from .sequencing import (
    SequencingSimulator,
    SequencingResult,
    ShortPairedEndSimulator,
    ShortUnpairedSimulator,
    LongReadSimulator,
)


from .dataset import SyntheticBenchmarkDataset


__all__ = [
    "SyntheticSequence",
    "generate_linear_chains",
    "generate_bubbles",
    "generate_deep_branching",
    "generate_high_degree_repeats",
    "generate_random_sequences",
    "SequencingSimulator",
    "SequencingResult",
    "ShortPairedEndSimulator",
    "ShortUnpairedSimulator",
    "LongReadSimulator",
    "SyntheticBenchmarkDataset",
]
