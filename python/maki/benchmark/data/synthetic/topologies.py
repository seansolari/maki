"""Synthetic linear sequence generation for graph topology benchmarks."""

from __future__ import annotations

from dataclasses import dataclass
import random


__all__ = [
    "generate_bubbles",
    "generate_deep_branching",
    "generate_high_degree_repeats",
    "generate_linear_chains",
    "generate_random_sequences"
]

_DNA = "ACGT"


@dataclass(frozen=True)
class SyntheticSequence:
    """A generated linear sequence and metadata describing its purpose."""

    name: str
    sequence: str
    topology: str

    @property
    def metadata(self) -> dict[str, str | int]:
        return {"name": self.name, "topology": self.topology, "length": len(self.sequence)}


def _rng(seed: int) -> random.Random:
    return random.Random(seed)


def _random_dna(length: int, rng: random.Random) -> str:
    return "".join(rng.choice(_DNA) for _ in range(length))


def _nonempty_chunks(sequence: str, size: int) -> list[str]:
    return [sequence[i : i + size] for i in range(0, len(sequence), size) if sequence[i : i + size]]


def generate_linear_chains(sequence_count: int, length: int, seed: int) -> list[SyntheticSequence]:
    """Generate long low-repeat chains that favour simple paths in a de Bruijn graph."""

    rng = _rng(seed)
    return [
        SyntheticSequence(
            name=f"linear_chain_{i}",
            sequence=_random_dna(length, rng),
            topology="linear_chain",
        )
        for i in range(sequence_count)
    ]


def generate_bubbles(k: int, sequence_count: int, branch_length: int, seed: int) -> list[SyntheticSequence]:
    """Generate references with shared flanks and divergent interiors.

    Each output is still a linear sequence. Across the set, shared left/right
    flanks plus different middle branches create bubble-like alternatives after
    graph construction.
    """

    rng = _rng(seed + 101)
    flank = _random_dna(k * 2, rng)
    right = _random_dna(k * 2, rng)
    seqs: list[SyntheticSequence] = []
    for i in range(sequence_count):
        branch = _random_dna(branch_length, rng)
        seqs.append(SyntheticSequence(f"bubble_branch_{i}", flank + branch + right, "bubble"))
    return seqs


def generate_deep_branching(k: int, sequence_count: int, branch_count: int, branch_length: int, seed: int) -> list[SyntheticSequence]:
    """Generate nested branching paths using recursively shared anchors."""

    rng = _rng(seed + 202)
    anchors = [_random_dna(k, rng) for _ in range(branch_count + 1)]
    seqs: list[SyntheticSequence] = []
    for i in range(sequence_count):
        pieces: list[str] = []
        for depth in range(branch_count):
            pieces.append(anchors[depth])
            pieces.append(_random_dna(branch_length + i + depth, rng))
        pieces.append(anchors[-1])
        seqs.append(SyntheticSequence(f"deep_branch_{i}", "".join(pieces), "deep_branching"))
    return seqs


def generate_high_degree_repeats(k: int, sequence_count: int, branch_length: int, repeat_count: int, seed: int) -> list[SyntheticSequence]:
    """Generate repeated cores with unique spacers to create high-degree merges."""

    rng = _rng(seed + 303)
    repeat = _random_dna(k * 2, rng)
    seqs: list[SyntheticSequence] = []
    for i in range(sequence_count):
        spacers = [_random_dna(branch_length + j, rng) for j in range(repeat_count)]
        sequence = repeat.join(spacers)
        seqs.append(SyntheticSequence(f"high_degree_repeat_{i}", sequence, "high_degree_repeat"))
    return seqs


def generate_random_sequences(sequence_count: int, length: int, seed: int) -> list[SyntheticSequence]:
    """Generate random controls with no deliberate graph motif."""

    rng = _rng(seed + 404)
    return [
        SyntheticSequence(f"random_{i}", _random_dna(length, rng), "random")
        for i in range(sequence_count)
    ]
