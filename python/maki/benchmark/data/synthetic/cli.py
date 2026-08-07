"""Typer CLI for synthetic benchmark data generation.

Each topology is exposed as a top-level subcommand. Commands generate one
synthetic topology, pass it through a selected sequencing simulator, and print a
JSON payload containing generated file roles and metadata. The JSON payload is
intended to be consumed by benchmark workflow entry points.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

import typer

from .dataset import SyntheticBenchmarkDataset, generate_synthetic_benchmark
from .sequencing import LongReadSimulator, SequencingSimulator, ShortPairedEndSimulator, ShortUnpairedSimulator
from .topologies import generate_bubbles, generate_deep_branching, generate_high_degree_repeats, generate_linear_chains, generate_random_sequences, SyntheticSequence


app = typer.Typer(
    help="Generate lightweight synthetic sequence data for graph benchmarks.",
    no_args_is_help=True,
)


SequencingKind = typer.Option("paired-end", "--sequencing", "-s", help="Sequencing output type: paired-end, unpaired, or long-read.")


def _build_simulator(
    *,
    sequencing: str,
    read_length: int,
    read_count: int,
    sequencing_seed: int,
    output_dir: Optional[Path],
    prefix: str,
    insert_size: int,
    min_read_length: int,
) -> SequencingSimulator:
    """Construct the requested sequencing simulator from CLI parameters."""

    normalized = sequencing.strip().lower().replace("_", "-")
    kwargs: dict[str, Any] = {
        "read_length": read_length,
        "read_count": read_count,
        "seed": sequencing_seed,
        "output_dir": output_dir,
        "prefix": prefix,
    }

    if normalized in {"paired-end", "paired", "pe"}:
        return ShortPairedEndSimulator(insert_size=insert_size, **kwargs)
    if normalized in {"unpaired", "single-end", "single", "se"}:
        return ShortUnpairedSimulator(**kwargs)
    if normalized in {"long-read", "long", "longread", "lr"}:
        return LongReadSimulator(min_read_length=min_read_length, **kwargs)

    raise typer.BadParameter(
        "sequencing must be one of: paired-end, unpaired, long-read",
        param_hint="--sequencing",
    )


def _generate_for_topology(
    *,
    sequences: list[SyntheticSequence],
    sequencing: str,
    read_length: int,
    read_count: int,
    sequencing_seed: int,
    insert_size: int,
    min_read_length: int,
    output_dir: Optional[Path],
    prefix: str,
    cleanup_on_exit=True,
) -> SyntheticBenchmarkDataset:
    """Generate one topology and print role-to-file mappings."""

    simulator = _build_simulator(
        sequencing=sequencing,
        read_length=read_length,
        read_count=read_count,
        sequencing_seed=sequencing_seed,
        output_dir=output_dir,
        prefix=prefix,
        insert_size=insert_size,
        min_read_length=min_read_length,
    )
    dataset = generate_synthetic_benchmark(
        sequences=sequences,
        simulator=simulator,
        cleanup_on_exit=cleanup_on_exit,
    )
    return dataset


@app.command("linear-chain")
def linear_chain(
    sequence_count: int = typer.Option(4, help="Number of linear sequences to generate."),
    length: int = typer.Option(2_000, help="Length of each generated linear sequence."),
    topology_seed: int = typer.Option(1, help="Random seed for topology generation."),
    sequencing: str = SequencingKind,
    read_length: int = typer.Option(150, help="Read length for generated FASTQ records."),
    read_count: int = typer.Option(10_000, help="Number of reads or read pairs to generate."),
    sequencing_seed: int = typer.Option(1, help="Random seed for sequencing simulation."),
    insert_size: int = typer.Option(350, help="Insert size for paired-end sequencing."),
    min_read_length: int = typer.Option(500, help="Minimum read length for long-read sequencing."),
    output_dir: Optional[Path] = typer.Option(None, help="Directory for outputs. If omitted, a temporary directory is created."),
    prefix: str = typer.Option("synthetic_linear_chain", help="Prefix for generated FASTQ files."),
) -> None:
    """Generate long linear chains that favour simple graph paths."""
    sequences = generate_linear_chains(sequence_count, length, topology_seed)
    data = _generate_for_topology(
        sequences=sequences,
        sequencing=sequencing,
        read_length=read_length,
        read_count=read_count,
        sequencing_seed=sequencing_seed,
        insert_size=insert_size,
        min_read_length=min_read_length,
        output_dir=output_dir,
        prefix=prefix
    )


@app.command("bubble")
def bubble(
    k: int = typer.Option(31, help="Intended graph k-mer size."),
    sequence_count: int = typer.Option(4, help="Number of alternative bubble branches to generate."),
    branch_length: int = typer.Option(120, help="Length of each divergent bubble branch."),
    topology_seed: int = typer.Option(1, help="Random seed for topology generation."),
    sequencing: str = SequencingKind,
    read_length: int = typer.Option(150, help="Read length for generated FASTQ records."),
    read_count: int = typer.Option(10_000, help="Number of reads or read pairs to generate."),
    sequencing_seed: int = typer.Option(1, help="Random seed for sequencing simulation."),
    insert_size: int = typer.Option(350, help="Insert size for paired-end sequencing."),
    min_read_length: int = typer.Option(500, help="Minimum read length for long-read sequencing."),
    output_dir: Optional[Path] = typer.Option(None, help="Directory for outputs. If omitted, a temporary directory is created."),
    prefix: str = typer.Option("synthetic_bubble", help="Prefix for generated FASTQ files."),
) -> None:
    """Generate shared-flank divergent branches that form graph bubbles."""

    sequences = generate_bubbles(k, sequence_count, branch_length, topology_seed)
    data = _generate_for_topology(
        sequences=sequences,
        sequencing=sequencing,
        read_length=read_length,
        read_count=read_count,
        sequencing_seed=sequencing_seed,
        insert_size=insert_size,
        min_read_length=min_read_length,
        output_dir=output_dir,
        prefix=prefix
    )


@app.command("deep-branching")
def deep_branching(
    k: int = typer.Option(31, help="Intended graph k-mer size."),
    sequence_count: int = typer.Option(4, help="Number of deeply branching paths to generate."),
    branch_count: int = typer.Option(4, help="Number of nested branch levels."),
    branch_length: int = typer.Option(120, help="Base length of each nested branch segment."),
    topology_seed: int = typer.Option(1, help="Random seed for topology generation."),
    sequencing: str = SequencingKind,
    read_length: int = typer.Option(150, help="Read length for generated FASTQ records."),
    read_count: int = typer.Option(10_000, help="Number of reads or read pairs to generate."),
    sequencing_seed: int = typer.Option(1, help="Random seed for sequencing simulation."),
    insert_size: int = typer.Option(350, help="Insert size for paired-end sequencing."),
    min_read_length: int = typer.Option(500, help="Minimum read length for long-read sequencing."),
    output_dir: Optional[Path] = typer.Option(None, help="Directory for outputs. If omitted, a temporary directory is created."),
    prefix: str = typer.Option("synthetic_deep_branching", help="Prefix for generated FASTQ files."),
) -> None:
    """Generate recursively shared anchors with nested branch segments."""

    sequences = generate_deep_branching(k, sequence_count, branch_count, branch_length, topology_seed)
    data = _generate_for_topology(
        sequences=sequences,
        sequencing=sequencing,
        read_length=read_length,
        read_count=read_count,
        sequencing_seed=sequencing_seed,
        insert_size=insert_size,
        min_read_length=min_read_length,
        output_dir=output_dir,
        prefix=prefix
    )


@app.command("high-degree-repeat")
def high_degree_repeat(
    k: int = typer.Option(31, help="Intended graph k-mer size."),
    sequence_count: int = typer.Option(4, help="Number of repeat-rich sequences to generate."),
    branch_length: int = typer.Option(120, help="Length of unique spacers between repeated cores."),
    repeat_count: int = typer.Option(8, help="Number of unique spacers joined by repeated cores."),
    topology_seed: int = typer.Option(1, help="Random seed for topology generation."),
    sequencing: str = SequencingKind,
    read_length: int = typer.Option(150, help="Read length for generated FASTQ records."),
    read_count: int = typer.Option(10_000, help="Number of reads or read pairs to generate."),
    sequencing_seed: int = typer.Option(1, help="Random seed for sequencing simulation."),
    insert_size: int = typer.Option(350, help="Insert size for paired-end sequencing."),
    min_read_length: int = typer.Option(500, help="Minimum read length for long-read sequencing."),
    output_dir: Optional[Path] = typer.Option(None, help="Directory for outputs. If omitted, a temporary directory is created."),
    prefix: str = typer.Option("synthetic_high_degree_repeat", help="Prefix for generated FASTQ files."),
) -> None:
    """Generate repeated cores with unique spacers to create high-degree merges."""

    sequences = generate_high_degree_repeats(k, sequence_count, branch_length, repeat_count, topology_seed)
    data = _generate_for_topology(
        sequences=sequences,
        sequencing=sequencing,
        read_length=read_length,
        read_count=read_count,
        sequencing_seed=sequencing_seed,
        insert_size=insert_size,
        min_read_length=min_read_length,
        output_dir=output_dir,
        prefix=prefix
    )


@app.command("random")
def random_sequences(
    sequence_count: int = typer.Option(4, help="Number of random control sequences to generate."),
    length: int = typer.Option(2_000, help="Length of each random sequence."),
    topology_seed: int = typer.Option(1, help="Random seed for topology generation."),
    sequencing: str = SequencingKind,
    read_length: int = typer.Option(150, help="Read length for generated FASTQ records."),
    read_count: int = typer.Option(10_000, help="Number of reads or read pairs to generate."),
    sequencing_seed: int = typer.Option(1, help="Random seed for sequencing simulation."),
    insert_size: int = typer.Option(350, help="Insert size for paired-end sequencing."),
    min_read_length: int = typer.Option(500, help="Minimum read length for long-read sequencing."),
    output_dir: Optional[Path] = typer.Option(None, help="Directory for outputs. If omitted, a temporary directory is created."),
    prefix: str = typer.Option("synthetic_random", help="Prefix for generated FASTQ files."),
) -> None:
    """Generate random control sequences with no deliberate graph motif."""

    sequences = generate_random_sequences(sequence_count, length, topology_seed)
    data = _generate_for_topology(
        sequences=sequences,
        sequencing=sequencing,
        read_length=read_length,
        read_count=read_count,
        sequencing_seed=sequencing_seed,
        insert_size=insert_size,
        min_read_length=min_read_length,
        output_dir=output_dir,
        prefix=prefix
    )


if __name__ == "__main__":
    app()
