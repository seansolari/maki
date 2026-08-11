"""Typer CLI for synthetic benchmark data generation.

Each topology is exposed as a top-level subcommand. Commands generate one
synthetic topology, pass it through a selected sequencing simulator, and print a
JSON payload containing generated file roles and metadata. The JSON payload is
intended to be consumed by benchmark workflow entry points.
"""

from __future__ import annotations

from pathlib import Path
from typing import Annotated
import typer

from maki.benchmark.data.enums import SyntheticSequencingDatasetType
from maki.benchmark.workflow.cli import Workflow


app = typer.Typer(
    help="Lightweight synthetic sequence data for graph benchmarks.",
    no_args_is_help=True,
)


SEQUENCING_OPTS_PANEL = "Sequencing Simulation Options"
SequencingKind = Annotated[
    SyntheticSequencingDatasetType,
    typer.Option("--sequencing", "-s",
                 help="Sequencing output type: paired-end, unpaired, or long-read.",
                 rich_help_panel=SEQUENCING_OPTS_PANEL)
]
ReadLength = Annotated[
    int,
    typer.Option("--read-length",
                 help="Read length for generated FASTQ records.",
                 rich_help_panel=SEQUENCING_OPTS_PANEL),
]
ReadCount = Annotated[
    int,
    typer.Option("--read-count",
                 help="Number of reads or read pairs to generate.",
                 rich_help_panel=SEQUENCING_OPTS_PANEL),
]
SequencingSeed = Annotated[
    int,
    typer.Option("--sequencing-seed",
                 help="Random seed for sequencing simulation.",
                 rich_help_panel=SEQUENCING_OPTS_PANEL),
]
InsertSize = Annotated[
    int,
    typer.Option("--insert-size",
                 help="Insert size for paired-end sequencing.",
                 rich_help_panel=SEQUENCING_OPTS_PANEL),
]
MinReadLength = Annotated[
    int,
    typer.Option("--min-readl-length",
                 help="Minimum read length for long-read sequencing.",
                 rich_help_panel=SEQUENCING_OPTS_PANEL),
]

GENERATION_OPTS_PANEL = "Sequence Generation Options"
SequenceCount = Annotated[
    int,
    typer.Option("--sequence-count",
                 help="Number of linear sequences to generate.",
                 rich_help_panel=GENERATION_OPTS_PANEL),
]
KmerSize = Annotated[
    int,
    typer.Option("-k", "--kmer-size",
                 help="Intended graph k-mer size.",
                 rich_help_panel=GENERATION_OPTS_PANEL),
]
SequenceLength = Annotated[
    int,
    typer.Option("--length",
                 help="Length of each generated linear sequence.",
                 rich_help_panel=GENERATION_OPTS_PANEL),
]
BranchCount = Annotated[
    int,
    typer.Option("--branch-count",
                 help="Number of nested branch levels.",
                 rich_help_panel=GENERATION_OPTS_PANEL)
]
BranchLength = Annotated[
    int,
    typer.Option("--branch-length",
                 help="Length of each divergent bubble branch.",
                 rich_help_panel=GENERATION_OPTS_PANEL)
]
RepeatCount = Annotated[
    int,
    typer.Option("--repeat-count",
                 help="Number of unique spacers joined by repeated cores.",
                 rich_help_panel=GENERATION_OPTS_PANEL)
]
TopologySeed = Annotated[
    int,
    typer.Option("--topology-seed",
                 help="Random seed for topology generation.",
                 rich_help_panel=GENERATION_OPTS_PANEL),
]

GENERAL_OPTONS = "General Options"
OutputDir = Annotated[
    Path | None,
    typer.Option("--output",
                 help="Directory for outputs. If omitted, a temporary directory is created.",
                 rich_help_panel=GENERAL_OPTONS),
]
OutputPrefix = Annotated[
    str,
    typer.Option("--prefix",
                 help="Prefix for generated FASTQ files.",
                 rich_help_panel=GENERAL_OPTONS),
]


@app.command("linear-chain")
def linear_chain(
    workflow: Workflow,
    sequence_count: SequenceCount = 4,
    length: SequenceLength = 2_000,
    topology_seed: TopologySeed = 1,
    sequencing: SequencingKind = SyntheticSequencingDatasetType.PairedEnd,
    read_length: ReadLength = 150,
    read_count: ReadCount = 10_000,
    sequencing_seed: SequencingSeed = 1,
    insert_size: InsertSize = 350,
    min_read_length: MinReadLength = 500,
    output_dir: OutputDir = None,
    prefix: OutputPrefix = "synthetic_linear_chain",
) -> None:
    from maki.benchmark.data.materialise.synthetic import SyntheticDatasetToken
    from maki.benchmark.data.synthetic.sequencing_registry import SequencingSimulatorOptions
    from maki.benchmark.data.synthetic import topologies
    from maki.benchmark.runner import run_workflow
    from maki.benchmark.workflow.base import BenchmarkWorkflow
    
    """Generate long linear chains that favour simple graph paths."""
    workflow_mgr = BenchmarkWorkflow.from_name(workflow)
    token = SyntheticDatasetToken(
        id="linear-chain",
        name="Linear Chain",
        description="Long linear chains that favour simple graph paths",
        domain="Synthetic Data",
        data_type=sequencing.value,
        sequences=topologies.generate_linear_chains(sequence_count, length, topology_seed),
        sequencing_opts=SequencingSimulatorOptions(
            sequencing,
            read_length,
            read_count,
            sequencing_seed,
            output_dir,
            prefix,
            insert_size,
            min_read_length
        )
    )
    result = run_workflow(workflow_mgr, token, cleanup_on_exit=True)
    print(result.to_json())


@app.command("bubble")
def bubble(
    workflow: Workflow,
    k: KmerSize = 31,
    sequence_count: SequenceCount = 4,
    branch_length: BranchLength = 120,
    topology_seed: TopologySeed = 1,
    sequencing: SequencingKind = SyntheticSequencingDatasetType.PairedEnd,
    read_length: ReadLength = 150,
    read_count: ReadCount = 10_000,
    sequencing_seed: SequencingSeed = 1,
    insert_size: InsertSize = 350,
    min_read_length: MinReadLength = 500,
    output_dir: OutputDir = None,
    prefix: OutputPrefix = "synthetic_linear_chain",
) -> None:
    from maki.benchmark.data.materialise.synthetic import SyntheticDatasetToken
    from maki.benchmark.data.synthetic.sequencing_registry import SequencingSimulatorOptions
    from maki.benchmark.data.synthetic import topologies
    from maki.benchmark.runner import run_workflow
    from maki.benchmark.workflow.base import BenchmarkWorkflow
    
    """Generate shared-flank divergent branches that form graph bubbles."""
    workflow_mgr = BenchmarkWorkflow.from_name(workflow)
    token = SyntheticDatasetToken(
        id="bubble",
        name="Bubble",
        description="Generate shared-flank divergent branches that form graph bubbles",
        domain="Synthetic Data",
        data_type=sequencing.value,
        sequences=topologies.generate_bubbles(k, sequence_count, branch_length, topology_seed),
        sequencing_opts=SequencingSimulatorOptions(
            sequencing,
            read_length,
            read_count,
            sequencing_seed,
            output_dir,
            prefix,
            insert_size,
            min_read_length
        )
    )
    result = run_workflow(workflow_mgr, token, cleanup_on_exit=True)
    print(result.to_json())


@app.command("deep-branching")
def deep_branching(
    workflow: Workflow,
    k: KmerSize = 31,
    sequence_count: SequenceCount = 4,
    branch_count: BranchCount = 4,
    branch_length: BranchLength = 120,
    topology_seed: TopologySeed = 1,
    sequencing: SequencingKind = SyntheticSequencingDatasetType.PairedEnd,
    read_length: ReadLength = 150,
    read_count: ReadCount = 10_000,
    sequencing_seed: SequencingSeed = 1,
    insert_size: InsertSize = 350,
    min_read_length: MinReadLength = 500,
    output_dir: OutputDir = None,
    prefix: OutputPrefix = "synthetic_linear_chain",
) -> None:
    from maki.benchmark.data.materialise.synthetic import SyntheticDatasetToken
    from maki.benchmark.data.synthetic.sequencing_registry import SequencingSimulatorOptions
    from maki.benchmark.data.synthetic import topologies
    from maki.benchmark.runner import run_workflow
    from maki.benchmark.workflow.base import BenchmarkWorkflow
        
    """Generate recursively shared anchors with nested branch segments."""
    workflow_mgr = BenchmarkWorkflow.from_name(workflow)
    token = SyntheticDatasetToken(
        id="deep-branching",
        name="Deep Branching",
        description="Generate recursively shared anchors with nested branch segments",
        domain="Synthetic Data",
        data_type=sequencing.value,
        sequences=topologies.generate_deep_branching(k, sequence_count, branch_count, branch_length, topology_seed),
        sequencing_opts=SequencingSimulatorOptions(
            sequencing,
            read_length,
            read_count,
            sequencing_seed,
            output_dir,
            prefix,
            insert_size,
            min_read_length
        )
    )
    result = run_workflow(workflow_mgr, token, cleanup_on_exit=True)
    print(result.to_json())


@app.command("high-degree-repeat")
def high_degree_repeat(
    workflow: Workflow,
    k: KmerSize = 31,
    sequence_count: SequenceCount = 4,
    branch_length: BranchLength = 120,
    repeat_count: RepeatCount = 8,
    topology_seed: TopologySeed = 1,
    sequencing: SequencingKind = SyntheticSequencingDatasetType.PairedEnd,
    read_length: ReadLength = 150,
    read_count: ReadCount = 10_000,
    sequencing_seed: SequencingSeed = 1,
    insert_size: InsertSize = 350,
    min_read_length: MinReadLength = 500,
    output_dir: OutputDir = None,
    prefix: OutputPrefix = "synthetic_linear_chain",
) -> None:
    from maki.benchmark.data.materialise.synthetic import SyntheticDatasetToken
    from maki.benchmark.data.synthetic.sequencing_registry import SequencingSimulatorOptions
    from maki.benchmark.data.synthetic import topologies
    from maki.benchmark.runner import run_workflow
    from maki.benchmark.workflow.base import BenchmarkWorkflow
    
    """Generate repeated cores with unique spacers to create high-degree merges."""
    workflow_mgr = BenchmarkWorkflow.from_name(workflow)
    token = SyntheticDatasetToken(
        id="high-degree-repeat",
        name="High Degree Repeat",
        description="Generate repeated cores with unique spacers to create high-degree merges",
        domain="Synthetic Data",
        data_type=sequencing.value,
        sequences=topologies.generate_high_degree_repeats(k, sequence_count, branch_length, repeat_count, topology_seed),
        sequencing_opts=SequencingSimulatorOptions(
            sequencing,
            read_length,
            read_count,
            sequencing_seed,
            output_dir,
            prefix,
            insert_size,
            min_read_length
        )
    )
    result = run_workflow(workflow_mgr, token, cleanup_on_exit=True)
    print(result.to_json())


@app.command("random")
def random_sequences(
    workflow: Workflow,
    sequence_count: SequenceCount = 4,
    length: SequenceLength = 2_000,
    topology_seed: TopologySeed = 1,
    sequencing: SequencingKind = SyntheticSequencingDatasetType.PairedEnd,
    read_length: ReadLength = 150,
    read_count: ReadCount = 10_000,
    sequencing_seed: SequencingSeed = 1,
    insert_size: InsertSize = 350,
    min_read_length: MinReadLength = 500,
    output_dir: OutputDir = None,
    prefix: OutputPrefix = "synthetic_linear_chain",
) -> None:
    from maki.benchmark.data.materialise.synthetic import SyntheticDatasetToken
    from maki.benchmark.data.synthetic.sequencing_registry import SequencingSimulatorOptions
    from maki.benchmark.data.synthetic import topologies
    from maki.benchmark.runner import run_workflow
    from maki.benchmark.workflow.base import BenchmarkWorkflow
    
    """Generate random control sequences with no deliberate graph motif."""
    workflow_mgr = BenchmarkWorkflow.from_name(workflow)
    token = SyntheticDatasetToken(
        id="random",
        name="Random",
        description="Generate random control sequences with no deliberate graph motif",
        domain="Synthetic Data",
        data_type=sequencing.value,
        sequences=topologies.generate_random_sequences(sequence_count, length, topology_seed),
        sequencing_opts=SequencingSimulatorOptions(
            sequencing,
            read_length,
            read_count,
            sequencing_seed,
            output_dir,
            prefix,
            insert_size,
            min_read_length
        )
    )
    result = run_workflow(workflow_mgr, token, cleanup_on_exit=True)
    print(result.to_json())
    

if __name__ == "__main__":
    app()
