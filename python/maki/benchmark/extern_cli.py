import json
from pathlib import Path
from typing import Annotated

import typer

from maki.benchmark.workflow.cli import Workflow


app = typer.Typer(help="Benchmarking on real datasets")


DATASET_OPTS_PANEL = "Dataset Parameters"
ForwardFile = Annotated[
    Path,
    typer.Option("--forward",
                 help="Forward read-pair file.",
                 rich_help_panel=DATASET_OPTS_PANEL)
]
ReverseFile = Annotated[
    Path,
    typer.Option("--reverse",
                 help="Reverse read-pair file.",
                 rich_help_panel=DATASET_OPTS_PANEL)
]
Subsample = Annotated[
    int | None,
    typer.Option("--subsample",
                 help="Subsample input reads to this number.",
                 rich_help_panel=DATASET_OPTS_PANEL)
]
Seed = Annotated[
    int,
    typer.Option("--seed",
                 help="Seed for subsampling reads.",
                 rich_help_panel=DATASET_OPTS_PANEL)
]

IO_OPTS_PANEL = "I/O Parameters"
StatsFile = Annotated[
    Path,
    typer.Option("--stats-out",
                 help="File to write stats.",
                 rich_help_panel=IO_OPTS_PANEL)
]


@app.command(
    "run-pe",
    context_settings={
        "allow_extra_args": True,
        "ignore_unknown_options": True
    },
    help="Run benchmark workflow on paired-end data."
)
def run_paired_end(
    ctx: typer.Context,
    forward: ForwardFile,
    reverse: ReverseFile,
    workflow: Workflow,
    stats_out: StatsFile = Path("stats.json"),
    subsample: Subsample = None,
    seed: int = 1
):
    from maki.benchmark.data.materialise.local import PairedEndDatasetToken
    from maki.benchmark.runner import run_workflow
    from maki.benchmark.workflow.base import BenchmarkWorkflow
    from maki.benchmark.workflow.cli import parse_workflow_args
    
    extra_args = parse_workflow_args(ctx)
    workflow_mgr = BenchmarkWorkflow.from_name(workflow, extra_args)
    
    if "paired-end" not in workflow_mgr.compatible_data_types:
        typer.echo(f"Workflow {workflow_mgr.name} does not accept paired-end data.", err=True)
        return 1
    
    token = PairedEndDatasetToken.from_local("metagenomics", forward, reverse, subsample, seed)
    result = run_workflow(workflow_mgr, token)
    
    with stats_out.open("wt") as f:
        json.dump(result.to_dict(), f)
    

if __name__ == "__main__":
    app()
