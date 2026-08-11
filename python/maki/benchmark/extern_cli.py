
from pathlib import Path
from typing import List, Optional

import typer


app = typer.Typer(help="Benchmarking on real datasets")


@app.command("list")
def list_datasets(
    manifest: Path = typer.Option(..., "--manifest", "-m", help="Path to the datasets YAML manifest."),
    data_type: Optional[str] = typer.Option(None, "--data-type", help="Filter by data type, for example shotgun_metagenomics."),
    domain: Optional[str] = typer.Option(None, "--domain", help="Filter by broad domain, for example metagenomics or genomics."),
    organism: Optional[str] = typer.Option(None, "--organism", help="Filter by organism or community label."),
    assay: Optional[str] = typer.Option(None, "--assay", help="Filter by assay, for example illumina_paired_end."),
    tags: Optional[List[str]] = typer.Option(None, "--tag", help="Filter by tags in metadata.tags."),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Show additional dataset details."),
) -> None:
    """
    List available datasets, optionally filtered by metadata fields.
    """
    from maki.benchmark.data.catalog.manager import DatasetManager
    
    lib = DatasetManager(manifest, Path.cwd())

    matches = lib.list_datasets(
        domain=domain,
        data_type=data_type,
        organism=organism,
        assay=assay,
        tags=tags,
    )

    if not matches:
        typer.echo("No datasets matched the provided filters.")
        raise typer.Exit(code=0)

    for dataset in matches:
        typer.echo(f"{dataset.id}\t{dataset.name}\t{dataset.data_type}")

        if verbose:
            typer.echo(f"  domain: {dataset.domain}")
            typer.echo(f"  assay: {dataset.assay}")
            typer.echo(f"  organism: {dataset.organism}")
            typer.echo(f"  source: {dataset.source}")
            typer.echo(f"  license: {dataset.license}")

            if dataset.description:
                typer.echo(f"  description: {dataset.description}")

            tags = dataset.metadata.get("tags", [])
            if tags:
                typer.echo(f"  tags: {', '.join(tags)}")

            typer.echo(f"  files: {len(dataset.files)}")
            typer.echo("")


@app.command("run")
def run(
    dataset_id: str = typer.Argument(..., help="Dataset ID to run the workflow on."),
    manifest: Path = typer.Option(..., "--manifest", "-m", help="Path to the datasets YAML manifest."),
    role: list[str] = typer.Option([], "--role", "-r", help="Restrict selected files by role. Can be supplied multiple times."),
    cache_dir: Optional[Path] = typer.Option(None, "--cache-dir", help="Directory to use for cached dataset files."),
    force: bool = typer.Option(False, "--force", help="Force re-download or re-processing, depending on your implementation."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Show what would be run without executing the workflow."),
) -> None:
    """
    Select a dataset and prepare to run a workflow.

    This is intentionally minimal. Extend this function with your own
    download, cache, validation, and analysis logic.
    """
    from maki.benchmark.data.catalog.manager import DatasetManager
    
    lib = DatasetManager(manifest, cache_dir or Path.cwd())
    dataset = lib.get_dataset(dataset_id)

    selected_files = list(dataset.files)

    if role:
        requested_roles = set(role)
        selected_files = [
            file_record
            for file_record in selected_files
            if file_record.role in requested_roles
        ]

    if role and not selected_files:
        typer.echo(
            f"Dataset '{dataset_id}' has no files matching role(s): {', '.join(role)}",
            err=True,
        )
        raise typer.Exit(code=1)

    typer.echo(f"Selected dataset: {dataset.id}")
    typer.echo(f"Name: {dataset.name}")
    typer.echo(f"Data type: {dataset.data_type}")

    if cache_dir is not None:
        typer.echo(f"Cache directory: {cache_dir}")

    if selected_files:
        typer.echo("Selected files:")
        for file_record in selected_files:
            file_role = file_record.role
            filename = file_record.filename
            typer.echo(f"  - {file_role}: {filename}")
    else:
        typer.echo("No files listed for this dataset.")

    if dry_run:
        typer.echo("Dry run requested. Exiting before workflow execution.")
        raise typer.Exit(code=0)

    # Extend this section with your own workflow.
    #
    # Suggested future steps:
    # 1. Resolve cache location.
    # 2. Download selected files if not cached.
    # 3. Verify checksums.
    # 4. Pass local paths into your analysis pipeline.
    # 5. Write results to an output directory.

    typer.echo("Run command scaffold complete. Add workflow logic here.")


if __name__ == "__main__":
    app()
