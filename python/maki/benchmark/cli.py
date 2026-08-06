from __future__ import annotations

from pathlib import Path
from typing import Optional
import json

import typer

from .config import BenchmarkConfig
from .datasets import DatasetCatalog
from .operations import DatasetOperation
from .runner import BenchmarkRequest, BenchmarkRunner
from .synthetic import synthetic_generators
from .workflows import benchmark_workflows


app = typer.Typer(
    name="maki-benchmark",
    help="Run maki graph benchmarks on synthetic or real datasets.",
)


synthetic_app = typer.Typer(help="Inspect synthetic benchmark generators.")
dataset_app = typer.Typer(help="Inspect real dataset catalogs.")
workflow_app = typer.Typer(help="Inspect benchmark workflows.")

app.add_typer(synthetic_app, name="synthetic")
app.add_typer(dataset_app, name="datasets")
app.add_typer(workflow_app, name="workflows")


def _parse_json_object(value: Optional[str], *, name: str) -> dict:
    if value is None or value == "":
        return {}
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError as exc:
        raise typer.BadParameter(f"{name} must be valid JSON") from exc
    if not isinstance(parsed, dict):
        raise typer.BadParameter(f"{name} must be a JSON object")
    return parsed


def _parse_operations(values: list[str]) -> list:
    operations: list[DatasetOperation] = []

    for raw in values:
        try:
            item = json.loads(raw)
        except json.JSONDecodeError as exc:
            raise typer.BadParameter("operation must be valid JSON") from exc

        if not isinstance(item, dict) or "name" not in item:
            raise typer.BadParameter("operation must be an object with a name field")

        operations.append(
            DatasetOperation(
                name=item["name"],
                params=dict(item.get("params", {})),
            )
        )

    return operations


@synthetic_app.command("list")
def list_synthetic() -> None:
    """List registered synthetic data generators."""
    for name, generator in synthetic_generators.items():
        typer.echo(
            json.dumps(
                {
                    "name": name,
                    "modality": generator.modality,
                    "parameters": generator.describe_parameters(),
                },
                sort_keys=True,
            )
        )


@dataset_app.command("list")
def list_datasets(
    catalog: Path = typer.Option(..., "--catalog", help="Path to dataset catalog JSON."),
) -> None:
    """List real datasets in a file-backed catalog."""
    ds_catalog = DatasetCatalog(catalog)
    ds_catalog.load()

    for entry in ds_catalog.entries():
        typer.echo(
            json.dumps(
                {
                    "name": entry.name,
                    "version": entry.version,
                    "modality": entry.modality,
                    "tags": list(entry.tags),
                    "description": entry.description,
                },
                sort_keys=True,
            )
        )


@workflow_app.command("list")
def list_workflows() -> None:
    """List registered benchmark workflows."""
    for name, workflow in benchmark_workflows.items():
        typer.echo(
            json.dumps(
                {
                    "name": name,
                    "compatible_modalities": list(workflow.compatible_modalities),
                    "parameters": workflow.describe_parameters(),
                },
                sort_keys=True,
            )
        )


@app.command("run")
def run_benchmark(
    benchmark: str = typer.Option(..., "--benchmark", "-b"),
    store: Path = typer.Option(..., "--store", "-s", help="Directory where results are written."),
    dataset_kind: str = typer.Option(..., "--dataset-kind", help="synthetic or real"),
    dataset: str = typer.Option(..., "--dataset", "-d", help="Synthetic generator or real dataset name."),
    catalog: Optional[Path] = typer.Option(None, "--catalog", help="Required for real datasets."),
    benchmark_params: Optional[str] = typer.Option(None, "--benchmark-params", help="JSON object."),
    dataset_params: Optional[str] = typer.Option(None, "--dataset-params", help="JSON object."),
    operation: list[str] = typer.Option(
        [],
        "--operation",
        help='JSON object, for example {"name":"subsample_reads","params":{"n_reads":100000}}',
    ),
    threads: int = typer.Option(1, "--threads", "-t"),
    seed: int = typer.Option(1, "--seed"),
    cache_dir: Optional[Path] = typer.Option(None, "--cache-dir"),
) -> None:
    """
    Run a registered benchmark workflow.

    This command is intentionally non-interactive and emits JSON to stdout,
    making it usable from CI and future regression test harnesses.
    """

    if dataset_kind not in {"synthetic", "real"}:
        raise typer.BadParameter("--dataset-kind must be synthetic or real")

    config = BenchmarkConfig.default()
    if cache_dir is not None:
        config = BenchmarkConfig(
            cache_dir=cache_dir,
            results_dir=config.results_dir,
            catalog_path=catalog,
        )

    request = BenchmarkRequest(
        benchmark=benchmark,
        dataset_kind=dataset_kind,
        dataset_name=dataset,
        store=store,
        benchmark_params=_parse_json_object(benchmark_params, name="benchmark_params"),
        dataset_params=_parse_json_object(dataset_params, name="dataset_params"),
        operations=_parse_operations(operation),
        threads=threads,
        seed=seed,
        catalog_path=catalog,
    )

    result = BenchmarkRunner(config).run(request)
    typer.echo(json.dumps(result.to_json(), indent=2, sort_keys=True))


if __name__ == "__main__":
    app()