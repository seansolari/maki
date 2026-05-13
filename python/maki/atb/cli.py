import sqlite3
from typing import Dict, List, Optional, Tuple

import typer
from rich.console import Console
from rich.table import Table
from pathlib import Path
from .config import DEFAULT_DB
from .db import connect, ensure_db
from .importer import import_csv_to_sqlite
from .schema import inspect_schema
from .columns import SampleColumns, resolve_columns
from .queries import count_taxon, count_all, list_sample_rows
from .download import plan_download as plan_download_impl


console = Console()
app = typer.Typer()


def load(db_path: Path) -> Tuple[sqlite3.Connection, Dict[str, List[str]], SampleColumns]:
    db_path = ensure_db(db_path)
    conn = connect(db_path)

    schema = inspect_schema(conn)
    mapping = resolve_columns(schema)

    return conn, schema, mapping


@app.command()
def init(db_path: Path = DEFAULT_DB):
    """
    Download and prepare metadata database.
    """
    ensure_db(db_path)
    typer.echo(f"Database ready at {db_path}")


@app.command()
def inspect(db_path: Path = DEFAULT_DB):
    """
    Inspect schema.
    """
    conn = connect(db_path)
    schema = inspect_schema(conn)

    for table, cols in schema.items():
      typer.echo(f"Table: {table}")
      typer.echo(f"Columns: {", ".join(cols)}\n")


@app.command()
def insert(
    file_path: Path,
    table_name: str,
    db_path: Path = DEFAULT_DB,
):
    """
    Import a CSV/TSV (.gz supported) into SQLite with deduplication.
    """
    conn = connect(db_path)
    table = import_csv_to_sqlite(conn, file_path, table_name)
    typer.echo(f"Imported into table: {table}")


@app.command()
def count(
    taxon: Optional[str] = typer.Option(None, help="Species name to count, otherwise count all unique species"),
    db_path: Path = DEFAULT_DB,
):
    conn, schema, mapping = load(db_path)

    if taxon:
        typer.echo(count_taxon(conn, mapping, taxon))
    else:
        typer.echo(count_all(conn, mapping))


@app.command()
def head(
    high_quality: bool = False,
    assemblies: bool = False,
    limit: int = 20,
    db_path: Path = DEFAULT_DB,
):
    conn, schema, mapping = load(db_path)
    rows = list_sample_rows(conn, mapping, high_quality, assemblies, limit)

    table = Table(*schema[mapping.table])
    for r in rows:
        table.add_row(*(str(v) for v in r))
        
    console.print(table)


@app.command()
def plan_download(
    taxon: str,
    high_quality: bool = True,
    assemblies: bool = True,
    db_path: Path = DEFAULT_DB,
):
    conn, schema, mapping = load(db_path)

    results = plan_download_impl(
        conn, schema, mapping, taxon, high_quality, assemblies
    )

    for archive, count in results:
        typer.echo(f"{archive} ({count} samples)")


def main():
  app()
  

if __name__ == "__main__":
    main()
