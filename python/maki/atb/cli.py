import csv
import gzip
from typing import Optional

import typer
from rich.console import Console
from rich.table import Table
from pathlib import Path
from .config import DEFAULT_DB
from .db import connect, ensure_db
from .importer import import_csv_to_sqlite
from .schema import inspect_schema
from .queries import QueryOptions, count_taxon, count_all, list_sample_rows, plan_download as plan_download_impl


console = Console()
app = typer.Typer()


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
    taxon: Optional[str] = None,
    high_quality: bool = False,
    has_assembly: bool = False,
    max_contamination: Optional[float] = None,
    min_completeness: Optional[float] = None,
    has_annotation: bool = False,
    limit: Optional[int] = None,
    db_path: Path = DEFAULT_DB,
):
    conn = connect(db_path)
    opts = QueryOptions(taxon=taxon, high_quality=high_quality, has_assembly=has_assembly, max_contamination=max_contamination, min_completeness=min_completeness, has_annotation=has_annotation)

    if taxon:
        typer.echo(count_taxon(conn, opts))
    else:
        table = Table("Species", "Count")
        
        rows = count_all(conn, opts)
        if limit:
          rows = rows[:limit]
        
        for sp, count in rows:
            table.add_row(sp, str(count))
      
        console.print(table)


@app.command()
def head(
    taxon: Optional[str] = None,
    high_quality: bool = False,
    has_assembly: bool = False,
    max_contamination: Optional[float] = None,
    min_completeness: Optional[float] = None,
    has_annotation: bool = False,
    limit: int = 20,
    db_path: Path = DEFAULT_DB,
):
    conn = connect(db_path)
    opts = QueryOptions(taxon=taxon, high_quality=high_quality, has_assembly=has_assembly, max_contamination=max_contamination, min_completeness=min_completeness, has_annotation=has_annotation)
    cols, rows = list_sample_rows(conn, opts, limit)
    
    table = Table(*cols)
    for row in rows:
      table.add_row(*(str(v) for v in row))
    
    console.print(table)


@app.command()
def plan_download(
    taxon: Optional[str] = None,
    high_quality: bool = True,
    has_assembly: bool = True,
    max_contamination: Optional[float] = 5.0,
    min_completeness: Optional[float] = 95.0,
    has_annotation: bool = True,
    outfile: Optional[str] = None,
    db_path: Path = DEFAULT_DB,
):
    conn = connect(db_path)
    opts = QueryOptions(taxon=taxon, high_quality=high_quality, has_assembly=has_assembly, max_contamination=max_contamination, min_completeness=min_completeness, has_annotation=has_annotation)
    
    cols, rows = plan_download_impl(conn, opts)
    
    if not outfile:
        outfile = db_path.with_suffix(".manifest.csv.gz").name
    elif not outfile.endswith(".gz"):
        outfile += ".gz"
    
    typer.echo(f"Writing plan to {outfile}")
    
    with gzip.open(outfile, 'wt', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(cols)
        writer.writerows(rows)


def main():
  app()
  

if __name__ == "__main__":
    main()
