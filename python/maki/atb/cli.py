import asyncio
import csv
from datetime import datetime
import gzip
from typing import Optional

import typer
from rich.console import Console
from rich.table import Table
from pathlib import Path
import maki.models.database as mdb

from .downloader import MakiAtbManifest
from .config import SQLITE_URL
from .importer import import_annotation_batches, import_assembly_batches
from .queries import QueryOptions, count_taxon, count_all, list_sample_rows, plan_download as plan_download_impl
from .sql import connect, ensure_db, inspect_schema
from ..db import build as build_impl, update as update_impl


console = Console()
app = typer.Typer(
    help="Helper scripts to interface with AllTheBacteria.\n\n"
         "Provides user-friendly interface to download genomes and "
         "annotation data."
)


@app.command()
def init(
    db_path: Path,
    remote_src: str = SQLITE_URL,
    force: bool = False
):
    """
    Download and prepare metadata database.
    """
    ensure_db(remote_src, db_path)
    typer.echo(f"SQLite database ready at {db_path}")
    
    db = connect(db_path)
    schema = inspect_schema(db.conn)
    
    import_assembly_batches(db, schema, force=force)
    import_annotation_batches(db, schema, force=force)
    
    typer.echo("Database initialisation complete")


@app.command()
def inspect(db_path: Path):
    """
    Inspect schema.
    """
    db = connect(db_path)
    schema = inspect_schema(db.conn)

    for table, cols in schema.items():
      typer.echo(f"Table: {table}")
      typer.echo(f"Columns: {", ".join(cols)}\n")


@app.command()
def count(
    db_path: Path,
    taxon: Optional[str] = None,
    high_quality: bool = False,
    has_assembly: bool = False,
    max_contamination: Optional[float] = None,
    min_completeness: Optional[float] = None,
    has_annotation: bool = False,
    limit: Optional[int] = None
):
    """
    Count assembly records matching search criteria.
    """
    db = connect(db_path)
    opts = QueryOptions(taxon=taxon, high_quality=high_quality, has_assembly=has_assembly, max_contamination=max_contamination, min_completeness=min_completeness, has_annotation=has_annotation)

    if taxon:
        typer.echo(count_taxon(db.conn, opts))
    else:
        table = Table("Species", "Count")
        
        rows = count_all(db.conn, opts)
        if limit:
          rows = rows[:limit]
        
        for sp, count in rows:
            table.add_row(sp, str(count))
      
        console.print(table)


@app.command()
def head(
    db_path: Path,
    taxon: Optional[str] = None,
    high_quality: bool = False,
    has_assembly: bool = False,
    max_contamination: Optional[float] = None,
    min_completeness: Optional[float] = None,
    has_annotation: bool = False,
    limit: int = 20
):
    """
    View the assembly table.
    """
    db = connect(db_path)
    opts = QueryOptions(taxon=taxon, high_quality=high_quality, has_assembly=has_assembly, max_contamination=max_contamination, min_completeness=min_completeness, has_annotation=has_annotation)
    cols, rows = list_sample_rows(db.conn, opts, limit)
    
    table = Table(*cols)
    for row in rows:
      table.add_row(*(str(v) for v in row))
    
    console.print(table)


@app.command()
def plan_download(
    db_path: Path,
    taxon: Optional[str] = None,
    high_quality: bool = True,
    has_assembly: bool = True,
    max_contamination: Optional[float] = 5.0,
    min_completeness: Optional[float] = 95.0,
    has_annotation: bool = True,
    outfile: Optional[str] = None
):
    """
    Create a download manifest.
    """
    db = connect(db_path)
    opts = QueryOptions(taxon=taxon, high_quality=high_quality, has_assembly=has_assembly, max_contamination=max_contamination, min_completeness=min_completeness, has_annotation=has_annotation)
    
    cols, rows = plan_download_impl(db.conn, opts)
    
    if not outfile:
        outfile = db_path.with_suffix(".manifest.csv.gz").name
    elif not outfile.endswith(".gz"):
        outfile += ".gz"
    
    typer.echo(f"Writing plan to {outfile}")
    
    with gzip.open(outfile, 'wt', newline='') as f:
        # Write metadata
        f.write(f"# {datetime.now()}\n")
        f.write(f"# database: {db_path}\n")
        if opts.taxon:
            f.write(f"# taxon: {opts.taxon}\n")
        f.write(f"# high quality filter: {opts.high_quality}\n")
        f.write(f"# assembly filter: {opts.has_assembly}\n")
        if opts.max_contamination:
            f.write(f"# max contamination: {opts.max_contamination}\n")
        if opts.min_completeness:
            f.write(f"# min completeness: {opts.min_completeness}\n")
        f.write(f"# annotation filter: {opts.has_annotation}\n")
        
        # Write data
        writer = csv.writer(f)
        writer.writerow(cols)
        writer.writerows(rows)


@app.command()
def download(
    manifest: Path,
    output_dir: Path = Path("downloads"),
    concurrency: int = 6
):
    """
    Download assembly and annotation data according to a download manifest created with `plan-download`.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    mani = MakiAtbManifest.from_csv(manifest, concurrency, output_dir)
    
    mani.download_batches(output_dir)
    
    mani.write_csv(output_dir / "manifest.csv.gz")
    typer.echo(f"Manifest exported to {output_dir / "manifest.csv.gz"}")
    
    
@app.command(help="""
Build a new database using a download manifest.

Modes:
- fixed: minimal size, cannot be updated
- updatable: stores compressed source data to allow updates
""")
def build(
    manifest_path: Path = typer.Option(..., help="Download package manifest"),
    db_path: Path = typer.Option(..., help="Database output path"),
    kmer_size: int = typer.Option(31),
    rank: str = typer.Option(..., help="Taxonomic rank"),
    threads: int = typer.Option(4),
    mode: str = typer.Option(
        "fixed",
        help="Database mode: fixed or updatable"
    ),
    force: bool = False
):
    manifest = MakiAtbManifest.from_csv(manifest_path, threads, db_path)
    
    return build_impl(manifest, db_path, kmer_size, rank, threads, mdb.UpdateMode[mode.lower()], mdb.TaxonomySource.gtdb, force)


@app.command(help="""
Update an existing database with new genomes from a download manifest.

- Performs incremental diff vs current manifest
- Only rebuilds affected clusters
- Preserves unchanged indices
""")
def update(
    manifest_path: Path = typer.Option(..., help="Download package manifest"),
    db_path: Path = typer.Option(..., help="Existing database"),
    threads: int = typer.Option(4, help="Parallel threads")
):
    manifest = MakiAtbManifest.from_csv(manifest_path, threads, db_path)
    
    return update_impl(manifest, db_path, threads)


def main():
    app()
  

if __name__ == "__main__":
    main()
