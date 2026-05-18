import asyncio
import csv
from datetime import datetime
import gzip
from typing import Optional

import typer
from rich.console import Console
from rich.table import Table
from pathlib import Path
from .config import DEFAULT_DB, SQLITE_URL
from .db import connect, ensure_db, inspect_schema
from .importer import import_annotation_batches, import_assembly_batches
from .queries import QueryOptions, count_taxon, count_all, list_sample_rows, plan_download as plan_download_impl
from .async_downloader import AssemblyItem, BaktaItem, Manifest, run_jobs_async
from .lists import AssemblyFileLists, BaktaFileLists


console = Console()
app = typer.Typer(
    help="Helper scripts to interface with AllTheBacteria.\n\n"
         "Provides user-friendly interface to download genomes and "
         "annotation data."
)


@app.command()
def init(
    db_path: Path = DEFAULT_DB,
    remote_src: str = SQLITE_URL,
    force: bool = False
):
    """
    Download and prepare metadata database.
    """
    ensure_db(remote_src, db_path)
    typer.echo(f"SQLite database ready at {db_path}")
    
    conn = connect(db_path)
    schema = inspect_schema(conn)
    
    import_assembly_batches(conn, schema, force=force)
    import_annotation_batches(conn, schema, force=force)
    
    typer.echo("Database initialisation complete")


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
    """
    Count assembly records matching search criteria.
    """
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
    """
    View the assembly table.
    """
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
    """
    Create a download manifest.
    """
    conn = connect(db_path)
    opts = QueryOptions(taxon=taxon, high_quality=high_quality, has_assembly=has_assembly, max_contamination=max_contamination, min_completeness=min_completeness, has_annotation=has_annotation)
    
    cols, rows = plan_download_impl(conn, opts)
    
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
    concurrency: int = 6,
    override_file_lists: bool = False
):
    """
    Download assembly and annotation data according to a download manifest created with `plan-download`.
    """    
    output_dir.mkdir(parents=True, exist_ok=True)

    async def main():
        mani = Manifest.from_csv(manifest)
        
        # assemblies
        typer.echo("Preparing assembly batches...")
        assembly_list = AssemblyFileLists(output_dir / "atb.assembly.list.csv.gz", override_file_lists)
        batches = mani.batch(AssemblyItem.from_row, assembly_list)
    
        total = await run_jobs_async(batches, output_dir, concurrency)
        typer.echo(f"Extracted {total} assemblies")
    
        if mani.has_annotations():
            typer.echo("Preparing annotation batches...")
            bakta_list = BaktaFileLists(output_dir / "atb.bakta.list.csv.gz", override_file_lists)
            batches = mani.batch(BaktaItem.from_row, bakta_list)
        
            total = await run_jobs_async(batches, output_dir, concurrency)
            typer.echo(f"Extracted {total} annotations")
    
    asyncio.run(main())


def main():
    app()
  

if __name__ == "__main__":
    main()
