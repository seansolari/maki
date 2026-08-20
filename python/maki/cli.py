from pathlib import Path
from typing import Annotated, Optional
import typer

from maki.models.enums import Ranks, TaxonomySource, UpdateMode
from maki.benchmark.cli import app as bmark_app
from maki.workflows.build import build_app


app = typer.Typer(
    help="Metagenomic database builder and classifier.\n\n"
         "Supports taxonomy-aware clustering, incremental updates, "
         "and parallel classification."
)

app.add_typer(build_app)


@app.command(help="""
Build a new database.

Modes:
- fixed: minimal size, cannot be updated
- updatable: stores compressed source data to allow updates
""")
def build(
    manifest_path: Path = typer.Option(..., help="Genome manifest CSV"),
    db_path: Path = typer.Option(..., help="Database output path"),
    kmer_size: int = typer.Option(31),
    suffix_size: int = typer.Option(6),
    rank: Annotated[Ranks, typer.Option(help="Taxonomic rank")] = Ranks.Species,
    threads: int = typer.Option(4),
    mode: Annotated[UpdateMode, typer.Option(help="Build mode")] = UpdateMode.fixed,
    taxonomy: Annotated[TaxonomySource, typer.Option(help="Taxonomy database")] = TaxonomySource.gtdb,
    releases: Optional[str] = typer.Option(None, help="Taxonomy release to use"),
    force: bool = False,
    dry_run: bool = False
):
    import maki.models.database as mdb
    from .db import build as build_impl
    
    manifest = mdb.read_manifest(manifest_path, mdb.GenomeSchema("accession", "taxonomy", "fasta", "gff"))
    typer.echo(f"{len(manifest)} records parsed from manifest {manifest_path}.")
    
    return build_impl(manifest, db_path, kmer_size, suffix_size, rank, threads, mode, taxonomy, releases, force, dry_run)


# @app.command(help="""
# Update an existing database with new genomes.
# 
# - Performs incremental diff vs current manifest
# - Only rebuilds affected clusters
# - Preserves unchanged indices
# """)
# def update(
#     manifest_path: Path = typer.Option(..., help="New genome manifest"),
#     db_path: Path = typer.Option(..., help="Existing database"),
#     suffix_size: int = typer.Option(6),
#     threads: int = typer.Option(4, help="Parallel threads")
# ):
#     import maki.models.database as mdb
#     from .db import update as update_impl
#     
#     manifest = mdb.read_manifest(manifest_path, mdb.GenomeSchema("accession", "taxonomy", "fasta", "gff"))
#     typer.echo(f"{len(manifest)} records parsed from manifest {manifest_path}.")
#     
#     return update_impl(manifest, db_path, suffix_size, threads)


@app.command(help="""
Classify paired-end sequencing samples.

- Runs filter cluster first
- Classifies reads against genome clusters
- Outputs per-sample results
""")
def classify(
    samples_path: Path = typer.Option(..., help="CSV sample manifest"),
    db_path: Path = typer.Option(..., help="Database path"),
    output: Path = typer.Option(..., help="Output directory"),
    threads: int = typer.Option(4, help="Parallel classification"),
):
    from maki.classify.sample_manifest import SampleManifest
    from maki.classify.classifier import Classifier
    import maki.models.database as mdb
    
    db = mdb.StaticDatabase.load(db_path)
    samples = SampleManifest.from_csv(samples_path)

    classifier = Classifier(db, threads)
    classifier.run(samples, output)


@app.command(hidden=True)
def to_graph(manifest: str, k: int, s: int, threads: int, database_path: Path):
    import maki.core as mx
    
    obj = mx.read_manifest(manifest, mx.FileType.GFF3)
    opts = mx.build_opts(k, s, database_path, threads)
    
    mx.construct_cdbg(obj, opts)
    
    return 0
    

@app.command(hidden=True)
def sample_to_debruijn(forward: str, reverse: str, k: int, s: int, threads: int, out: Path):
    import os
    import maki.core as mx
    
    if not os.path.exists(forward):
        typer.echo(f"Error: file does not exist: {forward}", err=True)
        return 1
    
    if not os.path.exists(reverse):
        typer.echo(f"Error: file does not exist: {reverse}", err=True)
        return 1
    
    rp = mx.read_pair(forward, reverse)
    opts = mx.build_opts(k, s, out, threads)
    mx.construct_wdbg(rp, opts)
    
    return 0


app.add_typer(bmark_app, name="benchmark")


def main():
    app()
    

if __name__ == "__main__":
    main()
