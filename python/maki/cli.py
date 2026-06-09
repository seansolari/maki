from pathlib import Path

import typer
import maki.models.database as mdb
from maki.classify.sample_manifest import SampleManifest
from maki.classify.classifier import Classifier

from .db import build as build_impl, update as update_impl


app = typer.Typer(
    help="Metagenomic database builder and classifier.\n\n"
         "Supports taxonomy-aware clustering, incremental updates, "
         "and parallel classification."
)


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
    rank: str = typer.Option(..., help="Taxonomic rank"),
    threads: int = typer.Option(4),
    mode: str = typer.Option(
        "fixed",
        help="Database mode: fixed or updatable"
    ),
    taxonomy: str = typer.Option(
        "gtdb",
        help="Database mode: gtdb or ncbi"
    ),
    force: bool = False
):
    manifest = mdb.read_manifest(manifest_path, mdb.GenomeSchema("accession", "taxonomy", "fasta", "gff"))
    typer.echo(f"{len(manifest)} records parsed from manifest {manifest_path}.")
    
    return build_impl(manifest, db_path, kmer_size, rank, threads, mdb.UpdateMode[mode.lower()], mdb.TaxonomySource[taxonomy.lower()], force)


@app.command(help="""
Update an existing database with new genomes.

- Performs incremental diff vs current manifest
- Only rebuilds affected clusters
- Preserves unchanged indices
""")
def update(
    manifest_path: Path = typer.Option(..., help="New genome manifest"),
    db_path: Path = typer.Option(..., help="Existing database"),
    threads: int = typer.Option(4, help="Parallel threads")
):
    manifest = mdb.read_manifest(manifest_path, mdb.GenomeSchema("accession", "taxonomy", "fasta", "gff"))
    typer.echo(f"{len(manifest)} records parsed from manifest {manifest_path}.")
    
    return update_impl(manifest, db_path, threads)


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
    db = mdb.StaticDatabase.load(db_path)
    samples = SampleManifest.from_csv(samples_path)

    classifier = Classifier(db, threads)
    classifier.run(samples, output)


@app.command(hidden=True)
def to_graph(manifest: str, k: int, s: int, threads: int, database_path: str):
    import maki.core as mx
    manifest = mx.read_manifest(manifest, mx.FileType.GFF3)
    opts = mx.build_opts(k, s, database_path, threads)
    mx.construct_cdbg(manifest, opts)


def main():
    app()
    

if __name__ == "__main__":
    main()
    