from pathlib import Path

import typer
import maki.models.database as mdb
from maki.classify.sample_manifest import SampleManifest
from maki.classify.classifier import Classifier

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
    if db_path.exists() and not force:
        print(f"[error] Build directory {db_path} already exists.")
        return 1
  
    opts = mdb.DatabaseOptions(db_path, kmer_size, rank, mdb.UpdateMode[mode.lower()], mdb.TaxonomySource[taxonomy.lower()])
    db = mdb.StaticDatabase.create(opts)
    typer.echo(f"Database initialised at {db_path}.")
    
    manifest = mdb.read_manifest(manifest_path, mdb.GenomeSchema("accession", "taxonomy", "fasta", "gff"))
    typer.echo(f"{len(manifest)} records parsed from manifest {manifest_path}.")
    
    db = db.decompress()
    db.insert(manifest, threads)
    
    db = db.compress()
    typer.echo("Database build complete.")
    
    return 0


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
    
    db = mdb.StaticDatabase.load(db_path).decompress()
    db.insert(manifest, threads)
    db.compress()
    typer.echo("Update complete.")


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
    confidence: float = typer.Option(0.1, help="Confidence threshold"),
    min_hits: int = typer.Option(5, help="Minimum hits per assignment")
):
    db = mdb.StaticDatabase.load(db_path)
    samples = SampleManifest.from_csv(samples_path)

    # classifier = Classifier(db, threads, confidence, min_hits)
    # classifier.run(samples, output)


def main():
    app()
    

if __name__ == "__main__":
    main()
    