import typer
from maki.models.manifest import Manifest
from maki.models.database import MetaGenomicDatabase
from maki.classify.sample_manifest import SampleManifest
from maki.classify.classifier import Classifier

app = typer.Typer(
    help="Metagenomic database builder and classifier.\n\n"
         "Supports taxonomy-aware clustering, incremental updates, "
         "and parallel classification."
)


@app.command(help="""
Download taxonomy databases.

Supports:
- ncbi: NCBI taxonomy dump
- gtdb: GTDB taxonomy files
""")
def download_taxonomy(
    source: str = typer.Option(..., help="ncbi or gtdb"),
    output: str = typer.Option("./taxonomy", help="Output directory")
):
    from maki.utils.taxonomy_download import download_taxonomy

    download_taxonomy(source, output)

    typer.echo(f"Downloaded {source} taxonomy to {output}")


@app.command(help="""
Build a new database.

Modes:
- fixed: minimal size, cannot be updated
- updatable: stores compressed source data to allow updates
""")
def build(
    manifest_path: str = typer.Option(..., help="Genome manifest CSV"),
    db_path: str = typer.Option(..., help="Database output path"),
    rank: str = typer.Option(..., help="Taxonomic rank"),
    kmer_size: int = typer.Option(31),
    threads: int = typer.Option(4),
    mode: str = typer.Option(
        "fixed",
        help="Database mode: fixed or updatable"
    ),
    force: bool = False
):
    manifest = Manifest.from_csv(manifest_path)
    db = MetaGenomicDatabase(db_path, kmer_size, mode = mode)
    db.build(manifest, rank, threads, force)
    typer.echo("✔ Database build complete.")


@app.command(help="""
Update an existing database with new genomes.

- Performs incremental diff vs current manifest
- Only rebuilds affected clusters
- Preserves unchanged indices
""")
def update(
    manifest_path: str = typer.Option(..., help="New genome manifest"),
    db_path: str = typer.Option(..., help="Existing database"),
    threads: int = typer.Option(4, help="Parallel threads"),
    strategy: str = typer.Option(
        "lazy",
        help="Rebuild strategy: lazy (default) or full"
    ),
    dry_run: bool = typer.Option(
        False,
        help="Show changes without rebuilding"
    )
):
    manifest = Manifest.from_csv(manifest_path)
    db = MetaGenomicDatabase.load(db_path)

    db.update(manifest, threads, strategy, dry_run)
    typer.echo("Update complete.")


@app.command(help="""
Classify paired-end sequencing samples.

- Runs filter cluster first
- Classifies reads against genome clusters
- Outputs per-sample results
""")
def classify(
    samples_path: str = typer.Option(..., help="CSV sample manifest"),
    db_path: str = typer.Option(..., help="Database path"),
    output: str = typer.Option(..., help="Output directory"),
    threads: int = typer.Option(4, help="Parallel classification"),
    confidence: float = typer.Option(0.1, help="Confidence threshold"),
    min_hits: int = typer.Option(5, help="Minimum hits per assignment")
):
    db = MetaGenomicDatabase.load(db_path)
    samples = SampleManifest.from_csv(samples_path)

    classifier = Classifier(db, threads, confidence, min_hits)
    classifier.run(samples, output)


def main():
    app()
    

if __name__ == "__main__":
    main()
    