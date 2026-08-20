from __future__ import annotations
from pathlib import Path

import typer


build_app = typer.Typer()


@build_app.command(help="""
Cluster input genome sequences.
""")
def cluster(
    manifest_path: Path = typer.Option(..., help="Genome manifest CSV"),
    output: Path = typer.Option(..., help="Clustering results prefix"),
    kmer_size: int = typer.Option(31),
    scaled: int = typer.Option(2000),
    seed: int = typer.Option(42),
    ani_threshold: float = typer.Option(0.95),
    lsh_hashes: int = typer.Option(64),
    band_size: int = typer.Option(4),
    max_candidates: int = typer.Option(128),
    processes: int = typer.Option(32)
):
    from maki.models.database import read_manifest, GenomeSchema
    from maki.models.database.sketch.cluster import (
        GenomeClusterBuilder,
        SketchParameters,
    )

    builder = GenomeClusterBuilder(
        output_dir=output,
        sketch_params=SketchParameters(
            ksize=kmer_size,
            scaled=scaled,
            seed=seed,
        ),
        ani_threshold=ani_threshold,
        lsh_hashes=lsh_hashes,
        band_size=band_size,
        max_candidates=max_candidates,
        processes=processes,
    )
    
    manifest = read_manifest(manifest_path, GenomeSchema("accession", "taxonomy", "fasta", "gff"))
    typer.echo(f"{len(manifest)} records parsed from manifest {manifest_path}.")
    
    records = builder.build(manifest)

    print(f"Clustered {len(records)} genomes")
    print(f"Clusters: {len({r.cluster_id for r in records})}")


@build_app.command(help="""
Construct sequence index.
""")
def index_clusters(
    db_path: Path = typer.Option(..., help="Database output path"),
    cluster_path: Path = typer.Option(..., help="Path to clustering results"),
    suffix_size: int = typer.Option(6),
    parallel: int = typer.Option(4),
    delete_cluster_data: bool = typer.Option(True),
    force: bool = False,
    dry_run: bool = False
):
    from maki.models.database.sketch.query import ClusterContainmentDatabase

    database = ClusterContainmentDatabase.build(
        phase1_dir=cluster_path,
        database_dir=db_path / "sketches",
        processes=parallel,
        delete_genome_sketches=delete_cluster_data,
        sbt_threshold=0.01,
    )
    
    # continue build
    ...
