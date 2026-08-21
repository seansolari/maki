from __future__ import annotations
import json
from pathlib import Path

import typer


build_app = typer.Typer()


@build_app.command(help="""
Add input sequences to sketch.
""")
def sketch(
    manifest_path: Path = typer.Option(..., help="Genome manifest CSV"),
    output: Path = typer.Option(..., help="Sketch base"),
    kmer_size: int = typer.Option(31),
    scaled: int = typer.Option(2000),
    seed: int = typer.Option(42),
    processes: int = typer.Option(32),
    skip_existing: bool = typer.Option(False, help="Skip existing records")
):
    from maki.models.database import read_manifest, GenomeSchema
    from maki.models.database.sketch.cluster import (
        SketchParameters, SourmashSketchStore
    )
    
    manifest = read_manifest(manifest_path, GenomeSchema("accession", "taxonomy", "fasta", "gff"))
    typer.echo(f"{len(manifest)} records parsed from manifest {manifest_path}.")
    
    store = SourmashSketchStore(
        output,
        params=SketchParameters(
            ksize=kmer_size,
            scaled=scaled,
            seed=seed,
        )
    )
    
    store.sketch_many(manifest, workers=processes, skip_existing=skip_existing)
    

@build_app.command(help="""
Cluster input genome sequences.
""")
def cluster(
    sketch_dir: Path = typer.Option(..., help="Path to sketch directory"),
    ani_threshold: float = typer.Option(0.95, help="Cluster similarity boundary."),
    downsample_factor: int = typer.Option(10),
    max_candidates: int = typer.Option(128),
):
    from maki.models.database.sketch.cluster import (
        GenomeClusterBuilder,
        SourmashSketchStore
    )

    builder = GenomeClusterBuilder(
        ani_threshold=ani_threshold,
        downsample_factor=downsample_factor,
        max_candidates=max_candidates,
    )
    
    sketches = SourmashSketchStore(sketch_dir)
    clusters, metadata = builder.cluster(sketches)
    
    with (sketch_dir / "clusters.csv").open("wt") as f:
        for cx in clusters:
            f.write(f"{cx.genome_id},{cx.cluster_id}\n")
    
    with (sketch_dir / "clustering-metadata.json").open("wt") as f:
        json.dump(metadata, f)


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
