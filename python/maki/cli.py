
from __future__ import annotations
import logging
from pathlib import Path
from typing import Annotated

from maki.benchmark import bmark_app
from maki.models.enums import Ranks, TaxonomySource

import typer


logger = logging.getLogger(__name__)


app = typer.Typer(
    help="Metagenomic database builder and classifier.\n\n"
         "Supports taxonomy-aware clustering, incremental updates, "
         "and parallel classification."
)


@app.command(help="""
Initialise a database and download taxonomy.
""")
def init(
    db_path: Annotated[Path, typer.Option("-db", "--db-path", help="Database output path")],
    kmer_size: Annotated[int, typer.Option("-k", "--kmer-size", help="K-mer size for database")],
    taxonomy_database: Annotated[TaxonomySource, typer.Option(help="Taxonomy database")] = TaxonomySource.gtdb,
    taxonomy_release: str = typer.Option("latest", help="Taxonomy release to use"),
    scale: int = typer.Option(1000, help="Sketching scale parameter"),
    seed: int = typer.Option(42, help="Seed for sketching"),
    permissive: Annotated[bool, typer.Option("--permissive", "-p", help="Do not fail if one already exists")] = False
):
    from maki.models.database import DatabaseHook, DatabaseParameters
    
    DatabaseHook.init(
        db_path,
        DatabaseParameters(
            kmer_size,
            taxonomy_database,
            taxonomy_release,
            scale,
            seed
        ),
        exist_ok=permissive
    )
    
    
@app.command(help="""
Add sketches to database.             
""")
def add_sketches(
    db_path: Annotated[Path, typer.Option("-db", "--db-path", help="Database output path")],
    manifest_path: Annotated[Path, typer.Option("-i", "--manifest", help="Genome manifest CSV")],
    workers: Annotated[int, typer.Option("-w", "--workers", help="Number of Sourmash sketching workers to run in parallel.")]
):
    from maki.models.database import DatabaseHook, GenomeSchema, read_manifest
    
    manifest = read_manifest(manifest_path, GenomeSchema("accession", "taxonomy", "fasta", "gff"))
    logger.info("%d records parsed from manifest %s.", len(manifest), manifest_path)
    
    db = DatabaseHook(db_path)
    sketched = db.sketch_package(manifest, fail_if_complete=True, workers=workers)
    logger.info("Inserted %d sketches into database %s.", len(sketched), db_path)
    

@app.command(help="""
Initialise genome clusters at taxonomic rank.
""")
def init_clusters(
    db_path: Annotated[Path, typer.Option("-db", "--db-path", help="Database output path")],
    rank: Annotated[Ranks, typer.Option(help="Taxonomic rank")] = Ranks.Species,
):
    from maki.models.database import DatabaseHook
    
    db = DatabaseHook(db_path)
    
    if db.clusters is not None:
        logging.error(f"Database %s already has clusters defined.", db_path)
        return 1
    
    db.initialise_clusters_by_rank(rank)


@app.command(help="""
Initialise genome clusters at taxonomic rank.
""")
def clear_clusters(
    db_path: Annotated[Path, typer.Option("-db", "--db-path", help="Database output path")],
):
    from maki.models.database import DatabaseHook
    
    DatabaseHook.clear_clusters(db_path)


@app.command(help="""
Refine clustering at multiple levels
""")
def refine_clusters(
    db_path: Annotated[
        Path,
        typer.Option("-db", "--db-path", help="Database output path")
    ],
    ani: Annotated[
        list[float],
        typer.Option(
            "-a", "--ani",
            help="Comma-separated clustering ANI thresholds"
        )
    ],
    parallel: Annotated[
        int,
        typer.Option("-p", "--parallel", help="Sourmash cores/workers parameter.")
    ],
    dryrun: Annotated[
        bool,
        typer.Option(help="Do not persist new clusters, just report statistics.")
    ] = False
):
    from maki.models.database import DatabaseHook
    
    db = DatabaseHook(db_path)
    
    if db.clusters is None:
        logger.info(f"Database %s has not had clusters initialised. Run `maki init-clusters...`.", db_path)
        return 1
    
    taxa_tag, taxa_cluster_ids = db.ensure_pairwise(parallel)
    
    # Isolate singleton clusters
    singletons = {
        cid
        for cid in taxa_cluster_ids
        if db.clusters.get_cluster(taxa_tag, cid).is_singleton
    }
    
    logger.info(
        "Ignoring %d singleton clusters (%d total) during hierarchical clustering.",
        len(singletons), len(taxa_cluster_ids)
    )
    
    non_singleton_clusters = taxa_cluster_ids - singletons
    
    # Hierarchical clustering on non-singleton clusters
    ani.sort()
    
    hierarchical = db.pairwise.refine_clusters(
        non_singleton_clusters,
        ani,
        parallel
    )
    
    for cutoff, clustering in zip(ani, hierarchical):
        mi, me, ma = clustering.stats()
        
        logger.info(
            "ANI %.3f:\t%d\t%d\t%d\t%d",
            ani, clustering.num_clusters, mi, me, ma
        )

    if not dryrun:
        for cutoff, clustering in zip(ani, hierarchical):
            db.clusters.add_anon_clustering(clustering, f"ani{cutoff}")
        
        db.clusters.save()


@app.command(help="""
Index clusters without
""")
def index_clusters(
    db_path: Annotated[Path, typer.Option("-db", "--db-path", help="Database output path")],
    manifest_path: Annotated[Path, typer.Option("-i", "--manifest", help="Genome manifest CSV")],
    tag: Annotated[list[float], typer.Option("-t", "--tag", help="Clustering to use, identified by tag")],
    parallel: Annotated[int, typer.Option("-p", "--parallel", help="Parallelism (threads)")],
):
    ...


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
