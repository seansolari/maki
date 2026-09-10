
from __future__ import annotations
from datetime import datetime
import logging
from pathlib import Path
from typing import Annotated, Optional

from maki.models.enums import Ranks, TaxonomySource
import typer


logger = logging.getLogger(__name__)


build_app = typer.Typer()


@build_app.command(help="""
Initialise a database and download taxonomy.
""")
def init(
    db_path: Annotated[Path, typer.Option("-db", "--db-path", help="Database output path")],
    kmer_size: Annotated[int, typer.Option("-k", "--kmer-size", help="K-mer size for database")],
    taxonomy_database: Annotated[TaxonomySource, typer.Option(help="Taxonomy database")] = TaxonomySource.gtdb,
    taxonomy_release: Annotated[Optional[str], typer.Option("--release", help="Taxonomy release to use")] = None,
    scale: Annotated[int, typer.Option(help="Sketching scale parameter")] = 1000,
    seed: Annotated[int, typer.Option(help="Seed for sketching")] = 42,
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
            seed,
            datetime.now().strftime("%Y-%m-%d")
        ),
        exist_ok=permissive
    )
    
    logger.info("Database initialised at %s.", db_path)
    
    return 0
    

@build_app.command(help="""
Add sketches to database.             
""")
def add_sketches(
    db_path: Annotated[Path, typer.Option("-db", "--db-path", help="Database output path")],
    manifest_path: Annotated[Path, typer.Option("-i", "--manifest", help="Genome manifest CSV")],
    workers: Annotated[int, typer.Option("-w", "--workers", help="Number of Sourmash sketching workers to run in parallel.")],
    permissive: Annotated[bool, typer.Option("-p", "--permissive", help="Skip existing sketches")] = False
):
    from maki.models.database import DatabaseHook, GenomeSchema, read_manifest
    
    manifest = read_manifest(manifest_path, GenomeSchema("accession", "taxonomy", "fasta", "gff"))
    logger.info("%d records parsed from manifest %s.", len(manifest), manifest_path)
    
    db = DatabaseHook(db_path)
    sketched = db.sketch_package(manifest, fail_if_complete=not permissive, workers=workers)
    logger.info("Inserted %d sketches into database %s.", len(sketched), db_path)
    
    return 0
    

@build_app.command(help="""
Initialise genome clusters at taxonomic rank.
""")
def init_clusters(
    db_path: Annotated[Path, typer.Option("-db", "--db-path", help="Database output path")],
    rank: Annotated[Ranks, typer.Option(help="Taxonomic rank")] = Ranks.Species,
):
    from maki.models.database import DatabaseHook
    
    db = DatabaseHook(db_path)
    
    if db.clusters is not None:
        logging.warning("Database %s already has clusters defined.", db_path)
        return 0
    
    db.initialise_clusters_by_rank(rank)
    
    return 0


@build_app.command(help="""
Clear all database clusters.
""")
def clear_clusters(
    db_path: Annotated[Path, typer.Option("-db", "--db-path", help="Database output path")],
):
    from maki.models.database import DatabaseHook
    
    DatabaseHook.clear_clusters(db_path)
    
    return 0


@build_app.command(help="""
Refine clustering based on (possibly multiple levels
of) ANI.
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
        logger.error("Database %s has not had clusters initialised. Run `maki init-clusters...`.", db_path)
        return 1
    
    # Ignore values that were already clustered
    
    new_ani = [v for v in ani if not db.clusters.has_tag(f"ani{v}")]
    
    if not new_ani:
        logger.warning("Clusters already exist for all requested thresholds.")
        return 0
    
    if len(new_ani) < len(ani):
        logger.info("Clustering at new thresholds: %s.", ", ".join(map(str, new_ani)))
    
    ani = sorted(new_ani)
    
    # Ensure pairwise data is loaded
    taxa_tag, taxa_cluster_ids = db.ensure_pairwise(parallel)
    
    # Isolate singleton clusters
    singletons = {
        cid
        for cid in taxa_cluster_ids
        if db.clusters.get_cluster(taxa_tag, cid).is_singleton
    }
    
    singleton_accessions = {
        db.clusters.get_cluster(taxa_tag, cid).accessions[0]
        for cid in singletons
    }
    
    logger.info(
        "Ignoring %d singleton clusters (%d total) during hierarchical clustering.",
        len(singletons), len(taxa_cluster_ids)
    )
    
    # Hierarchical clustering on non-singleton clusters
    non_singleton_clusters = taxa_cluster_ids - singletons
    
    hierarchical = db.pairwise.refine_clusters(
        non_singleton_clusters,
        ani,
        parallel
    )
    
    for cutoff, clustering in zip(ani, hierarchical):
        clustering.add_singletons(singleton_accessions)
        
        mi, me, ma = clustering.stats()
        
        logger.info(
            "ANI %.3f:\t%d\t%d\t%d\t%d",
            cutoff, clustering.num_clusters, mi, me, ma
        )

    if not dryrun:
        for cutoff, clustering in zip(ani, hierarchical):
            db.clusters.add_anon_clustering(clustering, f"ani{cutoff}")
        
        db.clusters.save()
    
    return 0


@build_app.command(help="""
Create de Bruijn graph indices based on a taxonomic
or refined clustering.
""")
def index_clusters(
    db_path: Annotated[Path, typer.Option("-db", "--db-path", help="Database output path")],
    manifest_path: Annotated[Path, typer.Option("-i", "--manifest", help="Genome manifest CSV")],
    tag: Annotated[str, typer.Option("-t", "--tag", help="Clustering to use, identified by tag")],
    suffix_size: Annotated[int, typer.Option("-s", "--suffix-size", help="Index construction suffix size")],
    parallel: Annotated[int, typer.Option("-p", "--parallel", help="Parallelism (threads)")],
    dryrun: Annotated[bool, typer.Option(help="Do not perform builds.")] = False,
    force: Annotated[bool, typer.Option(help="Force overwrite of existing clusters (except in `dryrun` mode where nothing happens).")] = False
):
    from maki.models.database import DatabaseHook, GenomeSchema, read_manifest
    
    # Prepare database
    db = DatabaseHook(db_path)
    
    if db.clusters is None:
        logger.error(
            "Database %s has no clusters. Run `maki init-clusters ...` "
            "and then possibly `maki refine-clusters ...`.",
            db_path
        )
        return 1
    
    if tag not in set(db.clusters.tags):
        logger.error(
            "Tag %s not found in database %s. Current clustering tags are: %s.",
            tag,
            db_path,
            ", ".join(db.clusters.tags)
        )
        return 1
    
    # Load manifest
    manifest = read_manifest(
        manifest_path,
        GenomeSchema("accession", "taxonomy", "fasta", "gff")
    )
    logger.info("%d records parsed from manifest %s.", len(manifest), manifest_path)
    
    taxids = {
        accn: taxid
        for (accn, taxid) in db.manifest.iter_taxonomy()
    }
    
    # Build clusters
    for cluster in db.clusters.iter_tag(tag):
        db.index_cluster(
            tag,
            cluster,
            taxids,
            manifest,
            suffix_size,
            parallel,
            dryrun=dryrun,
            force=force
        )
        
    db.indexes.set_params(index_tag=tag)
    
    logger.info("Index construction complete.")
    
    return 0


@build_app.command(
    help="Remove temporary file fragments from failed cluster indexing.",
    hidden=True
)
def clear_failed_indexes(
    db_path: Annotated[Path, typer.Option("-db", "--db-path", help="Database output path")]
):
    from maki.models.database import DatabaseHook
    
    # Prepare database
    db = DatabaseHook(db_path)
    
    if db.clusters is None:
        logger.warning("Database %s has no clusters. Nothing to clean.", db_path)
        return 0
    
    # Clean all tags
    for tag in db.clusters.tags:
        for cluster in db.clusters.iter_tag(tag):
            logger.info(
                "Cleaning any stale files for tag=%s:cluster=%s.",
                tag, cluster.cluster_id
            )
            
            db.indexes.get_handle(tag, cluster.cluster_id).cleanup_temp()

    logger.info("Cleaning complete.")
    
    return 0


@build_app.command(help="""
Construct a reverse index for Sourmash-based pre-filter.
""")
def construct_filter(
    db_path: Annotated[Path, typer.Option("-db", "--db-path", help="Database output path")],
    parallel: Annotated[int, typer.Option("-p", "--parallel", help="Parallelism (threads)")],
    tag: Annotated[Optional[str], typer.Option("-t", "--tag", help="Clustering to use, identified by tag")] = None,
    seed: Annotated[int, typer.Option("--seed, help=De-replication seed")] = 43,
    overwrite: Annotated[bool, typer.Option("--overwrite", help="Overwrite existing filter")] = False
):
    from maki.models.database import DatabaseHook
    
    db = DatabaseHook(db_path)
    
    if db.clusters is None:
        logger.error(
            "Database %s has no clusters. Run `maki init-clusters ...` "
            "and then possibly `maki refine-clusters ...`.",
            db_path
        )
        return 1
    
    if db.indexes.params.index_tag is None:
        logger.error(
            "Database %s has not been indexed. Run `maki index-clusters ...` "
            "before constructing a filter.",
            db_path
        )
        return 1
    
    # Validate tag
    
    if tag == db.indexes.params.index_tag:
        tag = None
        
    if tag:
        if db.indexes.params.index_tag in set(r.value for r in Ranks):
            valid_tags = {
                tag
                for tag in db.clusters.tags
                if tag.startswith("ani")
            }
        
        else:
            assert db.indexes.params.index_tag.startswith("ani")
            
            indexed_ani = float(db.indexes.params.index_tag[3:])
            
            valid_tags = {
                tag
                for tag in db.clusters.tags
                if tag.startswith("ani") and float(tag[3:]) >= indexed_ani
            }
    
        if tag not in valid_tags:
            logger.error(
                "Tag %s not suitable for constructing pre-filter - must be "
                "at least as specific as indexing tag. Valid filter tags "
                "are: %s.",
                tag,
                ",".join(valid_tags)
            )
            return 1
        
    # Construct index
    
    groups = db.clusters.group_by(
        tag=db.indexes.params.index_tag,
        derep_tag=tag,
        seed=seed
    )
    signatures = [
        signature
        for signature in (
            db.sketch_db.load_many([
                accn for accn_list in groups.values() for accn in accn_list
            ])
        ).values()
    ]
    
    db.indexes.rocksdb.construct(
        groups,
        signatures,
        db.sketch_db.params,
        threads=parallel,
        overwrite=overwrite
    )
    db.indexes.set_params(derep_tag=tag)
    
    return 1


def main():
    SystemExit(build_app())


if __name__ == "__main__":
    main()
