
from __future__ import annotations
import logging
from pathlib import Path
from typing import Annotated

from maki.benchmark import bmark_app
from maki.models.enums import Ranks, TaxonomySource
import typer


logger = logging.getLogger(__name__)


app = typer.Typer(help="""
Metagenomics using an Annotated K-mer Index
""")


# ---------------------------------------------------------------------
# Build Workflow Wrapper
# ---------------------------------------------------------------------

CORE_PANEL_NAME = "Core Arguments"

DatabaseArg = Annotated[
    Path,
    typer.Option("-db", "--db-path", help="Database output path", rich_help_panel=CORE_PANEL_NAME)
]
ManifestArg = Annotated[
    Path,
    typer.Option("-i", "--manifest", help="Genome manifest CSV", rich_help_panel=CORE_PANEL_NAME)
]
KmerSizeArg = Annotated[
    int,
    typer.Option("-k", "--kmer-size", help="K-mer size for database", rich_help_panel=CORE_PANEL_NAME)
]
ParallelismArg = Annotated[
    int,
    typer.Option("-p", "--parallel", help="Parallelism parameter.", rich_help_panel=CORE_PANEL_NAME)
]

TAXONOMY_PANEL_NAME = "Taxonomy Arguments"

TaxonomyDatabaseArg = Annotated[
    TaxonomySource,
    typer.Option(help="Taxonomy database", rich_help_panel=TAXONOMY_PANEL_NAME)
]
TaxonomyReleaseArg = Annotated[
    str,
    typer.Option("--release", help="Taxonomy release to use", rich_help_panel=TAXONOMY_PANEL_NAME)
]

CLUSTERING_PANEL_NAME = "Clustering & Pre-filter Arguments"

SketchScaleArg = Annotated[
    int,
    typer.Option(help="Sketching scale parameter", rich_help_panel=CLUSTERING_PANEL_NAME)
]
SketchSeedArg = Annotated[
    int,
    typer.Option(help="Seed for sketching", rich_help_panel=CLUSTERING_PANEL_NAME)
]
TaxonomicRankArg = Annotated[
    Ranks,
    typer.Option(help="Taxonomic rank", rich_help_panel=CLUSTERING_PANEL_NAME)
]
NoRefinementArgs = Annotated[
    bool,
    typer.Option("--no-refinement", help="Skip ANI clustering refinement", rich_help_panel=CLUSTERING_PANEL_NAME)
]
ClusterANIArg = Annotated[
    float,
    typer.Option("-a", "--ani", help="ANI clustering threshold", rich_help_panel=CLUSTERING_PANEL_NAME)
]
NoDerepArgs = Annotated[
    bool,
    typer.Option("--no-derep", help="Skip pre-filter de-replication step", rich_help_panel=CLUSTERING_PANEL_NAME)
]
DereplicationArg = Annotated[
    float,
    typer.Option("--derep", help="De-replication for Sourmash pre-filter index", rich_help_panel=CLUSTERING_PANEL_NAME)
]

BUILD_PANEL_NAME = "Database Build Arguments"

SuffixSizeArg = Annotated[
    int,
    typer.Option("-s", "--suffix-size", help="Index construction suffix size", rich_help_panel=BUILD_PANEL_NAME)
]


@app.command(help="""
Construct a database
""")
def build(
    db_path: DatabaseArg,
    manifest_path: ManifestArg,
    kmer_size: KmerSizeArg,
    parallel: ParallelismArg,
    taxonomy_database: TaxonomyDatabaseArg = TaxonomySource.gtdb,
    taxonomy_release: TaxonomyReleaseArg = "latest",
    scale: SketchScaleArg = 1000,
    seed: SketchSeedArg = 42,
    rank: TaxonomicRankArg = Ranks.Species,
    skip_refine: NoRefinementArgs = False,
    ani: ClusterANIArg = 0.95,
    skip_derep: NoDerepArgs = False,
    derep: DereplicationArg = 0.99,
    suffix_size: SuffixSizeArg = 6,
):
    import maki.workflows.buildstages as build_stages
    
    # Initialise database
    build_stages.init(db_path, kmer_size, taxonomy_database, taxonomy_release,
                      scale, seed, True)
    
    # Design index clusters
    build_stages.add_sketches(db_path, manifest_path, parallel, permissive=True)
    
    build_stages.init_clusters(db_path, rank)
    
    thresholds = []
    if not skip_refine:
        thresholds.append(ani)
        
        if (not skip_derep) and (derep != ani):
            thresholds.append(derep)
        
        rcode = build_stages.refine_clusters(db_path, thresholds, parallel)
        if rcode != 0:
            return rcode
    
    # Index database
    rcode = build_stages.index_clusters(db_path, manifest_path,
                                        f"ani{thresholds[0]}" if thresholds else rank,
                                        suffix_size, parallel)
    return rcode


# ---------------------------------------------------------------------
# Hidden Commands
# ---------------------------------------------------------------------

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
        logger.error("File does not exist: %s", forward)
        return 1

    if not os.path.exists(reverse):
        logger.error("File does not exist: %s", reverse)
        return 1
    
    rp = mx.read_pair(forward, reverse)
    opts = mx.build_opts(k, s, out, threads)
    mx.construct_wdbg(rp, opts)
    
    return 0


app.add_typer(bmark_app, name="benchmark")


def main():
    SystemExit(app())


if __name__ == "__main__":
    main()
