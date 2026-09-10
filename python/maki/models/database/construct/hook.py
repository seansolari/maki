
from __future__ import annotations
from dataclasses import asdict, dataclass
import json
import logging
import os
from pathlib import Path
from typing import Optional
import uuid

from maki.models.database.clustering import ClusterManager, ClusterNode, PairwiseManager
from maki.models.database.manifest import DatabasePackage
from maki.sketch import SketchParameters, pairwise_ani_comparison
from maki.models.enums import Ranks, TaxonomySource
from maki.models.taxonomy import BaseTaxonomy, GTDBTaxonomy, NCBITaxonomy, resolve_accession_taxids

from .annotations import AnnotationIndexBuilder, SeedTaxonomy
from .manifest import SQLiteManifest
from .sketch import GenomeSketchStore
from .store import IndexStoreHandle, compute_cluster_digest


logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class DatabaseParameters:
    # Index options
    kmer_size: int
    
    # Taxonomy options
    taxonomy_database: TaxonomySource
    taxonomy_release: Optional[str]
    
    # Sketching options
    sketch_scale: int
    sketch_seed: int
    
    # Metadata
    creation_date: str
    
    @property
    def sketch_params(self):
        return SketchParameters(
            self.kmer_size,
            self.sketch_scale,
            self.sketch_seed
        )
    
    @classmethod
    def from_dict(cls, **kwargs):
        kwargs["taxonomy_database"] = TaxonomySource[kwargs["taxonomy_database"]]
        return cls(**kwargs)

    def __post_init__(self) -> None:
        if self.kmer_size <= 0:
            raise ValueError("kmer_size must be greater than zero")

        if self.sketch_scale <= 0:
            raise ValueError("sketch_scale must be greater than zero")

        if not self.taxonomy_database:
            raise ValueError("taxonomy_database must not be empty")


class DatabaseHook:
    """Lightweight interface to a metagenomic database preparation directory.
    """

    METADATA_FILE = "metadata.json"
    MANIFEST_FILE = "manifest.sqlite3"
    CLUSTER_FILE = "clusters.json.bz2"
    PAIRWISE_FILE = "pairwise.json.bz2"
    DATA_DIRECTORY = "data"
    TAXONOMY_DIRECTORY = "taxonomy"
    SKETCH_DIRECTORY = "sketch"
    INDEX_DIRECTORY = "index"
    
    # Database creation and hook initialisation
    # -----------------------------------------

    def __init__(self, database_root: str | Path) -> None:
        self.root = Path(database_root)

        if not self.root.is_dir():
            raise NotADirectoryError(
                f"Database root does not exist: {self.root}"
            )

        self.metadata_path = self.root / self.METADATA_FILE
        self.manifest_path = self.root / self.MANIFEST_FILE
        self.cluster_file = self.root / self.CLUSTER_FILE
        self.pairwise_file = self.root / self.PAIRWISE_FILE
        self.data_root = self.root / self.DATA_DIRECTORY
        self.taxonomy_root = self.root / self.TAXONOMY_DIRECTORY
        
        self.sketch_root = self.data_root / self.SKETCH_DIRECTORY
        self.index_root = self.data_root / self.INDEX_DIRECTORY

        if not self.metadata_path.is_file():
            raise FileNotFoundError(
                f"Database metadata does not exist: {self.metadata_path}"
            )

        if not self.data_root.is_dir():
            raise FileNotFoundError(
                f"Database data directory does not exist: {self.data_root}"
            )
            
        if not self.taxonomy_root.is_dir():
            raise FileNotFoundError(
                "Database taxonomy directory does not exist: "
                f"{self.taxonomy_root}"
            )
            
        if not self.sketch_root.is_dir():
            raise FileNotFoundError(
                "Database sketch directory does not exist: "
                f"{self.sketch_root}"
            )
            
        if not self.index_root.is_dir():
            raise FileNotFoundError(
                "Database index directory does not exist: ",
                f"{self.index_root}"
            )

        with self.metadata_path.open("r", encoding="utf-8") as handle:
            metadata_data = json.load(handle)

        try:
            self.metadata = DatabaseParameters.from_dict(**metadata_data)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Invalid database metadata in {self.metadata_path}"
            ) from exc

        # Record metadata
        self.manifest = SQLiteManifest(self.manifest_path)
        self.taxonomy = self._load_taxonomy(self.taxonomy_root, self.metadata.taxonomy_database)
        
        # Sketch records
        self.sketch_db = GenomeSketchStore(self.sketch_root, params=self.metadata.sketch_params)
        
        # Clusters
        if self.cluster_file.exists():
            self.clusters = ClusterManager(self.cluster_file)
        else:
            self.clusters = None
        
        # Pairwise distances within root clusters
        self.pairwise = PairwiseManager(self.pairwise_file)
        
        # Index handle
        self.indexes = IndexStoreHandle(self.index_root)

    @classmethod
    def init(
        cls,
        database_root: str | Path,
        metadata: DatabaseParameters,
        *,
        exist_ok: bool = False,
    ) -> DatabaseHook:
        """
        Atomically initialize a new database root.

        If exist_ok is True, an existing, already initialized database is
        opened. Its metadata is not silently replaced.
        """
        root = Path(database_root)

        if root.exists():
            if not root.is_dir():
                raise NotADirectoryError(root)

            if not exist_ok:
                raise FileExistsError(root)

            hook = cls(root)

            if hook.metadata != metadata:
                raise ValueError(
                    "Existing database metadata does not match the "
                    "requested metadata"
                )

            return hook

        # Build the database beside the requested path and atomically rename it.
        root.parent.mkdir(parents=True, exist_ok=True)

        staging_root = root.with_name(
            f".{root.name}.initializing-{uuid.uuid4().hex}"
        )

        try:
            staging_root.mkdir()
            (staging_root / cls.DATA_DIRECTORY).mkdir()
            (staging_root / cls.DATA_DIRECTORY / cls.SKETCH_DIRECTORY).mkdir()
            (staging_root / cls.DATA_DIRECTORY / cls.INDEX_DIRECTORY).mkdir()

            metadata_path = staging_root / cls.METADATA_FILE
            temporary_metadata = staging_root / f".{cls.METADATA_FILE}.tmp"

            with temporary_metadata.open("w", encoding="utf-8") as handle:
                json.dump(
                    asdict(metadata),
                    handle,
                    indent=2,
                    sort_keys=True,
                )
                handle.write("\n")
                handle.flush()
                os.fsync(handle.fileno())

            os.replace(temporary_metadata, metadata_path)

            SQLiteManifest.create(
                staging_root / cls.MANIFEST_FILE
            )
            
            cls._ensure_taxonomy(staging_root / cls.TAXONOMY_DIRECTORY, metadata.taxonomy_database, metadata.taxonomy_release)

            # This succeeds atomically only if another process has not already
            # created root.
            staging_root.rename(root)

        except Exception:
            if staging_root.exists():
                import shutil

                shutil.rmtree(staging_root, ignore_errors=True)
            raise

        return cls(root)
    
    # Insert data packages
    # --------------------

    def sketch_package(
        self,
        package: DatabasePackage,
        *,
        fail_if_complete: bool = False,
        workers: int | None = None
    ) -> tuple[str, ...]:
        """
        Reserve, retrieve, publish, and finalize a package.

        Returns the accessions inserted by this call. When fail_if_complete is
        False, accessions already sketched are skipped.
        """
        
        # Reserve accessions
        
        records = self._assign_taxids(package)

        plan = self.manifest.reserve(
            records,
            fail_if_complete=fail_if_complete,
        )

        if not plan:
            return ()
        
        # Clean old failed accessions
        if plan.cleanup_accessions:
            try:
                self.sketch_db.cleanup_accessions(plan.cleanup_accessions)
                
                self.manifest.finish_cleanup(
                    plan.cleanup_accessions,
                    package_id=plan.package_id
                )
                
            except Exception as exc:
                cleanup_error = f"{type(exc).__name__}: {exc}"
                
                try:
                    self.manifest.fail_cleanup(
                        plan.cleanup_accessions,
                        package_id=plan.package_id,
                        error=cleanup_error,
                    )
                    
                except Exception as manifest_exc:
                    exc.add_note(
                        "Additionally failed to restore cleanup entries to "
                        f"'failed': {manifest_exc}"
                    )
                    
                # New accessions were reserved in the same initial transaction.
                # Since this package will not proceed, mark those as failed too.
                try:
                    self.manifest.mark_failed(
                        plan.reserved_accessions,
                        package_id=plan.package_id,
                        error=(
                            "Package aborted because cleanup of another accession "
                            f"failed: {cleanup_error}"
                        ),
                    )
                    
                except Exception as manifest_exc:
                    exc.add_note(
                        "Additionally failed to mark new reservations as failed: "
                        f"{manifest_exc}"
                    )
                    
                raise
        
        # At this point, both new entries and cleaned entries are reserved and
        # owned by this package.
        insertion_accessions = plan.all_accessions
        
        try:
            sketched_accessions = self.sketch_db.sketch_many(
                package,
                insertion_accessions,
                workers=workers
            )
            
            self.manifest.mark_complete(
                sketched_accessions,
                package_id=plan.package_id,
            )
            
        except Exception as exc:
            
            try:
                self.manifest.mark_failed(
                    insertion_accessions,
                    package_id=plan.package_id,
                    error=f"{type(exc).__name__}: {exc}",
                )
            except Exception as manifest_exc:
                exc.add_note(
                    "Additionally failed to mark the insertion as failed: "
                    f"{manifest_exc}"
                )
            
            raise

        return sketched_accessions
    
    # Cluster by taxonomy
    # -------------------
    
    def initialise_clusters_by_rank(self, rank: str):
        self.clusters = ClusterManager(self.cluster_file)
        
        # Refine by rank

        data = list(self.manifest.iter_taxonomy(complete_only=True))
        ancestors = self._group_by_rank(data, rank)
        
        self.clusters.add_named_clustering(
            ancestors,
            tag=rank
        )
        
        self.clusters.save()
        
    # Taxonomy helpers
    # ----------------
    
    @staticmethod
    def _load_taxonomy(root: Path, source: TaxonomySource) -> BaseTaxonomy:
        if source == TaxonomySource.ncbi:
            return NCBITaxonomy(root)
        else:
            return GTDBTaxonomy(root)

    @staticmethod
    def _ensure_taxonomy(root: Path, source: TaxonomySource, release_str: Optional[str] = None) -> None:
        if source == TaxonomySource.ncbi:
            NCBITaxonomy.ensure_release(root, release_str)
        else:
            GTDBTaxonomy.ensure_release(root, release_str)
    
    def _assign_taxids(self, package: DatabasePackage):
        taxids = resolve_accession_taxids(
            ((r.accession, r.taxonomy) for r in package.records() if not r.taxid),
            self.taxonomy
        )
        
        return [
            (r.accession, r.taxid or taxids[r.accession])
            for r in package.records()
        ]
        
    def _group_by_rank(self, recs: list[tuple[str, str]], rank: str):
        ancestors: dict[str, str] = {}
        
        for accn, taxid in recs:
            ancestor = self.taxonomy.get_ancestor_at_rank(taxid, rank)
            ancestors[accn] = ancestor
        
        return ancestors
    
    # Pairwise comparison within taxonomy
    # -----------------------------------
    
    def ensure_pairwise(self, workers: int):
        assert isinstance(self.clusters, ClusterManager)
        
        self.pairwise.try_load()
        
        # Get base taxonomic clusters
        taxa_tag = self._extract_taxa_tag()
        
        taxa_clusters = {
            c.cluster_id
            for c in self.clusters.iter_tag(taxa_tag)
        }
        
        # Identify taxa that have not undergone sketch-based pairwise comparison
        existing_pairwise = set(self.pairwise.group_ids())
        to_do = taxa_clusters - existing_pairwise
        
        if to_do:
            _changed = False
            
            for cluster_id in to_do:
                cluster = self.clusters.get_cluster(taxa_tag, cluster_id)
                
                # No need to compare singletons
                if not cluster.is_singleton:
                    sigs = self.sketch_db.load_many(cluster.accessions)
                    
                    pw = pairwise_ani_comparison(
                        sigs.values(),
                        self.sketch_db.params,
                        self.root,
                        workers
                    )
                    
                    self.pairwise.insert_results(cluster_id, pw)
                    _changed = True
            
            if _changed:
                self.pairwise.save()
                
        return taxa_tag, set(taxa_clusters)
        
    def _extract_taxa_tag(self):
        assert isinstance(self.clusters, ClusterManager)
        
        candidates = set(r.value for r in Ranks) & set(self.clusters.tags)
        
        if not candidates:
            raise RuntimeError("Taxonomic clustering tag could not be identified.")
        elif len(candidates) > 1:
            raise RuntimeError("Multiple taxonomic clustering tags identified.")
        
        return candidates.pop()
    
    # Regression
    # ----------
    
    @staticmethod
    def clear_clusters(db_path: Path):
        (db_path / DatabaseHook.CLUSTER_FILE).unlink(missing_ok=True)
        (db_path / DatabaseHook.PAIRWISE_FILE).unlink(missing_ok=True)
        
    # Build indexes
    # -------------
    
    def index_cluster(
        self,
        tag: str,
        cluster: ClusterNode,
        taxids: dict[str, str],
        package: DatabasePackage,
        suffix_size: int,
        parallel: int,
        dryrun: bool = False,
        force: bool = False
    ):  
        # Compare data with existing cluster
        this_digest = compute_cluster_digest(cluster)
        existing_digest = self.indexes.current_digest(tag, cluster.cluster_id)
        
        if existing_digest:
            if this_digest == existing_digest:
                logger.info("Skipping build for tag=%s:cluster=%s, identical accession digests.", tag, cluster.cluster_id)
                return
            
            elif force and not dryrun:
                logger.info(
                    "Removing existing index for for tag="
                    "%s:cluster=%s, stale digest.",
                    tag,
                    cluster.cluster_id
                )
                
                self.indexes.remove_index(tag, cluster.cluster_id)
                
        # Start cluster build
        logger.info("Starting index build for tag=%s:cluster=%s.", tag, cluster.cluster_id)

        if dryrun:
            logger.info(
                "(dryrun) %d accessions to insert for tag=%s:cluster=%s: %s%s",
                len(cluster.accessions),
                tag,
                cluster.cluster_id,
                ",".join(cluster.accessions[:5]),
                "..." if len(cluster.accessions) > 5 else ""
            )
            return
        
        with package.retrieve_data(cluster.accessions) as data:
            dbh = self.indexes.get_handle(tag, cluster.cluster_id)
            
            logger.info("Inserting sequences tag=%s:cluster=%s.", tag, cluster.cluster_id)
            
            # Process sequence data to disk
            annots = dbh.pinsert(data.genomes(), parallel)
            
            # Save annotation and taxonomy metadata
            with AnnotationIndexBuilder(dbh.sqlite_file) as table:
                table.add_annotations(annots)
                
                table.update_taxonomy(
                    SeedTaxonomy(accn, taxids[accn]) for accn in cluster.accessions
                )
                
                table.finalize()
            
            # Construct index
            try:
                logger.info("Constructing cluster index for tag=%s:cluster=%s.", tag, cluster.cluster_id)
                
                dbh.build(
                    self.metadata.kmer_size,
                    suffix_size,
                    parallel
                )
                
                dbh.write_digest(this_digest)
            
            finally:
                logger.info("Clearing sources for cluster index for tag=%s:cluster=%s.", tag, cluster.cluster_id)
                dbh.clear_sources()
                
            logger.info("Index construction complete for tag=%s:cluster=%s.", tag, cluster.cluster_id)
