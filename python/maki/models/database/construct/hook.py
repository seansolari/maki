
from __future__ import annotations
from dataclasses import asdict, dataclass
import json
import logging
import os
from pathlib import Path
from typing import Optional
import uuid

from maki.models.database.manifest import DatabasePackage
from maki.models.database.sketch import SourmashSketchStore, SketchParameters
from maki.models.enums import TaxonomySource
from maki.models.taxonomy import BaseTaxonomy, GTDBTaxonomy, NCBITaxonomy, resolve_accession_taxids
from maki.models.taxonomy.gtdb import GTDBRelease
from .manifest import SQLiteManifest


logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class DatabaseParameters:
    # Index options
    kmer_size: int
    
    # Taxonomy options
    taxonomy_database: TaxonomySource
    taxonomy_release: str
    
    # Sketching options
    sketch_scale: int
    sketch_seed: int
    
    @property
    def sketch_params(self):
        return SketchParameters(
            self.kmer_size,
            self.sketch_scale,
            self.sketch_seed
        )
    
    @classmethod
    def from_dict(cls, data: dict):
        data["taxonomy_database"] = TaxonomySource[data["taxonomy_database"]]
        return cls(**data)

    def __post_init__(self) -> None:
        if self.kmer_size <= 0:
            raise ValueError("kmer_size must be greater than zero")

        if self.sketch_scale <= 0:
            raise ValueError("sketch_scale must be greater than zero")

        if not self.taxonomy_database:
            raise ValueError("taxonomy_database must not be empty")

        if not self.taxonomy_release:
            raise ValueError("taxonomy_release must not be empty")


class DatabaseHook:
    """
    Lightweight interface to a metagenomic database preparation directory.

    Directory layout:

        database-root/
            metadata.json
            manifest.sqlite3
            data/
                GCF/
                    GCF_000001405.40/
                        ...
            working/
    """

    METADATA_FILE = "metadata.json"
    MANIFEST_FILE = "manifest.sqlite3"
    DATA_DIRECTORY = "data"
    TAXONOMY_DIRECTORY = "taxonomy"
    SKETCH_DIRECTORY = "sketch"
    
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
        self.data_root = self.root / self.DATA_DIRECTORY
        self.taxonomy_root = self.root / self.TAXONOMY_DIRECTORY
        self.sketch_root = self.data_root / self.SKETCH_DIRECTORY

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
        self.taxonomy = self._ensure_taxonomy(
            self.taxonomy_root,
            self.metadata.taxonomy_database,
            self.metadata.taxonomy_release
        )
        
        # Sketch records
        self.sketch_db = SourmashSketchStore(self.sketch_root, params=self.metadata.sketch_params)

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
        fail_if_exists: bool = False,
        workers: int | None = None
    ) -> tuple[str, ...]:
        """
        Reserve, retrieve, publish, and finalize a package.

        Returns the accessions inserted by this call. When fail_if_exists is
        False, accessions already present in the manifest are skipped.
        """
        
        # Reserve accessions
        
        records = self._assign_taxids(package)

        package_id, reserved_accessions = self.manifest.reserve(
            records,
            fail_if_exists=fail_if_exists,
        )

        if not reserved_accessions:
            return ()
        
        sketched_accessions = self.sketch_db.sketch_many(
            package,
            reserved_accessions,
            workers=workers
        )
        
        self.manifest.mark_complete(
            sketched_accessions,
            package_id=package_id,
        )

        return sketched_accessions
    
    # Manifest I/O
    # ------------

    def export_taxonomy_tsv(
        self,
        destination: str | Path | None = None,
        *,
        include_header: bool = False,
        overwrite: bool = False,
    ) -> Path:
        """
        Export mappings for complete accessions.

        By default this creates database-root/taxonomy.tsv.
        """
        if destination is None:
            destination = self.root / "taxonomy.tsv"

        return self.manifest.export_taxonomy_tsv(
            destination,
            complete_only=True,
            include_header=include_header,
            overwrite=overwrite,
        )

    def export_manifest_tsv(
        self,
        destination: str | Path | None = None,
        *,
        complete_only: bool = False,
        overwrite: bool = False,
    ) -> Path:
        """
        Export the full manifest, including processing state and audit fields.
        """
        if destination is None:
            destination = self.root / "manifest.tsv"

        return self.manifest.export_manifest_tsv(
            destination,
            complete_only=complete_only,
            include_header=True,
            overwrite=overwrite,
        )
                
    # Taxonomy helpers
    # ----------------

    @staticmethod
    def _ensure_taxonomy(root: Path, source: TaxonomySource, release_str: Optional[str] = None) -> BaseTaxonomy:
        if source == TaxonomySource.ncbi:
            taxonomy = NCBITaxonomy(root)
        else:
            versions = []
            
            if release_str:
                release_str = release_str.strip()
                if ',' in release_str:
                    requested_releases = release_str.split(',')
                elif ' ' in release_str:
                    requested_releases = release_str.split()
                else:
                    requested_releases = [release_str]
                
                for request in requested_releases:
                    try:
                        versions.append(GTDBRelease[request])
                    except KeyError:
                        err = f"Unrecognised taxonomy release ({request}), available values are: {','.join(m.name for m in GTDBRelease)}"
                        logger.error(err)
                        raise ValueError(err)
            
            taxonomy = GTDBTaxonomy(root, releases=versions or None)
        
        return taxonomy
    
    def _assign_taxids(self, package: DatabasePackage):
        taxids = resolve_accession_taxids(
            ((r.accession, r.taxonomy) for r in package.records() if not r.taxid),
            self.taxonomy
        )
        
        return [
            (r.accession, r.taxid or taxids[r.accession])
            for r in package.records()
        ]
