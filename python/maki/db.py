
import logging
from pathlib import Path
from typing import Optional
import maki.models.database as mdb
from maki.models.database.manifest import DatabasePackage


logger = logging.getLogger(__name__)


def build(manifest: DatabasePackage, db_path: Path, kmer_size: int, suffix_size: int, rank: str, threads: int, mode: mdb.UpdateMode, taxonomy: mdb.TaxonomySource, tax_releases: Optional[str] = None, force: bool = False, dry_run: bool = False):
    if db_path.exists() and not force:
        logger.error("Build directory %s already exists.", db_path)
        return 1
  
    opts = mdb.DatabaseOptions(db_path, kmer_size, rank, mode, taxonomy, tax_releases)
    db = mdb.StaticDatabase.create(opts)
    logger.info("Database initialised at %s.", db_path)
    
    db = db.decompress()
    db.insert(manifest, suffix_size, threads, dry_run)
    
    logger.info("Compressing database.")
    db = db.compress()
    logger.info("Database build complete.")
    
    return 0


def update(manifest: DatabasePackage, db_path: Path, suffix_size: int, threads: int):
    db = mdb.StaticDatabase.load(db_path).decompress()
    
    db.insert(manifest, suffix_size, threads)
    db.compress()
    logger.info("Update complete.")
    
    return 0
