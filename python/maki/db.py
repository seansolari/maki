
from pathlib import Path
import maki.models.database as mdb
from maki.models.database.manifest import DatabasePackage


def build(manifest: DatabasePackage, db_path: Path, kmer_size: int, rank: str, threads: int, mode: mdb.UpdateMode, taxonomy: mdb.TaxonomySource, force: bool = False):
    if db_path.exists() and not force:
        print(f"[error] Build directory {db_path} already exists.")
        return 1
  
    opts = mdb.DatabaseOptions(db_path, kmer_size, rank, mode, taxonomy)
    db = mdb.StaticDatabase.create(opts)
    print(f"Database initialised at {db_path}.")
    
    db = db.decompress()
    db.insert(manifest, threads)
    
    print("Compressing database.")
    db = db.compress()
    print("Database build complete.")
    
    return 0


def update(manifest: DatabasePackage, db_path: Path, threads: int):
    db = mdb.StaticDatabase.load(db_path).decompress()
    
    db.insert(manifest, threads)
    db.compress()
    print("Update complete.")
    
    return 0