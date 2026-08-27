
import hashlib
from pathlib import Path
import shutil
from typing import Optional

from maki.models.database.clustering.manager import ClusterNode
from maki.utils.hashing import hash_strings

from .index import PreIndex


def _compute_storage_path(tag: str, index_id: str) -> tuple[str, str]:
    """Create a deterministic sharded path with .bz2 extension.
    """
    h = hashlib.sha256((tag + index_id).encode()).hexdigest()

    shard1 = h[:2]
    shard2 = h[2:4]

    relpath = f"tag/{shard1}/{shard2}/{index_id}"

    return h, relpath


def compute_cluster_digest(cluster: ClusterNode):
    return hash_strings(cluster.accessions)


class IndexStoreHandle:
    def __init__(self, root: Path) -> None:
        self.root = root
    
    def current_digest(self, tag: str, index_id: str) -> Optional[str]:
        _, relpath = _compute_storage_path(tag, index_id)
        index_root = self.root / relpath
        
        if not index_root.exists():
            return None
        
        else:
            with (index_root / "digest").open("rt") as f:
                digest = f.read().strip()
            
            return digest
        
    def remove_index(self, tag: str, index_id: str):
        _, relpath = _compute_storage_path(tag, index_id)
        index_root = self.root / relpath
        
        shutil.rmtree(index_root)
        
    def get_handle(self, tag: str, index_id: str):
        _, relpath = _compute_storage_path(tag, index_id)
        index_root = self.root / relpath
        return PreIndex(index_root)
