
from dataclasses import asdict, dataclass
import hashlib
import json
import logging
import os
from pathlib import Path
import shutil
from typing import Optional

from maki.models.database.clustering.manager import ClusterNode
from maki.sketch.core import build_rocksdb_index
from maki.utils.hashing import hash_strings
from sourmash import SourmashSignature

from .index import PreIndex


logger = logging.getLogger(__name__)


def _compute_storage_path(tag: str, index_id: str) -> tuple[str, str]:
    """Create a deterministic sharded path with .bz2 extension.
    """
    h = hashlib.sha256((tag + index_id).encode()).hexdigest()

    shard1 = h[:2]
    shard2 = h[2:4]

    relpath = f"{tag}/{shard1}/{shard2}/{index_id}"

    return h, relpath


def compute_cluster_digest(cluster: ClusterNode):
    return hash_strings(cluster.accessions)


@dataclass(slots=True)
class IndexParameters:
    index_tag: Optional[str] = None
    derep_tag: Optional[str] = None
    
    
class ReverseIndexHandle:
    def __init__(self, prefix: Path) -> None:
        self.index_file = prefix.parent / f"{prefix.name}.rocksdb"
        self.manifest_file = prefix.parent / f"{prefix.name}-manifest.txt"
        
    def construct(
        self,
        groups: dict[str, list[str]],
        signatures: list[SourmashSignature],
        *,
        threads: int,
        overwrite: bool = False
    ):
        if self.index_file.exists():
            if overwrite:
                logger.info("Removing existing index at %s", self.index_file)
                self.index_file.unlink(missing_ok=True)
                self.manifest_file.unlink(missing_ok=True)

            else:
                raise RuntimeError(
                    "Index already exists at %s",
                    self.index_file
                )

        # Staging files
        tmp_index = self.index_file.with_suffix(".tmprocksdb")
        tmp_manifest = self.manifest_file.with_suffix(".tmp")
        
        build_rocksdb_index(signatures, tmp_index, threads)
        self._dump_groups(groups, tmp_manifest)
        
        tmp_manifest.replace(self.manifest_file)
        tmp_index.replace(self.index_file)
    
    @staticmethod
    def _dump_groups(
        groups: dict[str, list[str]],
        outfile: Path
    ):
        with outfile.open("w") as f:
            for cluster_id, accession_list in groups.items():
                for accession in accession_list:
                    f.write(f"{accession}\t{cluster_id}\n")


class IndexStoreHandle:
    
    METADATA_FILE = "parameters.json"
    
    def __init__(self, root: Path) -> None:
        self.root = root
        
        meta_file = self.root / self.METADATA_FILE
        if meta_file.exists():
            with meta_file.open("r") as fh:
                param_dict = json.load(fh)
            
            self.params = IndexParameters(**param_dict)
            
        else:
            self.params = IndexParameters()
            
        self.rocksdb = ReverseIndexHandle(self.root / "reverse-index")
            
    def set_params(
        self,
        *,
        index_tag: Optional[str] = None,
        derep_tag: Optional[str] = None
    ):
        if index_tag or derep_tag:
            if index_tag:
                self.params.index_tag = index_tag
            
            if derep_tag:
                self.params.derep_tag = derep_tag
                
            self._dump_params()
    
    def _dump_params(self):
        meta_file = self.root / self.METADATA_FILE
        tmp_meta = meta_file.with_suffix(".tmp")
        
        with tmp_meta.open("w", encoding="utf-8") as fh:
            json.dump(asdict(self.params), fh, indent=2)
        
        os.replace(tmp_meta, meta_file)
    
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
