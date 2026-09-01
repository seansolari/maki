from __future__ import annotations
from abc import ABC, abstractmethod
import hashlib
from pathlib import Path
from typing import Iterable

from .core import SketchParameters


class SourmashSketchStore(ABC):
    """
    Large-scale sourmash signature manager.

    Features
    --------
    * SHA256 sharded filesystem layout
    * gzip-compressed signatures
    * SQLite index
    * multiprocessing batch sketch generation
    * random signature loading
    """

    def __init__(
        self,
        root_dir: str | Path,
        *,
        params: SketchParameters
    ):
        self.root = Path(root_dir)
        self.params = params
        self.params.validate()
        
    @staticmethod
    @abstractmethod
    def signature_name(accession: str) -> str: ...

    # --------------------------------------------------
    # sketch helpers
    # --------------------------------------------------
    
    @classmethod
    def _compute_storage_path(cls, accession: str) -> str:
        """Create a deterministic sharded path.
        """
        h = hashlib.sha256(accession.encode()).hexdigest()

        shard1 = h[:2]
        shard2 = h[2:4]

        relpath = f"{shard1}/{shard2}/{cls.signature_name(accession)}"

        return relpath
    
    def filter_existing(self, accessions: Iterable[str]):
        return {
            accession
            for accession in accessions
            if not self._sketch_exists(accession)
        }
    
    def _sketch_exists(self, accession: str):
        return (self.root / self._compute_storage_path(accession)).exists()
    
    def cleanup_accession(self, accession: str):
        sigpath = self.root / self._compute_storage_path(accession)
        sigpath.unlink(missing_ok=True)
        sigpath.with_suffix(sigpath.suffix + ".tmp").unlink(missing_ok=True)
    
    def cleanup_accessions(self, accessions: Iterable[str]):
        for accession in accessions:
            self.cleanup_accession(accession)
