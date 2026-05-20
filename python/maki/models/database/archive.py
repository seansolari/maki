
from __future__ import annotations
from abc import ABC, abstractmethod
from pathlib import Path
from typing import List

from .cluster import GenomeCluster


class Archive(ABC):
    """
    Abstract archive interface (e.g. tar.xz, zip, etc.)
    """

    @abstractmethod
    def list(self) -> List[str]:
        pass
      
    @abstractmethod
    def get(self, cluster_id: str) -> GenomeCluster:
      pass


class RawArchive(Archive):
    def __init__(self, root: Path) -> None:
        super().__init__()
        
        self.root = root
        self.root.mkdir(parents=True, exist_ok=True)
  
    def list(self) -> List[str]: ...
    
    def get(self, cluster_id: str) -> GenomeCluster: ...
    
    def compress(self) -> XzArchive:
        xz_file = self.root.parent / f"{self.root.name}.tar.xz"
        ...
        return XzArchive(xz_file)


class XzArchive(Archive):
    def __init__(self, xz_file: Path) -> None:
        super().__init__()
        
        self.xz_file = xz_file
        if not self.xz_file.exists():
            self._init_archive()
            
        self.ar = ...
  
    def list(self) -> List[str]: ...
    
    def get(self, cluster_id: str) -> GenomeCluster: ...
    
    def _init_archive(self): ...
    
    def decompress(self) -> RawArchive:
        dir = self.xz_file.parent / self.xz_file.name.removesuffix(".tar.xz")
        if dir.exists():
            ...
        return RawArchive(dir)
