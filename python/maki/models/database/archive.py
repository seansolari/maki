
from __future__ import annotations
from pathlib import Path
import shutil
import tarfile
from tempfile import TemporaryDirectory
from typing import Iterator, Set

from .cluster import ClusterHandle, ReadWriteCluster, TemporaryClusterData


class Archive:
    """
    Abstract archive interface (e.g. tar.xz, zip, etc.)
    """
    def __init__(self) -> None:
        self.manifest: Set[str] = set()

    def __iter__(self) -> Iterator[str]:
        return self.manifest.__iter__()
      
    def __contains__(self, qry: str) -> bool:
        return self.manifest.__contains__(qry)


class RawArchive(Archive):
    def __init__(self, root: Path) -> None:
        super().__init__()
        
        self.root = root
        self.root.mkdir(parents=True, exist_ok=True)
        
        self.manifest_file = self.root / "manifest.txt"
        if not self.manifest_file.exists():
            self.manifest_file.touch()
        else:
            with self.manifest_file.open() as f:
                for line in f:
                    self.manifest.add(line.strip())
    
    def get_or_create(self, cluster_id: str) -> ReadWriteCluster:
        self.manifest.add(cluster_id)
        return ReadWriteCluster(self.root / Path(cluster_id))
    
    def compress(self) -> XzArchive:
        xz_file = self.root.parent / f"{self.root.name}.tar.xz"
        
        with tarfile.open(xz_file, "w:xz") as tar:
            tar.add(self.manifest_file, arcname="manifest.txt")
            
            for p in self.manifest:
                tar.add(self.root / Path(p), arcname=p)
        
        shutil.rmtree(self.root)
        
        return XzArchive(xz_file)


class XzArchive(Archive):
    def __init__(self, xz_file: Path) -> None:
        super().__init__()
        
        self.xz_file = xz_file
        if not self.xz_file.exists():
            self._init_archive()
            
        with tarfile.open(self.xz_file, "r:xz") as tar, TemporaryDirectory(dir=self.xz_file.parent) as tmp:
            tmp = Path(tmp)
            tar.extract("manifest.txt", path=tmp)
            
            with (tmp / "manifest.txt").open() as f:
                for line in f:
                    self.manifest.add(line.strip())
    
    def get(self, cluster_id: str) -> ClusterHandle:
        if cluster_id not in self.manifest:
            raise KeyError(f"Unrecognised cluster id: {cluster_id}")
        
        def _extract(p: Path) -> Path:
            with tarfile.open(self.xz_file, "r:xz") as tar:
                tar.extract(cluster_id, p)
            return p / Path(cluster_id)
          
        return ClusterHandle(_extract, dir=self.xz_file.parent)
    
    def _init_archive(self):
        raw = RawArchive(self.root())
        raw.compress()
    
    def root(self) -> Path:
        return self.xz_file.parent / self.xz_file.name.removesuffix(".tar.xz")
    
    def decompress(self) -> RawArchive:
        dir = self.root()
        dir.mkdir(parents=True, exist_ok=True)

        with tarfile.open(self.xz_file, "r:xz") as tar:
            tar.extract("manifest.txt", path=dir)
            
            for p in self.manifest:
                tar.extract(p, path=dir)

        self.xz_file.unlink()
        
        return RawArchive(dir)
