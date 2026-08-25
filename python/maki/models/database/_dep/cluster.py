
from __future__ import annotations
import gzip
from itertools import repeat
import logging
import os
from pathlib import Path
import shutil
import tarfile
from tempfile import NamedTemporaryFile, TemporaryDirectory
from typing import Callable, Iterable, Optional

from tqdm.contrib.concurrent import process_map
from maki.models.database.manifest import GenomeData
from maki.utils.io import open_maybe_gzip
import maki.core as mx


logger = logging.getLogger(__name__)


def try_write_record(source_dir: Path, record: GenomeData) -> Optional[Path]:
    ext = "gff" if record.gff else "fna"
    dest = source_dir / f"{record.accession}.{ext}.gz"
    if dest.exists():
        return None
    else:
        tmp_dest = dest.with_suffix(".tmp")
        try:
            write_record_to(record, tmp_dest)
            os.replace(tmp_dest, dest)
        finally:
            tmp_dest.unlink(missing_ok=True)
        return dest


def write_record_to(record: GenomeData, file: Path):
    with gzip.open(file, "wt") as f:
        # Write GFF
        if record.gff:
            with open_maybe_gzip(Path(record.gff), "rt") as gff:
                for line in gff:
                    if "\tbakta\tregion\t" in line:
                        f.write(f"#{line}")
                    elif line.startswith("##FASTA"):
                        if record.fasta:
                            logger.warning("Record supplies %s sequence in both GFF and FNA, preferring GFF.", record.accession)
                        f.write(line)
                        f.writelines(gff)
                        return
                    else:
                        f.write(line)
            f.write("##FASTA\n")
        
        # Write FNA
        if not record.fasta:
            raise RuntimeError(f"No FASTA sequence supplied for {record.accession}")
        else:
            with open_maybe_gzip(Path(record.fasta), "rt") as fna:
                f.writelines(fna)


class Cluster:
    def __init__(self, root: Path) -> None:
        self.root = root
        self.index_dir = self.root / "index"
        self.index_dir.mkdir(parents=True, exist_ok=True)


class XzClusterIndex(Cluster):
    def __init__(self, root: Path) -> None:
        super().__init__(root)
        
        self._db = mx.cdbg()
        mx.cdbg.disk_load(self._db, self.index_dir)


class ClusterHandle:
    def __init__(self, _initter: Callable[[Path], Path], dir: Optional[Path] = None) -> None:
        if dir:
            dir.mkdir(parents=True, exist_ok=True)
        
        self._tmp = TemporaryDirectory(dir=dir)
        self.cluster = XzClusterIndex(_initter(Path(self._tmp.name)))
        
    def __del__(self):
        self._tmp.cleanup()
        
    @property
    def name(self) -> str:
        return self.cluster.root.name
      
    @property
    def db(self) -> mx.cdbg:
        return self.cluster._db
        

class SequenceSourceDir:
    def __init__(self, source_dir: Path) -> None:
        self.source_dir = source_dir
        
        self._changed = False
        self._prepare_source_dir()
      
    def insert_genome(self, record: GenomeData):
        result = try_write_record(self.source_dir, record)
        if result:
            self._changed = True
        else:
            logger.warning("skipping writing %s as it already exists", record.accession)
    
    def insert(self, data: Iterable[GenomeData]):
        for record in data:
            self.insert_genome(record)
    
    def pinsert(self, data: Iterable[GenomeData], concurrency: int):
        for result in process_map(try_write_record, repeat(self.source_dir), data, max_workers=concurrency, chunksize=1, desc=f"Insert@{self.source_dir.name}"):
            if result:
                self._changed = True
    
    @property
    def xz_file(self) -> Path:
        return self.source_dir.with_suffix(".tar.xz")
      
    @property
    def updateable(self) -> bool:
        return self.xz_file.exists()
      
    @property
    def changed(self) -> bool:
        return self._changed
      
    def persist_sources(self):
        with tarfile.open(self.xz_file, "w:xz") as tar:
            for file in filter(lambda p: p.is_file(), self.source_dir.iterdir()):
                tar.add(file, arcname=file.name)
          
        self._changed = False
        
    def remove_sources(self):
        self.xz_file.unlink(missing_ok=True)
        
    def sources(self):
        for p in self.source_dir.iterdir():
            if p.name.endswith(".fna.gz") or p.name.endswith(".gff.gz"):
                yield p
        
    def accessions(self):
        for p in self.sources():
            yield p.name.rsplit(".", 2)[0]
        
    def write_manifest(self, fh):
        for p in self.sources():
            fh.write(f"{p}\n")
        fh.flush()
        
    def _prepare_source_dir(self):
        self.source_dir.mkdir(parents=True, exist_ok=True)
        
        if self.xz_file.exists():
            with tarfile.open(self.xz_file, "r:xz") as tar:
                tar.extractall(self.source_dir)


class ReadWriteCluster(Cluster, SequenceSourceDir):
    def __init__(self, root: Path) -> None:
        Cluster.__init__(self, root)
        SequenceSourceDir.__init__(self, self.root / "source")
        
    def __del__(self):
        if self.updateable and self.changed:
            self.persist_sources()

        shutil.rmtree(self.source_dir)
    
    def build(self, k: int, s: int, threads: int):
        with NamedTemporaryFile(suffix=".txt", dir=self.root, mode="wt") as fh:
            self.write_manifest(fh)
          
            manifest = mx.read_manifest(fh.name, mx.FileType.GFF3)
            opts = mx.build_opts(k, s, self.index_dir, threads)
            
            mx.construct_cdbg(manifest, opts)

                