
from pathlib import Path
import shutil
import tarfile
from tempfile import NamedTemporaryFile, TemporaryDirectory
from typing import Optional

from maki.models.database.manifest import GenomeData, SequencePackage
from maki.utils.io import open_maybe_gzip
import maki.core as mx


class Cluster:
    def __init__(self, root: Path) -> None:
        self.root = root
        self.index_dir = self.root / "index"
        self.index_dir.mkdir(parents=True, exist_ok=True)


class TemporaryClusterData:
    def __init__(self, dir: Optional[Path] = None):
        self._fh = TemporaryDirectory(dir=dir)
        
        self.root: Optional[Path] = None
        
    @property
    def dir(self) -> Path:
        return Path(self._fh.name)
    
    def set_root(self, rel_path: Path):
        self.root = self.dir / rel_path
        
    def __del__(self):
        self._fh.cleanup()


class ReadOnlyCluster(Cluster):
    def __init__(self, handle: TemporaryClusterData):
        if not handle.root:
            raise ValueError("TemporaryClusterData not properly initialised.")
      
        super().__init__(handle.root)
        self._dh = handle


class ReadWriteCluster(Cluster):
    def __init__(self, root: Path) -> None:
        super().__init__(root)
        self.source_dir = self.root / "source"
        self._changed = False
        self._prepare_source_dir()
        
    def __del__(self):
        if self.updateable and self._changed:
            self.persist_sources()

        shutil.rmtree(self.source_dir)
    
    @property
    def xz_file(self) -> Path:
        return self.source_dir.with_suffix(".tar.xz")
      
    @property
    def updateable(self) -> bool:
        return self.xz_file.exists()
      
    # ======================
    # SOURCE MANAGEMENT
    # ======================
        
    def insert(self, data: SequencePackage):
        for record in data.records():
            trg = self.source_dir / f"{record.accession}.{"gff" if record.gff else "fa"}"
            if not trg.exists():
                self._write_record_to(record, trg)
                self._changed = True
            else:
                print(f"[warning] skipping writing {trg} as it already exists")
    
    def persist_sources(self):
        with tarfile.open(self.xz_file, "w:xz") as tar:
            for file in filter(lambda p: p.is_file(), self.source_dir.iterdir()):
                tar.add(file, arcname=file.name)
          
        self._changed = False
        
    def remove_sources(self):
        self.xz_file.unlink(missing_ok=True)
        
    # ======================
    # INDEX
    # ======================
        
    def build(self, k: int, threads: int):
        with NamedTemporaryFile(suffix=".txt", dir=self.root) as fh:
            for p in self.source_dir.iterdir():
                fh.write(f"{p}\n".encode("utf-8"))
            fh.flush()
          
            manifest = mx.read_manifest(fh.name, 0, ",", mx.FileType.GFF3)
            opts = mx.build_opts(k, 7, self.index_dir, threads)
            
            mx.construct_cdbg(manifest, opts)
    
    # ======================
    # DETAILS
    # ======================
        
    def _prepare_source_dir(self):
        self.source_dir.mkdir(parents=True, exist_ok=True)
        
        if self.xz_file.exists():
            with tarfile.open(self.xz_file, "r:xz") as tar:
                tar.extractall(self.source_dir)

    def _write_record_to(self, record: GenomeData, file: Path):
        with file.open("w") as f:
            if record.gff:
                with open_maybe_gzip(Path(record.gff)) as gff:
                    for line in gff:
                        if "\tbakta\t" in line:
                            f.write(f"#{line}")
                        else:
                            f.write(line)
                f.write("##FASTA\n")
            
            with open_maybe_gzip(Path(record.fasta)) as fna:
                f.writelines(fna)
                