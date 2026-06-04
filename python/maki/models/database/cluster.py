
from pathlib import Path
import shutil
import tarfile
from tempfile import NamedTemporaryFile, TemporaryDirectory
from typing import Callable, Optional

from maki.models.database.manifest import GenomeData, SequencePackage
from maki.utils.io import open_maybe_gzip
import maki.core as mx


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
        result = self._try_write_record(record)
        if result:
            self._changed = True
        else:
            print(f"[warning] skipping writing {record.accession} as it already exists")
    
    def insert(self, data: SequencePackage):
        for record in data.genomes():
            self.insert_genome(record)
    
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
        
    def write_manifest(self, fh):
        for p in self.source_dir.iterdir():
            fh.write(f"{p}\n".encode("utf-8"))
        fh.flush()
        
    def _prepare_source_dir(self):
        self.source_dir.mkdir(parents=True, exist_ok=True)
        
        if self.xz_file.exists():
            with tarfile.open(self.xz_file, "r:xz") as tar:
                tar.extractall(self.source_dir)
        
    def _try_write_record(self, record: GenomeData) -> Optional[Path]:
        ext = "gff" if record.gff else "fna"
        dest = self.source_dir / f"{record.accession}.{ext}"
        if dest.exists():
            return None
        else:
            self._write_record_to(record, dest)
            return dest

    def _write_record_to(self, record: GenomeData, file: Path):
        with file.open("w") as f:
            # Write GFF
            if record.gff:
                with open_maybe_gzip(Path(record.gff)) as gff:
                    for line in gff:
                        if "\tbakta\tregion\t" in line:
                            f.write(f"#{line}")
                        elif line.startswith("##FASTA"):
                            if record.fasta:
                                print(f"[warning] record {record.accession} supplies sequence in both GFF and FNA, preferring GFF.")
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
                with open_maybe_gzip(Path(record.fasta)) as fna:
                    f.writelines(fna)


class ReadWriteCluster(Cluster, SequenceSourceDir):
    def __init__(self, root: Path) -> None:
        Cluster.__init__(self, root)
        SequenceSourceDir.__init__(self, self.root / "source")
        
    def __del__(self):
        if self.updateable and self.changed:
            self.persist_sources()

        shutil.rmtree(self.source_dir)
    
    def build(self, k: int, threads: int):
        with NamedTemporaryFile(suffix=".txt", dir=self.root) as fh:
            self.write_manifest(fh)
          
            manifest = mx.read_manifest(fh.name, 0, ",", mx.FileType.GFF3)
            opts = mx.build_opts(k, 7, self.index_dir, threads)
            
            mx.construct_cdbg(manifest, opts)

                