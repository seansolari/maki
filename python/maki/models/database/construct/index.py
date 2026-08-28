
from __future__ import annotations
from concurrent.futures import ProcessPoolExecutor
import gzip
from itertools import repeat
import logging
import os
from pathlib import Path
import shutil
from typing import Iterable, Optional

from maki.models.database.manifest import GenomeData
from maki.utils.io import open_maybe_gzip
import maki.core as mx


logger = logging.getLogger(__name__)


def _record_worker(args: tuple[Path, GenomeData]):
    source_dir, record = args
    return _try_write_record(source_dir, record)


def _try_write_record(source_dir: Path, record: GenomeData) -> Optional[Path]:
    ext = "gff" if record.gff else "fna"
    dest = source_dir / f"{record.accession}.{ext}.gz"
    if dest.exists():
        return None
    else:
        tmp_dest = dest.with_suffix(".tmp")
        try:
            _write_record_to(record, tmp_dest)
            os.replace(tmp_dest, dest)
        finally:
            tmp_dest.unlink(missing_ok=True)
        return dest


def _write_record_to(record: GenomeData, file: Path):
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
        

class SequenceSourceDir:
    def __init__(self, source_dir: Path) -> None:
        self.source_dir = source_dir
        self.source_dir.mkdir(parents=True, exist_ok=True)
    
    def insert_genome(self, record: GenomeData):
        result = _try_write_record(self.source_dir, record)
        if not result:
            logger.warning("skipping writing %s as it already exists", record.accession)
    
    def insert(self, data: Iterable[GenomeData]):
        for record in data:
            self.insert_genome(record)
    
    def pinsert(self, data: Iterable[GenomeData], concurrency: int):
        tasks = zip(repeat(self.source_dir), data)
        
        with ProcessPoolExecutor(max_workers=concurrency) as executor:
            list(executor.map(_record_worker, tasks))
        
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


class PreIndex(SequenceSourceDir):
    def __init__(self, root: Path) -> None:
        self.root = root
        self.root.parent.mkdir(parents=True, exist_ok=True)
        
        SequenceSourceDir.__init__(
            self,
            self.root.parent / f"{self.root.name}-source"
        )
    
    def clear_sources(self):
        shutil.rmtree(self.source_dir)
        
    def write_digest(self, digest: str):
        if not self.root.exists():
            raise RuntimeError(f"Index does not exist at {self.root}")
        
        (self.root / "digest").write_text(digest)
        
    def _temp_build_file_locs(self):
        tmp_build_dir = self.root.parent / f".{self.root.name}-temp"
        tmp_manifest_file = self.root.parent /  f".{self.root.name}-temp-manifest.txt"
        return tmp_build_dir, tmp_manifest_file
    
    def cleanup_temp(self):
        tmp_build_dir, tmp_manifest_file = self._temp_build_file_locs()
        
        if tmp_build_dir.exists():
            shutil.rmtree(tmp_build_dir)
        
        tmp_manifest_file.unlink(missing_ok=True)
    
    def build(self, k: int, s: int, threads: int):
        tmp_build_dir, tmp_manifest_file = self._temp_build_file_locs()
        tmp_build_dir.mkdir(parents=True)
        
        # Write build manifest
        with tmp_manifest_file.open("wt") as fh:
            self.write_manifest(fh)
        
        # Build
        try:
            manifest = mx.read_manifest(str(tmp_manifest_file), mx.FileType.GFF3)
            opts = mx.build_opts(k, s, tmp_build_dir, threads)
            
            mx.construct_cdbg(manifest, opts)
            
            # Finalise
            os.replace(tmp_build_dir, self.root)

        except Exception:
            self.cleanup_temp()
            
            raise
