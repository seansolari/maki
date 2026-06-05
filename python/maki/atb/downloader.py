from __future__ import annotations
from abc import ABC, abstractmethod
from copy import deepcopy
from dataclasses import dataclass
import logging
import multiprocessing
import os
import shutil
from tempfile import TemporaryDirectory
from typing import Callable, Dict, Iterable, List, Literal, Optional, Tuple, overload
from urllib.request import urlretrieve

import csv
import tarfile
from maki.atb.gff import pjson_to_gff
from maki.models.database.manifest import DatabasePackage, GenomeData, GenomeRecord, SequencePackage
from pathlib import Path

from tqdm import tqdm

from .models import ASSEMBLY_SCHEMA, ASSEMBLY_BATCH_SCHEMA, BAKTA_BATCH_SCHEMA, BAKTA_SCHEMA
from ..utils.io import open_maybe_gzip, md5sum, write_maybe_gzip
from .lists import RemoteBatchFile


logger = logging.getLogger(__name__)


class Manifest:
    def __init__(self, columns: Iterable[str]):
        self.columns = list(columns)
        self.rows: List[Dict[str, str]] = []
        
    def __iter__(self):
        return self.rows.__iter__()
      
    def __len__(self):
        return self.rows.__len__()
      
    def records(self):
        for row in self:
            yield GenomeRecord(row[ASSEMBLY_SCHEMA.sample], row[ASSEMBLY_SCHEMA.species], None)
        
    def insert(self, row: Dict[str, str]):
        assert all(c in row for c in self.columns), f"Row missing columns: {", ".join(c for c in self.columns if c not in row)}"
        self.rows.append(row)
        
    def append(self, rhs: Manifest):
        assert self.columns == rhs.columns
        self.rows.extend(rhs.rows)
        
    @property
    def has_assemblies(self) -> bool:
        return all(c in self.columns for c in (ASSEMBLY_SCHEMA.tar_xz, ASSEMBLY_BATCH_SCHEMA.tar_xz_url, ASSEMBLY_BATCH_SCHEMA.tar_xz_md5, ASSEMBLY_BATCH_SCHEMA.tar_xz_size_MB, ASSEMBLY_SCHEMA.sample))
        
    @property
    def has_annotations(self) -> bool:
        return all(c in self.columns for c in (BAKTA_SCHEMA.tar_xz, BAKTA_BATCH_SCHEMA.tar_xz_url, BAKTA_BATCH_SCHEMA.tar_xz_md5, BAKTA_BATCH_SCHEMA.tar_xz_size_MB, BAKTA_SCHEMA.file_name, BAKTA_SCHEMA.file_md5))
        
    @classmethod
    def from_csv(cls, path: Path, *args, **kwargs):
        with open_maybe_gzip(path) as f:
            fieldnames: Optional[List[str]] = None
            
            # Parse metadata
            for line in f:
                if not line.startswith("#"):
                    fieldnames = line.strip().split(",")
                    break
            
            assert fieldnames, f"Empty csv: {path}"
            
            reader = csv.DictReader(f, fieldnames=fieldnames)
            assert reader.fieldnames
            
            mf = cls(reader.fieldnames, *args, **kwargs)
            for r in reader:
                mf.insert(r)

            return mf
    
    def write_csv(self, file: Path) -> None:
        columns = deepcopy(self.columns)
        columns.append("fa")
        if self.has_annotations:
            columns.append("gff")
        
        def row_iter():
            do_annot = self.has_annotations
            for row in self.rows:
                res = [row[c] for c in self.columns]
                res.append(str(file.parent / "fa" / AssemblyItem.dest(row)))
                if do_annot:
                    res.append(str((file.parent / "gff" / BaktaItem.dest(row)).with_suffix(".gff")))
                yield res
        
        with write_maybe_gzip(file) as f:
            writer = csv.writer(f)
            writer.writerow(columns)
            writer.writerows(row_iter())
          
    def batch(self, dispatch: Callable[[Dict[str, str], ], Tuple[RemoteBatchFile, BatchItem]]):
        batches: Dict[RemoteBatchFile, List[BatchItem]] = {}
        unidentified_batches: List[Dict[str, str]] = []
        
        for r in self.rows:
            batch_key, batch_item = dispatch(r)
            
            if not all((batch_key.tar_xz, batch_key.tar_xz_url, batch_item.filename)):
                unidentified_batches.append(r)
                continue
            
            batches.setdefault(batch_key, []).append(batch_item)
                
        if unidentified_batches:
            cols = list(unidentified_batches[0].keys())
            logger.error(
                "could not resolve the following records:\n%s\n%s\n",
                ",".join(cols),
                "\n".join(",".join(r[c] for c in cols) for r in unidentified_batches)
            )
            raise RuntimeError("Unidentified batch records.")
          
        return batches


@dataclass(frozen=False)
class BatchItem:
    filename: str
    md5: Optional[str] = None
    
    @staticmethod
    def dest(r: Dict[str, str]) -> str: ...
    

class AssemblyItem(BatchItem):
    @staticmethod
    def dest(r: Dict[str, str]) -> str:
        return f"{r[ASSEMBLY_SCHEMA.tar_xz].removesuffix(".tar.xz")}/{r[ASSEMBLY_SCHEMA.sample]}.fa"
  
    @classmethod
    def from_row(cls, r: Dict[str, str]) -> Tuple[RemoteBatchFile, AssemblyItem]:
        return RemoteBatchFile(r[ASSEMBLY_SCHEMA.tar_xz], r[ASSEMBLY_BATCH_SCHEMA.tar_xz_url], r[ASSEMBLY_BATCH_SCHEMA.tar_xz_md5], r[ASSEMBLY_BATCH_SCHEMA.tar_xz_size_MB]), AssemblyItem(filename=cls.dest(r))


class BaktaItem(BatchItem):
    @staticmethod
    def dest(r: Dict[str, str]) -> str:
          return f"{r[BAKTA_SCHEMA.tar_xz].removesuffix(".tar.xz")}/{r[BAKTA_SCHEMA.file_name]}"
    
    @classmethod
    def from_row(cls, r: Dict[str, str]) -> Tuple[RemoteBatchFile, BaktaItem]:
        return RemoteBatchFile(r[BAKTA_SCHEMA.tar_xz], r[BAKTA_BATCH_SCHEMA.tar_xz_url], r[BAKTA_BATCH_SCHEMA.tar_xz_md5], r[BAKTA_BATCH_SCHEMA.tar_xz_size_MB]), BaktaItem(filename=cls.dest(r), md5=r[BAKTA_SCHEMA.file_md5])


def missing_batch_items(items: List[BatchItem], output_dir: Path):
    queries: Dict[str, BatchItem] = {}
        
    for item in items:
        target = output_dir / item.filename
        
        if target.exists():
            if item.md5:
                if md5sum(target) == item.md5:
                    continue
            else:
                continue
        
        queries[item.filename] = item

    return queries
  
  
def ensure_archive(batch_rec: RemoteBatchFile, dest: Path):
    if dest.exists():
        if md5sum(dest) != batch_rec.tar_xz_md5:
            logger.error("Corrupted archive %s, removing.", dest)
            os.remove(dest)
        else:
            logger.info("Archive %s already exists with valid MD5. Continuing...", dest)
            return dest

    logger.info("Downloading remote archive %s to %s.", batch_rec.tar_xz_url, dest)
    urlretrieve(batch_rec.tar_xz_url, dest)
    
    if md5sum(dest) != batch_rec.tar_xz_md5:
        logger.error("Download %s failed, invalid MD5.", batch_rec.tar_xz_url)
        raise RuntimeError(f"Download {batch_rec.tar_xz_url} failed, invalid MD5.")
    
    return dest


def download_batch_local(batch_rec: RemoteBatchFile, items: List[BatchItem], output_dir: Path):
    queries = missing_batch_items(items, output_dir)
    
    if not queries:
        return 0
    
    local_archive = ensure_archive(batch_rec, output_dir / batch_rec.tar_xz)
    extracted = 0
    
    with tarfile.open(local_archive, mode="r:*") as tar:
        tar_recnames: List[str] = []
        
        for member in tar:
            tar_recnames.append(member.name)
            
            if member.name not in queries:
                continue

            item = queries.pop(member.name)
                    
            target = output_dir / item.filename
            target.parent.mkdir(parents=True, exist_ok=True)

            extracted_file = tar.extractfile(member)
            if extracted_file is None:
                continue
            
            try:
                with target.open("wb") as f:
                    while chunk := extracted_file.read(8192):
                        f.write(chunk)
                        
                if item.md5 and md5sum(target) != item.md5:
                    raise RuntimeError(f"MD5 mismatch: {target}")
            except Exception as e:
                logging.warning("Exception raised while writing %s, deleting...", target)
                target.unlink(missing_ok=True)
                raise e

            extracted += 1
    
    
    # check for unresolved queries
    if queries:
        logger.error("Could not find requested files in batch %s @ %s: %s", batch_rec.tar_xz, batch_rec.tar_xz_url, ", ".join(queries.keys()))
        logger.error("Example keys are: %s...", ", ".join(tar_recnames[:10]))
        raise RuntimeError(f"Could not find requested files in batch {batch_rec.tar_xz}: {", ".join(queries.keys())}")
    else:
        logger.info("Completed extraction of files from archive %s, deleting...", local_archive)
        os.remove(local_archive)
      
    return extracted


def download_batch_local_handle(args):
    rec, items, output_dir = args
    return download_batch_local(rec, items, output_dir)


def run_jobs_parallel(batches: Dict[RemoteBatchFile, List[BatchItem]], output_dir: Path, concurrency: int = 6):
    if not batches:
        logger.info("Nothing to download")
        return 0

    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Pre-allocate worker slots
    total = 0
    
    with multiprocessing.Pool(concurrency) as pool:
      
        with tqdm(total=len(batches), position=0, desc="Batch") as pbar:
            for result in pool.imap_unordered(
                download_batch_local_handle,
                ((rec, items, output_dir) for rec, items in batches.items()),
                chunksize=1
            ):
              
                total += result
                pbar.update(1)

    return total


class _MakiAtbData(SequencePackage):
    def __init__(self) -> None:
        super().__init__()
        self._genomes: Dict[str, GenomeData] = {}
        
    def insert(self, rec: GenomeData):
        self._genomes[rec.accession] = rec
        
    def __getitem__(self, *args, **kwargs):
        return self._genomes.__getitem__(*args, **kwargs)
        
    def genomes(self):
        return iter(self._genomes.values())

        
class TempMakiAtbData(_MakiAtbData):
    def __init__(self, dir: Optional[Path] = None) -> None:
        super().__init__()
        self._thd = TemporaryDirectory(dir=dir)

    @property
    def path(self):
        return Path(self._thd.name)
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self._thd.cleanup()


class PersistentMakiAtbData(_MakiAtbData):
    """Only cleanup on graceful exit."""
    
    def __init__(self, dir: Optional[Path] = None) -> None:
        super().__init__()
        self._path = (dir or Path.cwd()) / ".maki"
        self._path.mkdir(parents=True, exist_ok=True)

    @property
    def path(self):
        return self._path
    
    def __exit__(self, *args):
        if all(a is None for a in args):
            self._rm()
    
    def _rm(self):
        shutil.rmtree(self._path)


class _MakiAtbManifest(Manifest):
    def __init__(self, columns: Iterable[str], concurrency: int, base: Path):
        super().__init__(columns)
        self.concurrency = concurrency
        self.base = base

    def empty(self):
        return type(self)(self.columns, self.concurrency, self.base)
      
    @overload
    def _retrieve_data_dispatch(self, accessions: Iterable[str], temp: Literal[True]) -> TempMakiAtbData: ...
    
    @overload
    def _retrieve_data_dispatch(self, accessions: Iterable[str], temp: Literal[False]) -> PersistentMakiAtbData: ...
      
    def _retrieve_data_dispatch(self, accessions: Iterable[str], temp: bool):
        # retrieve query rows
        submanifest = self.select(accessions)
        
        # download data
        result = TempMakiAtbData(self.base) if temp else PersistentMakiAtbData(self.base)
        submanifest.download_batches(result.path)
        
        do_annot = self.has_annotations
        for row in self:
            result.insert(GenomeData(row[ASSEMBLY_SCHEMA.sample], row[ASSEMBLY_SCHEMA.species], None, str(result.path / "fa" / AssemblyItem.dest(row)), None if not do_annot else str((result.path / "gff" / BaktaItem.dest(row)).with_suffix(".gff"))))
        
        return result
    
    def select(self, accessions: Iterable[str]):
        submanifest = self.empty()
        queries = set(accessions)
        
        for row in self:
            if row[ASSEMBLY_SCHEMA.sample] in queries:
                submanifest.insert(row)
        
        assert len(submanifest) == len(queries)
        return submanifest
  
    def download_batches(self, output_dir: Path):
        # assemblies
        assert self.has_assemblies, "Missing assembly metadata"
        logger.info("Preparing assembly batches...")
        batches = self.batch(AssemblyItem.from_row)

        total = run_jobs_parallel(batches, output_dir / "fa", self.concurrency)
        logger.info("Extracted %d assemblies", total)

        if self.has_annotations:
            logger.info("Preparing annotation batches...")
            batches = self.batch(BaktaItem.from_row)
        
            total = run_jobs_parallel(batches, output_dir / "bakta", self.concurrency)
            logger.info("Extracted %d annotations", total)
            
            pjson_to_gff([f.filename for files in batches.values() for f in files], output_dir, self.concurrency)


class LightAtbManifest(_MakiAtbManifest, DatabasePackage):
    def retrieve_data(self, accessions: Iterable[str]):
        return self._retrieve_data_dispatch(accessions, True)


class DiskAtbManifest(_MakiAtbManifest, DatabasePackage):
    def retrieve_data(self, accessions: Iterable[str]):
        return self._retrieve_data_dispatch(accessions, False)
