from __future__ import annotations
from copy import deepcopy
from dataclasses import dataclass
from itertools import repeat
import logging
import multiprocessing
import os
import shutil
from tempfile import TemporaryDirectory
import time
from typing import Callable, Dict, Iterable, List, Literal, Optional, Tuple, overload

import csv
import tarfile
from maki.atb.gff import pjson_to_gff
from maki.models.database.manifest import DatabasePackage, GenomeData, GenomeRecord, SequencePackage
from pathlib import Path

import requests
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map

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
          
    def batch(self, dispatch_: Callable[[Dict[str, str]], Tuple[RemoteBatchFile, BatchItem]], prefix_: Path, concurrency: int, desc: str):
        batches: Dict[RemoteBatchFile, List[BatchItemPromise]] = {}
        unidentified_batches: List[Dict[str, str]] = []
        
        for result in process_map(batch_row, self.rows, repeat(dispatch_), repeat(prefix_), max_workers=concurrency, chunksize=100, desc=desc):
            if isinstance(result, dict):
                unidentified_batches.append(result)
            else:
                batch_key, promise = result
                batches.setdefault(batch_key, []).append(promise)
            
        if unidentified_batches:
            cols = list(unidentified_batches[0].keys())
            logger.error(
                "could not resolve the following records:\n%s\n%s\n",
                ",".join(cols),
                "\n".join(",".join(r[c] for c in cols) for r in unidentified_batches)
            )
            raise RuntimeError("Unidentified batch records.")
          
        return batches


def batch_row(r: Dict[str, str], dispatch_: Callable[[Dict[str, str]], Tuple[RemoteBatchFile, BatchItem]], prefix_: Path):
    batch_key, batch_item = dispatch_(r)
                
    if not all((batch_key.tar_xz, batch_key.tar_xz_url, batch_item.filename)):
        return r
    
    target = prefix_ / batch_item.filename
    target_exists = False if not target.exists() else (True if not batch_item.md5 else md5sum(target) == batch_item.md5)
    return batch_key, BatchItemPromise(batch_item, target, target_exists)


@dataclass(frozen=False)
class BatchItem:
    filename: str
    md5: Optional[str] = None
    
    @staticmethod
    def dest(r: Dict[str, str]) -> str: ...
    

@dataclass(frozen=False)
class BatchItemPromise:
    item: BatchItem
    target: Path
    exists: bool
    

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
  

def try_download(url: str, dest_path: Path, urlmd5: str, retries: int = 3, timeout: int = 10):
    """
    Download a file from a URL with retry support and safe cleanup.

    Args:
        url (str): Source URL.
        dest_path (str): Destination file path.
        retries (int): Maximum number of retry attempts.
        timeout (int): Request timeout (seconds).

    Raises:
        Exception: If all retries fail.
    """

    temp_path = dest_path.with_suffix(f"{dest_path.suffix}.part")

    # Always ensure no leftover partial file exists before starting
    if os.path.exists(temp_path):
        os.remove(temp_path)

    for attempt in range(1, retries + 1):
        logger.info("Downloading remote archive %s to %s (attempt=%d).", url, dest_path, attempt)
        
        try:
            with requests.get(url, stream=True, timeout=timeout) as response:
                response.raise_for_status()

                with temp_path.open("wb") as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)

            # Check MD5
            if md5sum(temp_path) != urlmd5:
                logger.error("Download %s failed, invalid MD5 (attempt=%d).", url, attempt)
                raise RuntimeError(f"Download {url} failed, invalid MD5.")

            # Download completed — atomically move into place
            os.replace(temp_path, dest_path)
            return

        except Exception as e:
            # Clean up partial file if something went wrong
            temp_path.unlink(missing_ok=True)

            if attempt == retries:
                raise Exception(
                    f"Failed to download after {retries} attempts: {url}"
                ) from e

            # Small backoff before retrying
            time.sleep(2 ** (attempt - 1))
  

def ensure_archive(batch_rec: RemoteBatchFile, dest: Path):
    if dest.exists():
        if md5sum(dest) != batch_rec.tar_xz_md5:
            logger.error("Corrupted archive %s, removing.", dest)
            os.remove(dest)
        else:
            logger.info("Archive %s already exists with valid MD5. Continuing...", dest)
            return dest

    try_download(batch_rec.tar_xz_url, dest, batch_rec.tar_xz_md5)
    return dest


def download_batch_local(batch_rec: RemoteBatchFile, items: List[BatchItemPromise], output_dir: Path) -> int:      
    queries = {pm.item.filename: pm for pm in items if not pm.exists}
    
    if not queries:
        return 0
    
    local_archive = ensure_archive(batch_rec, output_dir / batch_rec.tar_xz)
    
    extracted: int = 0
    
    with tarfile.open(local_archive, mode="r:*") as tar:
        tar_recnames: List[str] = []
        
        for member in tar:
            tar_recnames.append(member.name)
            
            if member.name not in queries:
                continue

            pm = queries.pop(member.name)
            pm.target.parent.mkdir(parents=True, exist_ok=True)

            extracted_file = tar.extractfile(member)
            if extracted_file is None:
                continue
            
            target_tmp = pm.target.with_suffix(f"{pm.target.suffix}.part")
            try:
                with target_tmp.open("wb") as f:
                    while chunk := extracted_file.read(8192):
                        f.write(chunk)
                        
                if pm.item.md5 and md5sum(target_tmp) != pm.item.md5:
                    raise RuntimeError(f"MD5 mismatch: {pm.target}")
                  
                os.replace(target_tmp, pm.target)
            except Exception:
                logging.warning("Exception raised while writing %s, deleting...", pm.target)
                target_tmp.unlink(missing_ok=True)
                raise

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


def run_jobs_parallel(batches: Dict[RemoteBatchFile, List[BatchItemPromise]], output_dir: Path, concurrency: int, desc: str):
    extracted: int = 0
    
    if not batches:
        logger.info("Nothing to download")
        return extracted

    output_dir.mkdir(parents=True, exist_ok=True)
    
    with multiprocessing.Pool(concurrency) as pool:
      
        with tqdm(total=len(batches), position=0, desc=desc) as pbar:
            for result in pool.imap_unordered(
                download_batch_local_handle,
                ((rec, items, output_dir) for rec, items in batches.items()),
                chunksize=1
            ):
              
                extracted += result
                pbar.update(1)

    return extracted


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
        batches = self.batch(AssemblyItem.from_row, output_dir / "fa", self.concurrency, "Batching FASTA")

        extracted = run_jobs_parallel(batches, output_dir / "fa", self.concurrency, "Downloading FASTA")
        logger.info("Extracted %d assemblies", extracted)

        if self.has_annotations:
            logger.info("Preparing annotation batches...")
            batches = self.batch(BaktaItem.from_row, output_dir / "bakta", self.concurrency, "Batching Bakta")
        
            extracted = run_jobs_parallel(batches, output_dir / "bakta", self.concurrency, "Downloading Bakta")
            logger.info("Extracted %d annotations", extracted)
            
            self._extract_json(batches, output_dir / "gff")
            
    def _extract_json(self, batches: Dict[RemoteBatchFile, List[BatchItemPromise]], base: Path):
        jobs = [
            (pm.target, (base / pm.item.filename).with_suffix(".gff"))
            for promises in batches.values()
            for pm in promises
        ]
        pjson_to_gff(jobs, self.concurrency)


class LightAtbManifest(_MakiAtbManifest, DatabasePackage):
    def retrieve_data(self, accessions: Iterable[str]):
        return self._retrieve_data_dispatch(accessions, True)


class DiskAtbManifest(_MakiAtbManifest, DatabasePackage):
    def retrieve_data(self, accessions: Iterable[str]):
        return self._retrieve_data_dispatch(accessions, False)
