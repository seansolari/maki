from __future__ import annotations
from copy import deepcopy
from dataclasses import dataclass
import logging
import multiprocessing
from tempfile import TemporaryDirectory
import time
from typing import Callable, Dict, Iterable, List, Optional, Tuple

import csv
import tarfile
from maki.atb.gff import pjson_to_gff
from maki.models.database.manifest import DatabasePackage, GenomeData, GenomeRecord, SequencePackage
import requests
from requests.adapters import HTTPAdapter
from pathlib import Path

from tqdm import tqdm
import urllib3

from .models import ASSEMBLY_SCHEMA, ASSEMBLY_BATCH_SCHEMA, BAKTA_BATCH_SCHEMA, BAKTA_SCHEMA
from ..utils.io import open_maybe_gzip, md5sum, write_maybe_gzip
from .lists import RemoteBatchFile


logger = logging.getLogger(__name__)
_worker_id: Optional[int] = None


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


def init_downloader(lock, worker_ids: List[int]):
    tqdm.set_lock(lock)
    global _worker_id
    _worker_id = worker_ids.pop()


class TqdmStream:
    def __init__(self, raw, *args, **kwargs):
        self.raw = raw
        self.pbar = tqdm(*args, **kwargs)
        
    def __enter__(self):
        return self
      
    def __exit__(self, *args, **kwargs):
        self.pbar.__exit__(*args, **kwargs)

    def read(self, size=-1, *args, **kwargs):
        data = self.raw.read(size, *args, **kwargs)
        self.pbar.update(len(data))
        return data

    def readable(self):
        return True
      
    def write(self, b: bytes, *args, **kwargs):
        return self.raw.write(b, *args, **kwargs)
  
    def tell(self, *args, **kwargs) -> int:
        return self.raw.tell(*args, **kwargs)
      
    def seek(self, pos: int, *args, **kwargs):
        return self.raw.seek(pos, *args, **kwargs)
      
    def close(self, *args, **kwargs):
        return self.raw.close(*args, **kwargs)
      
    @property
    def name(self) -> str | bytes:
        return self.raw.name
      
    @property
    def mode(self):
        return self.raw.mode


def download_batch(batch_rec: RemoteBatchFile, items: List[BatchItem], output_dir: Path, max_retries: int = 3):
    global _worker_id
    assert _worker_id is not None, "Worker ID not set"
    
    extracted = 0

    # ignore already existing results
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

    if queries:
        session = requests.Session()
        retries = urllib3.Retry(total=max_retries, connect=max_retries, read=max_retries, backoff_factor=1, allowed_methods=None)
        session.mount('http://', HTTPAdapter(max_retries=retries))
        session.mount('https://', HTTPAdapter(max_retries=retries))
      
        with session.get(batch_rec.tar_xz_url, stream=True) as r:
            r.raise_for_status()
            total_size = int(r.headers.get('content-length', 0))
            
            with TqdmStream(r.raw, total=total_size, unit="B", unit_scale=True, desc=batch_rec.tar_xz_url, position=_worker_id + 1, leave=False) as stream,\
                tarfile.open(fileobj=stream, mode="r|*") as tar:
                
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
                    
                    with target.open("wb") as f:
                        while chunk := extracted_file.read(8192):
                            f.write(chunk)

                    # validate
                    if item.md5 and md5sum(target) != item.md5:
                        raise RuntimeError(f"MD5 mismatch: {target}")

                    extracted += 1

                # check for unresolved queries
                if queries:
                    logger.error("Could not find requested files in batch %s @ %s: %s", batch_rec.tar_xz, batch_rec.tar_xz_url, ", ".join(queries.keys()))
                    logger.error("Example keys are: %s...", ", ".join(tar_recnames[:10]))
                    raise RuntimeError(f"Could not find requested files in batch {batch_rec.tar_xz}: {", ".join(queries.keys())}")
    
    return extracted

  
def download_batch_resilient(args, max_tries: int = 3):
    batch_rec, items, output_dir = args
    
    last_exception: Optional[Exception] = None
    
    for attempt in range(1, max_tries + 1):
        try:
            return download_batch(batch_rec, items, output_dir, max_tries)
        except (requests.RequestException,
                urllib3.exceptions.ProtocolError,
                IOError) as e:
            
            last_exception = e
            
            logger.warning("Download failed (attempt %d/%d) for %s: %s", attempt, max_tries, batch_rec.tar_xz_url, e)
            
            if attempt < max_tries:
                time.sleep(2 ** attempt)
    
    logger.error("Download permanently failed after %d attempts: %s", max_tries, batch_rec.tar_xz_url)
    assert last_exception
    raise last_exception


def run_jobs_parallel(batches: Dict[RemoteBatchFile, List[BatchItem]], output_dir: Path, concurrency: int = 6):
    if not batches:
        logger.info("Nothing to download")
        return 0

    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Pre-allocate worker slots
    manager_list = list(range(concurrency))
    total = 0
    
    with multiprocessing.Pool(
        concurrency,
        initializer=init_downloader,
        initargs=(multiprocessing.RLock(), manager_list)
    ) as pool:
      
        with tqdm(total=len(batches), position=0, desc="Batch") as pbar:
            for result in pool.imap_unordered(
                download_batch_resilient,
                ((rec, items, output_dir) for rec, items in batches.items()),
                chunksize=1
            ):
              
                total += result
                pbar.update(1)

    return total


class MakiAtbData(SequencePackage):
    def __init__(self, dir: Optional[Path] = None) -> None:
        super().__init__()
        self._thd = TemporaryDirectory(dir=dir)
        self._genomes: Dict[str, GenomeData] = {}
        
    def insert(self, rec: GenomeData):
        self._genomes[rec.accession] = rec
        
    def __getitem__(self, *args, **kwargs):
        return self._genomes.__getitem__(*args, **kwargs)
        
    @property
    def path(self):
        return Path(self._thd.name)
      
    def genomes(self):
        return iter(self._genomes.values())
        
    def cleanup(self):
        self._thd.cleanup()


class MakiAtbManifest(Manifest, DatabasePackage):
    def __init__(self, columns: Iterable[str], concurrency: int, tmpdir: Path):
        super().__init__(columns)
        self.concurrency = concurrency
        self.tmp = tmpdir
        
    def select(self, accessions: Iterable[str]) -> MakiAtbManifest:
        submanifest = MakiAtbManifest(self.columns, self.concurrency, self.tmp)
        queries = set(accessions)
        
        for row in self:
            if row[ASSEMBLY_SCHEMA.sample] in queries:
                submanifest.insert(row)
        
        assert len(submanifest) == len(queries)
        return submanifest
  
    def retrieve_data(self, accessions: Iterable[str]) -> MakiAtbData:
        # retrieve query rows
        submanifest = self.select(accessions)
        
        # download data
        result = MakiAtbData(self.tmp)
        submanifest.download_batches(result.path)
        
        do_annot = self.has_annotations
        for row in self:
            result.insert(GenomeData(row[ASSEMBLY_SCHEMA.sample], row[ASSEMBLY_SCHEMA.species], None, str(result.path / "fa" / AssemblyItem.dest(row)), None if not do_annot else str((result.path / "gff" / BaktaItem.dest(row)).with_suffix(".gff"))))
        
        return result
  
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
