from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Optional, Tuple

import aiohttp
import asyncio
import csv
import hashlib
import tarfile
from pathlib import Path
from collections import defaultdict
from tqdm.asyncio import tqdm_asyncio

from .utils import open_maybe_gzip
from .lists import RemoteBatchFile, _FileLists


class Manifest:
    def __init__(self, columns: Iterable[str], annots: bool = False):
        self.columns = list(columns)
        self.rows: List[Dict[str, str]] = []
        self._annotations = annots
        
    def __iter__(self):
        return self.rows.__iter__()
        
    def insert(self, row: Dict[str, str]):
        assert all(c in row for c in self.columns), f"Row missing columns: {", ".join(c for c in self.columns if c not in row)}"
        self.rows.append(row)
        
    def has_annotations(self) -> bool:
        return self._annotations
        
    @classmethod
    def from_csv(cls, path: Path) -> Manifest:
        with open_maybe_gzip(path) as f:
            fieldnames = None
            has_annotations = False
            
            # Parse metadata
            for line in f:
                if line.startswith("#"):
                    if line.strip() == "# annotation filter: True":
                        has_annotations = True
                else:
                    fieldnames = line.strip().split(",")
                    break
          
            reader = csv.DictReader(f, fieldnames=fieldnames)
            
            assert reader.fieldnames
            mf = Manifest(reader.fieldnames, has_annotations)
            
            for r in reader:
                mf.insert(r)

            return mf
          
    def batch(self, dispatch: Callable[[Dict[str, str], ], Tuple[str, BatchItem]], lists: _FileLists):
        batches: Dict[RemoteBatchFile, List[BatchItem]] = defaultdict(list)
        unidentified_batches: List[str] = []
        
        for r in self.rows:
            batch_name, item = dispatch(r)
            try:
                batch_rec = lists[batch_name]
                batches[batch_rec].append(item)
            except KeyError:
                unidentified_batches.append(batch_name)
                
        if unidentified_batches:
            raise RuntimeError(f"Unidentified batch records: {", ".join(sorted(unidentified_batches))}")
          
        return batches


def md5sum(path: Path, chunk_size=8192):
    h = hashlib.md5()
    with open(path, "rb") as f:
        while chunk := f.read(chunk_size):
            h.update(chunk)
    return h.hexdigest()


@dataclass(frozen=False)
class BatchItem:
    filename: str
    md5: Optional[str] = None
    

class AssemblyItem(BatchItem):
    @classmethod
    def from_row(cls, r: Dict[str, str]) -> Tuple[str, AssemblyItem]:
        return r["osf_tarball_filename"], AssemblyItem(filename=f"{r["sample_accession"]}.fa")


class BaktaItem(BatchItem):
    @classmethod
    def from_row(cls, r: Dict[str, str]) -> Tuple[str, BaktaItem]:
        return r["tar_xz"], BaktaItem(filename=r["file_name"], md5=r["file_md5"])


def group_batches(rows, output_dir: Path, kind="assembly"):
    batches = defaultdict(list)

    for r in rows:
        url = r[f"{kind}_url"]

        item = {
            "path": r[f"{kind}_path"],
            "md5": r[f"{kind}_md5"],
            "sample": r["sample"],
        }

        target = output_dir / item["path"]

        # skip already valid files
        if target.exists():
            try:
                if md5sum(target) == item["md5"]:
                    continue
            except Exception:
                pass

        batches[url].append(item)

    return {k: v for k, v in batches.items() if v}


class AsyncBatchJob:
    def __init__(self, batch_rec: RemoteBatchFile, items: List[BatchItem], output_dir: Path, session):
        self.batch_rec = batch_rec
        self.items = items
        self.output_dir = output_dir / self.batch_rec.tar_name.removesuffix(".tar.xz")
        self.session = session

    async def run(self):
        extracted = 0

        # ignore already existing results
        queries: Dict[str, BatchItem] = {}
        
        for item in self.items:
            target = self.output_dir / item.filename
            
            if target.exists() and item.md5 and md5sum(target) == item.md5:
                continue
            
            queries[item.filename] = item

        if queries:
            async with self.session.get(self.batch_rec.download_url) as resp:
                resp.raise_for_status()

                # streaming tar read
                tar = tarfile.open(fileobj=resp.content, mode="r|*")

                for member in tar:
                    if member.name not in queries:
                        continue

                    item = queries.pop(member.name)
                    
                    target = self.output_dir / item.filename
                    target.parent.mkdir(parents=True, exist_ok=True)

                    extracted_file = tar.extractfile(member)
                    if extracted_file is None:
                        continue

                    with open(target, "wb") as f:
                        while chunk := extracted_file.read(8192):
                            f.write(chunk)

                    # validate
                    if item.md5 and md5sum(target) != item.md5:
                        raise RuntimeError(f"MD5 mismatch: {target}")

                    extracted += 1
            
            # check for unresolved queries
            if queries:
                raise RuntimeError(f"Could not find requested files in batch {self.batch_rec.tar_name}: {", ".join(queries.keys())}")

        return extracted


async def run_jobs_async(batches: Dict[RemoteBatchFile, List[BatchItem]], output_dir: Path, concurrency: int = 6):
    if not batches:
        print("Nothing to download")
        return 0

    connector = aiohttp.TCPConnector(limit=concurrency * 2)

    async with aiohttp.ClientSession(connector=connector) as session:
        sem = asyncio.Semaphore(concurrency)

        async def worker(batch_rec: RemoteBatchFile, items: List[BatchItem]):
            async with sem:
                job = AsyncBatchJob(batch_rec, items, output_dir, session)
                return await job.run()

        tasks = [worker(batch_rec, items) for batch_rec, items in batches.items()]

        results = []
        for coro in tqdm_asyncio.as_completed(tasks, total=len(tasks)):
            try:
                r = await coro
                results.append(r)
            except Exception as e:
                print(f"Batch failed: {e}")

        return sum(results)
