
from concurrent.futures import ProcessPoolExecutor
import logging
import os
from pathlib import Path
from typing import Iterable, Iterator

from maki.sketch.core import SketchParameters, singlesketch
from maki.sketch.store import SourmashSketchStore
from maki.utils.hashing import hash_strings

from .files import PairedReads, SampleData


logger = logging.getLogger(__name__)


def _singlesketch_worker(
    data: tuple[str | Path, Iterator[PairedReads], SketchParameters]
):
    outfile, file_iter, params = data
    
    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    # sketch data
    files = list(file_iter)
    
    singlesketch(
        (f for rp in files for f in (rp.r1.path, rp.r2.path)),
        outfile,
        params
    )
    
    # write hash
    tmp_hash = outfile.with_suffix(".tmpmd5")
    tmp_hash.write_text(MetagenomeSketchStore.md5(files))

    os.replace(tmp_hash, tmp_hash.with_suffix(".md5"))


class MetagenomeSketchStore(SourmashSketchStore):
    @staticmethod
    def signature_name(accession: str) -> str:
        return f"{accession}.sig"
    
    @staticmethod
    def md5(files: Iterable[PairedReads]):
        return hash_strings(str(f) for pr in files for f in (pr.r1.path, pr.r2.path))
    
    @classmethod
    def _compute_storage_path(cls, accession: str) -> str:
        return f"{accession}/{cls.signature_name(accession)}"
    
    # --------------------------------------------------
    # sketching
    # --------------------------------------------------
    
    def sketch_many(
        self,
        samples: list[SampleData],
        *,
        clear_old: bool = False,
        workers: int | None = None
    ) -> list[str]:
        """Sketch many metagenomes in parallel.
        """
        
        new_samples: set[str] = set()

        for sample in samples:
            if not self.valid_hash(sample.sample_name, sample.flatten(), clear_old=clear_old):
                new_samples.add(sample.sample_name)
            
        if new_samples:
            jobs = (
                (
                    self.root / self._compute_storage_path(sample.sample_name)[1],
                    sample.flatten(),
                    self.params
                )
                for sample in samples if sample.sample_name in new_samples
            )
            
            with ProcessPoolExecutor(max_workers=workers) as executor:
                list(executor.map(_singlesketch_worker, jobs))

        return list(new_samples)
    
    def valid_hash(
        self,
        sample_name: str,
        files: Iterator[PairedReads],
        *,
        clear_old: bool = False,
    ):
        _, relpath = self._compute_storage_path(sample_name)
        path = self.root / relpath
        
        if path.exists():
            existing_hash = path.with_suffix(".md5").read_text()
            this_hash = self.md5(files)
            
            if existing_hash == this_hash:
                logger.info("Skipping sample %s, identical hashes.", sample_name)
                return True
            
            elif not clear_old:
                logger.warning(
                    "Skipping sample %s but hash mismatch indicates "
                    "stale sketch. Call again with `--overwrite` to "
                    "overwrite old data.",
                    sample_name
                )
                return False
            
            else:
                logger.info("Removing data for sample %s due to stale hash.", sample_name)
                self.cleanup_accession(sample_name)
                return False
            
        else:
            return False
