from __future__ import annotations
import hashlib
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from pathlib import Path

from maki.utils.io import read_fasta, read_fasta_from_gff
from maki.models.database.manifest import DatabasePackage, GenomeData
from sourmash import SourmashSignature

from .core import SketchParameters, load_one_signature, new_minhash, save_one_signature, validate_signature


def _compute_storage_path(accession: str) -> tuple[str, str]:
    """Create a deterministic sharded path with .bz2 extension.
    """
    h = hashlib.sha256(accession.encode()).hexdigest()

    shard1 = h[:2]
    shard2 = h[2:4]

    relpath = f"{shard1}/{shard2}/{accession}.sig.bz2"

    return h, relpath


def _worker_create_signature(
    data: tuple[GenomeData, SketchParameters, str | Path],
) -> tuple[str, str]:
    """
    Worker process.
    """
    item, params, storage_root = data
    
    root = Path(storage_root)

    _, relpath = _compute_storage_path(item.accession)

    outfile = root / relpath
    outfile.parent.mkdir(parents=True, exist_ok=True)
    
    # setch file
    
    mh = new_minhash(params)
    if item.fasta:
        fname = item.fasta
        it = read_fasta(item.fasta)
    elif item.gff:
        fname = item.gff
        it = read_fasta_from_gff(item.gff)
    else:
        raise RuntimeError(f"Genome record {item.accession} did not supply FASTA or GFF data.")
    
    for _, sequence in it:
        # force=True skips k-mers containing ambiguous characters rather than
        # failing the complete genome.
        mh.add_sequence(sequence, force=True)

    if len(mh) == 0:
        raise ValueError(
            f"Genome {item.accession!r} produced an empty sketch"
        )

    signature = SourmashSignature(
        mh,
        name=item.accession,
        filename=str(fname),
    )

    save_one_signature(signature, outfile)
    
    return item.accession, str(relpath)


class SourmashSketchStore:
    """
    Large-scale sourmash signature manager.

    Features
    --------
    * SHA256 sharded filesystem layout
    * gzip-compressed signatures
    * SQLite index
    * multiprocessing batch sketch generation
    * random signature loading
    """

    def __init__(
        self,
        root_dir: str | Path,
        *,
        params: SketchParameters
    ):
        self.root = Path(root_dir)
        self.params = params
        self.params.validate()

    # --------------------------------------------------
    # batch sketching
    # --------------------------------------------------

    def sketch_many(
        self,
        package: DatabasePackage,
        accessions: tuple[str, ...],
        *,
        workers: int | None = None
    ):
        """Sketch many FASTA files in parallel.
        """
        new_accessions = {
            accession
            for accession in accessions
            if not self._sketch_exists(accession)
        }
        
        if new_accessions:
            with package.retrieve_data(new_accessions) as data:
                sketch_jobs = zip(data.genomes(), repeat(self.params), repeat(self.root))
                
                with ProcessPoolExecutor(max_workers=workers) as executor:
                    list(executor.map(_worker_create_signature, sketch_jobs))

        return accessions
    
    def _sketch_exists(self, accession: str):
        _, relpath = _compute_storage_path(accession)
        return (self.root / relpath).exists()

    # --------------------------------------------------
    # loading
    # --------------------------------------------------

    def load_signature(
        self,
        accession: str,
    ) -> SourmashSignature:
        """
        Load a SourmashSignature object.
        """
        _, relpath = _compute_storage_path(accession)
        sketch_path = self.root / relpath
        
        if not sketch_path.exists():
            raise ValueError(
                f"No sketch found for {accession} within {self.root}."
            )
        
        sig = load_one_signature(sketch_path)
        validate_signature(sig, self.params, allow_finer_scaled=False)
        
        return sig

    def load_many(
        self,
        accessions: list[str],
    ) -> dict[str, SourmashSignature]:
        return {
            acc: self.load_signature(acc)
            for acc in accessions
        }
