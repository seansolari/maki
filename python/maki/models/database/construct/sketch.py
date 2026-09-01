
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from pathlib import Path

from maki.models.database.manifest import DatabasePackage, GenomeData
from maki.sketch.core import (
    new_minhash,
    load_one_signature,
    save_one_signature,
    validate_signature,
    SketchParameters
)
from maki.sketch.store import SourmashSketchStore
from maki.utils.io import read_fasta, read_fasta_from_gff
from sourmash import SourmashSignature


def _worker_create_signature(
    data: tuple[GenomeData, SketchParameters, str | Path],
):
    """
    Worker process.
    """
    item, params, storage_root = data
    
    root = Path(storage_root)

    outfile = root / GenomeSketchStore._compute_storage_path(item.accession)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    
    # sketch file
    
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


class GenomeSketchStore(SourmashSketchStore):
    @staticmethod
    def signature_name(accession: str) -> str:
        return f"{accession}.sig.bz2"
    
    # --------------------------------------------------
    # sketching
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
        new_accessions = self.filter_existing(accessions)
        
        if new_accessions:
            with package.retrieve_data(new_accessions) as data:
                sketch_jobs = zip(data.genomes(), repeat(self.params), repeat(self.root))
                
                with ProcessPoolExecutor(max_workers=workers) as executor:
                    list(executor.map(_worker_create_signature, sketch_jobs))

        return accessions

    # --------------------------------------------------
    # loading
    # --------------------------------------------------

    def load_signature(
        self,
        accession: str,
    ) -> SourmashSignature:
        """Load a SourmashSignature object.
        """
        sketch_path = self.root / self._compute_storage_path(accession)
        
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
