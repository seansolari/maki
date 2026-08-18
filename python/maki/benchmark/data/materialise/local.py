
from dataclasses import dataclass
from pathlib import Path
import random
import tempfile
from typing import Optional

from maki.benchmark.data.local.dataset import LocalBenchmarkDataset
from maki.benchmark.data.materialise.materialiser import DatasetMaterialiser
from maki.benchmark.data.models import DatasetFile, DatasetReference
from maki.classify.utils import fastq_sample_name
from maki.utils.io import open_maybe_gzip


def read_fastq_record(handle):
    """
    Read one FASTQ record as a 4-line tuple.
    Returns None at EOF.
    """
    header = handle.readline()
    if not header:
        return None

    seq = handle.readline()
    plus = handle.readline()
    qual = handle.readline()

    if not seq or not plus or not qual:
        raise ValueError("Truncated FASTQ record detected")

    return header, seq, plus, qual


def subsample_paired_fastq(
    r1_path: Path,
    r2_path: Path,
    n: int,
    seed: int,
    out_r1: Path,
    out_r2: Path,
    check_names=True,
):
    """
    Randomly select exactly n read pairs from paired FASTQ files using reservoir sampling.

    Parameters
    ----------
    r1_path : str or Path
        Forward reads FASTQ file. Can be .fastq, .fq, .fastq.gz, or .fq.gz.

    r2_path : str or Path
        Reverse reads FASTQ file. Can be compressed or uncompressed.

    n : int
        Number of read pairs to sample.

    seed : int
        Random seed for reproducibility.

    output_prefix : str
        Prefix for output files.

    check_names : bool
        If True, checks that paired read names appear to match.

    Returns
    -------
    tuple
        Paths to written R1 and R2 output files.
    """

    if n <= 0:
        raise ValueError("n must be a positive integer")

    rng = random.Random(seed)

    reservoir = []
    total_pairs = 0

    with open_maybe_gzip(r1_path, "rt") as f1, open_maybe_gzip(r2_path, "rt") as f2:
        while True:
            rec1 = read_fastq_record(f1)
            rec2 = read_fastq_record(f2)

            if rec1 is None and rec2 is None:
                break

            if rec1 is None or rec2 is None:
                raise ValueError("R1 and R2 files contain different numbers of records")

            if check_names:
                name1 = rec1[0].split()[0]
                name2 = rec2[0].split()[0]

                # Remove common paired-end suffixes
                clean1 = name1.replace("/1", "").replace(" 1:", "")
                clean2 = name2.replace("/2", "").replace(" 2:", "")

                if clean1 != clean2:
                    raise ValueError(
                        f"Read name mismatch at pair {total_pairs + 1}: "
                        f"{name1} vs {name2}"
                    )

            total_pairs += 1
            pair = (rec1, rec2)

            if len(reservoir) < n:
                reservoir.append(pair)
            else:
                j = rng.randint(0, total_pairs - 1)
                if j < n:
                    reservoir[j] = pair

    if total_pairs < n:
        raise ValueError(
            f"Requested {n} pairs, but only found {total_pairs} pairs"
        )

    with open_maybe_gzip(out_r1, "wt") as o1, open_maybe_gzip(out_r2, "wt") as o2:
        for rec1, rec2 in reservoir:
            o1.write("".join(rec1))
            o2.write("".join(rec2))

    return out_r1, out_r2


@dataclass
class PairedEndDatasetToken(DatasetReference):
    forward: Path
    reverse: Path
    suffix: str
    
    subsample: Optional[int] = None
    seed: int = 1
    
    @classmethod
    def from_local(
        cls,
        domain: str,
        forward: Path,
        reverse: Path,
        subsample: Optional[int] = None,
        seed: int = 1
    ):
        name_data = fastq_sample_name(forward)
        
        return cls(
            id=name_data.sample_name,
            name=name_data.sample_name,
            description=None,
            domain=domain,
            data_type="paired-end",
            forward=forward,
            reverse=reverse,
            suffix=name_data.suffix,
            subsample=subsample,
            seed=seed
        )
    
    
class LocalDatasetMaterialiser(DatasetMaterialiser):
    input_type = PairedEndDatasetToken
    
    def __init__(self) -> None:
        super().__init__()
    
    def materialise(self, token: PairedEndDatasetToken) -> LocalBenchmarkDataset:
        if token.subsample:
            tmp = Path(tempfile.mkdtemp(prefix="subsample_", dir=token.forward.parent))
            
            fwd = tmp / f"{token.id}_{token.subsample}.R1{token.suffix}"
            rev = tmp / f"{token.id}_{token.subsample}.R2{token.suffix}"
            
            subsample_paired_fastq(token.forward, token.reverse, token.subsample, token.seed, fwd, rev, check_names=False)
        else:
            tmp, fwd, rev = None, token.forward, token.reverse
            
        return LocalBenchmarkDataset(
            id=token.id,
            name=token.name,
            description=token.description,
            domain=token.domain,
            data_type=token.data_type,
            organism=None,
            assay=None,
            source=None,
            license=None,
            files=[
                DatasetFile("forward", str(fwd), fwd.name.removesuffix(token.suffix)),
                DatasetFile("reverse", str(rev), rev.name.removesuffix(token.suffix))
            ],
            temporary_directory=tmp
        )
        