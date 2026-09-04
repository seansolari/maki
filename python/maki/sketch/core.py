from __future__ import annotations
import bz2
import csv
from dataclasses import dataclass
import hashlib
import os
from pathlib import Path
from tempfile import TemporaryDirectory, TemporaryFile
from typing import Iterable, Iterator, List

from maki.utils.sp import run_in_current_env
from sourmash import (
    MinHash,
    SourmashSignature,
    load_signatures_from_json,
    save_signatures
)
from sourmash.sourmash_args import SaveSignaturesToLocation


FORMAT_VERSION = 1


# ---------------------------------------------------------------------------
# Shared data structures
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SketchParameters:
    ksize: int
    scaled: int
    seed: int = 42
    moltype: str = "DNA"
    track_abundance: bool = False

    def validate(self) -> None:
        if self.ksize <= 0:
            raise ValueError("ksize must be positive")
        if self.scaled <= 0:
            raise ValueError("scaled must be positive")
        if self.moltype != "DNA":
            raise ValueError("This implementation currently supports DNA only")
        if self.track_abundance:
            raise ValueError(
                "Reference sketches should not track abundance; "
                "cluster sketches represent sets of reference k-mers."
            )
            
    def to_str(self):
        return (
            f"k={self.ksize},scaled={self.scaled},"
            f"{'abund' if self.track_abundance else 'noabund'},"
            f"seed={self.seed}"
        )


def signature_params(signature: SourmashSignature) -> SketchParameters:
    mh = signature.minhash

    # sourmash exposes moltype on recent versions. These fallbacks make the
    # validation error clearer on earlier 4.x releases.
    if getattr(mh, "is_protein", False):
        moltype = "protein"
    elif getattr(mh, "dayhoff", False):
        moltype = "dayhoff"
    elif getattr(mh, "hp", False):
        moltype = "hp"
    else:
        moltype = "DNA"

    return SketchParameters(
        ksize=mh.ksize,
        scaled=mh.scaled,
        seed=mh.seed,
        moltype=moltype,
        track_abundance=mh.track_abundance,
    )


@dataclass(frozen=True)
class PairwiseResults:
    fieldnames: list[str]
    data: list[dict[str, str]]


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------

def safe_filename(value: str) -> str:
    digest = hashlib.sha256(value.encode("utf-8")).hexdigest()[:16]
    cleaned = "".join(
        c if c.isalnum() or c in "._-" else "_"
        for c in value
    ).strip("._")

    if not cleaned:
        cleaned = "genome"

    return f"{cleaned[:80]}-{digest}"
  

def new_minhash(params: SketchParameters) -> MinHash:
    return MinHash(
        n=0,
        ksize=params.ksize,
        scaled=params.scaled,
        seed=params.seed,
        track_abundance=params.track_abundance,
    )


def save_one_signature(
    signature: SourmashSignature,
    output_path: str | Path
) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.suffix != ".bz2":
        raise RuntimeError(f"Filename must have .bz2 suffix: {output_path}")

    # Write atomically so an interrupted worker does not leave a valid-looking
    # partial signature.
    tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")

    with bz2.open(tmp_path, "wt") as fp:
        save_signatures([signature], fp=fp)

    os.replace(tmp_path, output_path)


def load_one_signature(path: str | Path) -> SourmashSignature:
    with bz2.open(path, "rt") as fp:
        data = fp.read()
    
    signatures = list(load_signatures_from_json(data))

    if len(signatures) != 1:
        raise ValueError(
            f"Expected exactly one signature in {path}; "
            f"found {len(signatures)}"
        )

    return signatures[0]


def validate_signature(
    signature: SourmashSignature,
    expected: SketchParameters,
    *,
    allow_finer_scaled: bool = True,
) -> None:
    observed = signature_params(signature)

    if observed.ksize != expected.ksize:
        raise ValueError(
            f"Incompatible k-mer size: query={observed.ksize}, "
            f"database={expected.ksize}"
        )

    if observed.seed != expected.seed:
        raise ValueError(
            f"Incompatible hash seed: query={observed.seed}, "
            f"database={expected.seed}"
        )

    if observed.moltype != expected.moltype:
        raise ValueError(
            f"Incompatible molecule type: query={observed.moltype}, "
            f"database={expected.moltype}"
        )

    if observed.scaled == 0:
        raise ValueError("The query must use a scaled/FracMinHash sketch")

    if allow_finer_scaled:
        # Lower scaled means a denser sketch. It can be downsampled to the
        # database scale. A coarser query cannot recover omitted hashes.
        if observed.scaled > expected.scaled:
            raise ValueError(
                f"Query sketch is too coarse: query scaled={observed.scaled}, "
                f"database scaled={expected.scaled}. Build the query at "
                f"scaled <= {expected.scaled}."
            )
    elif observed.scaled != expected.scaled:
        raise ValueError(
            f"Incompatible scaled value: query={observed.scaled}, "
            f"database={expected.scaled}"
        )


def downsample_signature(
    signature: SourmashSignature,
    scaled: int,
) -> SourmashSignature:
    mh = signature.minhash

    if mh.scaled == scaled:
        return signature

    if mh.scaled > scaled:
        raise ValueError(
            f"Cannot upsample scaled={mh.scaled} to the denser scaled={scaled}"
        )

    downsampled = mh.downsample(scaled=scaled)

    return SourmashSignature(
        downsampled,
        name=signature.name,
        filename=signature.filename,
    )


def singlesketch(files: Iterator[Path], outfile: Path, params: SketchParameters):
    tmp_out = outfile.with_suffix(".sig.tmp")
    
    run_in_current_env(
        [
            "sourmash", "scripts", "singlesketch",
            "-o", tmp_out,
            "-p", params.to_str(),
            "-I", params.moltype,
            *files
        ]
    )
    
    os.replace(tmp_out, outfile)


def zip_signatures(signatures: Iterable[SourmashSignature], zip_file: Path):
    assert zip_file.suffix == ".zip"
    
    tmp_zip = zip_file.with_suffix(".zip.tmp")
    
    try:
        with SaveSignaturesToLocation(str(tmp_zip)) as zar:
            for sig in signatures:
                zar.add(sig)
    
    except Exception:
        tmp_zip.unlink(missing_ok=True)
        
        raise

    else:
        os.replace(tmp_zip, zip_file)


def build_standalone_manifest(
    signatures: Iterable[str],
    outfile: str | Path
):
    outfile = Path(outfile).resolve()
    outfile.mkdir(parents=True, exist_ok=True)
    
    tmp_outfile = outfile.with_suffix(".tmp")
    
    with TemporaryDirectory(dir=outfile.parent) as tmp:
        # Create filelist
        filelist = Path(tmp) / "filelist.txt"

        with filelist.open("wt") as f:
            for sigpath in signatures:
                f.write(f"{sigpath}\n")
        
        # Create standalone manifest
        # sourmash sig collect pathlist.txt -o summary-manifest.csv -F csv
        run_in_current_env(
            [
                "sourmash", "sig", "collect", filelist,
                "-o", tmp_outfile,
                "-F", "csv"
            ]
        )

    tmp_outfile.rename(outfile)

        
def build_rocksdb_index(
    signatures: List[SourmashSignature],
    params: SketchParameters,
    index_dir: str | Path,
    threads: int = 8,
):
    """
    Build a RocksDB sourmash index from a zip containing signatures.

    Parameters
    ----------
    sig_zip : str
        Path to zip containing *.sig files.
    index_dir : str
        Output RocksDB database directory.
    threads : int
        Parallel threads.
    """
    
    index_dir = Path(index_dir)
    assert index_dir.suffix == ".rocksdb"
    
    tmp_zip = index_dir.with_suffix(".zip")
    zip_signatures(signatures, tmp_zip)
    
    # Build database
    try:
        run_in_current_env(
            [
                "sourmash", "scripts", "index",
                "-F", "rocksdb",
                str(index_dir),
                tmp_zip,
                "-k", str(params.ksize),
                "-s", str(params.scaled),
                "--cores", str(threads),
            ]
        )
    
    finally:
        tmp_zip.unlink()

    return index_dir


def pairwise_ani_comparison(
    signatures: Iterable[SourmashSignature],
    params: SketchParameters,
    tmp_prefix: str | Path,
    processes: int = 8
):
    with TemporaryDirectory(dir=tmp_prefix) as tmp:
        root = Path(tmp)
        
        zip_file = root / "signatures.zip"
        zip_signatures(signatures, zip_file)
        
        # pairwise comparisons
        stat_file = root / "pairwise.csv"

        run_in_current_env(
            [
                "sourmash", "scripts", "pairwise",
                zip_file,
                "-o", stat_file,
                "-k", str(params.ksize),
                "-s", str(params.scaled),
                "--cores", str(processes),
                "-a"
            ]
        )
        
        with stat_file.open("rt", encoding='utf-8') as fh:
            reader = csv.DictReader(fh)
            assert reader.fieldnames
            return PairwiseResults(list(reader.fieldnames), list(reader))
        

def cluster_from_pairwise(
    csv_file: str | Path,
    ani: float,
    cores: int = 8
) -> list[list[str]]:
    csv_file = Path(csv_file)
    
    assert csv_file.name.endswith(".csv")

    result_file = csv_file.parent / f"{csv_file.name[:-4]}-clusters_{int(round(ani * 1.0e6, 6))}.csv"
    
    run_in_current_env(
        [
            "sourmash", "scripts", "cluster",
            "-o", result_file,
            "--similarity-column", "average_containment_ani",
            "-t", str(ani),
            "-c", str(cores),
            csv_file,
        ]
    )
    
    with result_file.open("rt", encoding='utf-8') as fh:
        reader = csv.DictReader(fh)
        return [
            r["nodes"].split(";")
            for r in reader
        ]
