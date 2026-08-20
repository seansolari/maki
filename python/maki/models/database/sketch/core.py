from __future__ import annotations
from dataclasses import dataclass
import hashlib
import os
from pathlib import Path

from sourmash import (
    MinHash,
    SourmashSignature,
    load_file_as_signatures,
    save_signatures
)


FORMAT_VERSION = 1


# ---------------------------------------------------------------------------
# Shared data structures
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SketchParameters:
    ksize: int = 31
    scaled: int = 2000
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
            
            
@dataclass(frozen=True)
class GenomeSketchRecord:
    genome_id: str
    sketch_path: str
    cluster_id: int


@dataclass(frozen=True)
class ContainmentHit:
    cluster_id: int
    containment: float
    intersect_hashes: int | None = None
    query_hashes: int | None = None
    cluster_hashes: int | None = None


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
    output_path: str | Path,
) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Write atomically so an interrupted worker does not leave a valid-looking
    # partial signature.
    tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")

    with tmp_path.open("wt") as fp:
        save_signatures([signature], fp=fp)

    os.replace(tmp_path, output_path)


def load_one_signature(path: str | Path) -> SourmashSignature:
    signatures = list(load_file_as_signatures(str(path)))

    if len(signatures) != 1:
        raise ValueError(
            f"Expected exactly one signature in {path}; "
            f"found {len(signatures)}"
        )

    return signatures[0]


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

