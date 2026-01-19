
"""GFF3 annotated genome generation helpers.

Build non-overlapping annotations across a genome and write GFF3 files that
include a trailing `##FASTA` block with the actual sequence. Also supports
per-segment mutation while preserving gene IDs.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import List, Mapping, Optional

@dataclass
class Annotation:
    """Represents a single annotated region on a genome."""
    start: int      # 1-based inclusive
    end: int        # 1-based inclusive
    strand: str     # '+' or '-'
    gene_id: str    # 'gene00001', used in Parent and ID attributes

def _poisson_knuth(lam: float, rng_random) -> int:
    """Stdlib Poisson sampler (Knuth). Suitable for moderate lam (~≤ 300)."""
    import math
    L = math.exp(-lam)
    k = 0
    p = 1.0
    while p > L:
        k += 1
        p *= rng_random()
    return k - 1

def build_nonoverlapping_annotations(L: int, mean_len: float, rng) -> List[Annotation]:
    """Cover genome [1..L] with sequential, non-overlapping annotations."""
    if L <= 0 or mean_len <= 0:
        return []
    anns: List[Annotation] = []
    pos = 1
    idx = 1
    # Strand balancing: alternate '+' and '-'
    strand_flag = True

    # Detect NumPy availability
    has_numpy = hasattr(rng, "poisson")
    import random
    r = rng if isinstance(rng, random.Random) else random.Random(42)  # deterministic fallback if needed

    while pos <= L:
        size = 1
        if has_numpy:
            size = max(1, int(rng.poisson(mean_len)))
        else:
            # Knuth Poisson fallback; for large mean_len, consider capping
            size = max(1, _poisson_knuth(mean_len, r.random))
        end = min(L, pos + size - 1)
        strand = "+" if strand_flag else "-"
        strand_flag = not strand_flag
        gene_id = f"gene{idx:05d}"
        anns.append(Annotation(start=pos, end=end, strand=strand, gene_id=gene_id))
        pos = end + 1
        idx += 1
    return anns

def write_gff3_with_fasta(
    path: Path,
    accession: str,
    seq: str,
    annotations: List[Annotation],
    extra_attrs_for_gene: Optional[Mapping[str, Mapping[str, str]]] = None,
) -> None:
    """Write GFF3 file with trailing FASTA block for `accession`."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("##gff-version 3\n")
        for ann in annotations:
            attrs = [f"ID={accession}-gene-{ann.gene_id}", f"Parent=gene-{ann.gene_id}"]
            if extra_attrs_for_gene and ann.gene_id in extra_attrs_for_gene:
                for k, v in extra_attrs_for_gene[ann.gene_id].items():
                    attrs.append(f"{k}={v}")
            attr_str = ";".join(attrs)
            # seqid, source, type, start, end, score, strand, phase, attributes
            fh.write(f"{accession}\tmaki-test\tgene\t{ann.start}\t{ann.end}\t.\t{ann.strand}\t.\t{attr_str}\n")
        fh.write("##FASTA\n")
        fh.write(f">{accession}\n")
        fh.write(seq + "\n")

def mutate_annotated_template_per_segment(
    template_seq: str,
    annotations: List[Annotation],
    segment_mu_rate_map: Mapping[str, float],
    rng,
) -> str:
    """Mutate `template_seq` per annotation using its specific mutation rate."""
    if not template_seq or not annotations:
        return template_seq
    # Avoid circular import at module import time
    from .genomes import mutate_genome
    seq_list = list(template_seq)
    for ann in annotations:
        mu = segment_mu_rate_map.get(ann.gene_id, 0.0)
        # Convert to 0-based slice
        start0 = ann.start - 1
        end0 = ann.end       # slicing end is exclusive
        segment = "".join(seq_list[start0:end0])
        mutated = mutate_genome(segment, mu, rng)
        seq_list[start0:end0] = list(mutated)
    return "".join(seq_list)
