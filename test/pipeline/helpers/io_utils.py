
"""I/O utilities for manifests and TSV parsing."""
from __future__ import annotations
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple
import gzip

def write_manifest(path: Path, fasta_paths: Iterable[Path]) -> None:
    """Write a plain-text manifest (one path per line)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for p in fasta_paths:
            fh.write(str(Path(p).resolve()) + "\n")

def read_tsv_gz(path: Path) -> List[List[str]]:
    """Read a gzipped TSV file into a list of rows (no header expected)."""
    rows: List[List[str]] = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            rows.append(line.split("\t"))
    return rows

# Parsers return typed tuples for convenience in assertions
class JaccardResult:
    def __init__(self, row: Sequence[str]) -> None:
        self.id1 = int(row[0])
        self.id2 = int(row[1])
        self.kmer_matches = int(row[2])
        self.length1 = int(row[3])
        self.length2 = int(row[4])
        self.jaccard = float(row[5])
        self.pval = float(row[6])
        self.r_lo = float(row[7])
        self.r_hi = float(row[8])
        self.lca = row[9]

def parse_jaccard_row(row: Sequence[str]) -> JaccardResult:
    """Parse a jaccard summary row (no header)."""
    return JaccardResult(row)

class RandistResult(JaccardResult):
    def __init__(self, row: Sequence[str]) -> None:
        super().__init__(row)
        self.cluster_1 = row[10]
        self.cluster_2 = row[11]

def parse_feature_row(row: Sequence[str]) -> RandistResult:
    """Parse a feature-table row (no header)."""
    return RandistResult(row)

def parse_prevalence_row(row: Sequence[str]) -> Tuple[float, int, int, str]:
    """Parse a k-mer prevalence row (no header)."""
    signature_prev = float(row[0])
    node_size = int(row[1])
    signature_size = int(row[2])
    lca = row[3]
    return (signature_prev, node_size, signature_size, lca)
