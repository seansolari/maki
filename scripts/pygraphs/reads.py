from __future__ import annotations
from dataclasses import dataclass
from typing import List
import numpy as np

from genome import ALPHA_SET, Genome, basic_species


@dataclass(frozen=True)
class SequencingError:
    pos: int
    new: str


@dataclass
class ReadPair:
    tag: str
    fwd_seq: str
    rev_seq: str
    fwd_errs: List[SequencingError]
    rev_errs: List[SequencingError]
    
    @staticmethod
    def _put_errs(_seq: str, _errs: List[SequencingError]):
        data = list(_seq)
        for err in _errs:
            data[err.pos] = err.new
        return "".join(data)
    
    @property
    def fwd(self) -> str:
        return self._put_errs(self.fwd_seq, self.fwd_errs)
    
    @property
    def rev(self) -> str:
        return self._put_errs(self.rev_seq, self.rev_errs)
    
    
def sample_errors(seq: str, nerrs: int) -> List[SequencingError]:
    if not nerrs:
        return []
    else:
        rng = np.random.default_rng(seed=abs(hash(seq)))
        return [
            SequencingError(i, rng.choice([c for c in ALPHA_SET if c != seq[i]]))
            for i in rng.integers(0, len(seq), nerrs)
        ]
    
    
def sample_reads(genome: Genome, n: int, r: int, insert: int, e: float, seed: int = 1) -> List[ReadPair]:
    rng = np.random.default_rng(seed=seed)
    
    fwd_starts = rng.integers(0, len(genome) - ((2 * r) + insert), n, dtype=np.uint64)
    fwd_starts.sort()
    rev_starts = fwd_starts + (r + insert)
    
    breaks = genome.allele_ends()
    fwd_starts = np.split(fwd_starts, np.searchsorted(fwd_starts, breaks[:-1]))
    rev_starts = np.split(rev_starts, np.searchsorted(rev_starts, breaks[:-1]))
    
    for i in range(1, len(fwd_starts)):
        fwd_starts[i] -= breaks[i-1]
    for i in range(1, len(rev_starts)):
        rev_starts[i] -= breaks[i-1]
    
    reads = []
    fi, fj = 0, 0
    ri, rj = 0, 0
    while len(reads) < n:
        while fj >= fwd_starts[fi].shape[0]:
            fi += 1
            fj = 0
        assert fi < len(fwd_starts)
        fwd_seq = genome.seq_from(fi, fj, r)
        fwd_errs = sample_errors(fwd_seq, round(e * len(fwd_seq)))
        fj += 1
        
        while rj >= rev_starts[ri].shape[0]:
            ri += 1
            rj = 0
        assert ri < len(rev_starts)
        rev_seq = genome.seq_from(ri, rj, r)
        rev_errs = sample_errors(rev_seq, round(e * len(rev_seq)))
        rj += 1
        
        reads.append(ReadPair(f"read{len(reads)+1}", fwd_seq, rev_seq, fwd_errs, rev_errs))
    
    return reads


def test_reads():
    sp = basic_species(5, 1000, 0.01)
    g1 = sp.generate()
    
    for read in sample_reads(g1, 10, 150, 500, 0.01):
        print(read.tag, f"{read.fwd[:10]}...", f"{read.rev[:10]}...", sep = "\t")
    
    
if __name__ == "__main__":
    test_reads()
