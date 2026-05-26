from collections import defaultdict
from typing import List, Tuple

def kmer_counts(fastq: List[Tuple[str, str]], k: int):
    counts = defaultdict(int)
    for rp in fastq:
        for seq in rp:
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i+k]
                counts[kmer] += 1
    return counts
