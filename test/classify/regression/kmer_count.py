
from typing import Dict, List, Tuple
from generate_templates import COMPLEMENT

def kmer_counts(fastq: List[Tuple[str, str]], k: int) -> Dict[str, int]:
    counts = {}
    for rp in fastq:
        for seq in rp:
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i+k]
                try:
                    counts[kmer] += 1
                except KeyError:
                    counts[kmer] = 1
    return counts
    
def match_kmers(kmers: Dict[str, int], templates: List[str]):
    k = len(next(iter(kmers)))
    indexes = [kmer_counts([(seq, "".join(COMPLEMENT[c] for c in seq[::-1]))], k) for seq in templates]
    
    results = {}
    for kmer in kmers:
        idxs = [i for i, index in enumerate(indexes) if kmer in index]
        if idxs:
          results[kmer] = idxs
    return results
