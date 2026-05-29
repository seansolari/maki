from dataclasses import dataclass
import random
from typing import List, Tuple

from generate_templates import DNA, COMPLEMENT

def mutate(seq: str, n_mut: int, seed: int):
    random.seed(seed)
  
    template = list(seq)
    for _ in range(n_mut):
        i = random.randrange(len(seq))
        bases = [b for b in DNA if b != template[i]]
        template[i] = random.choice(bases)
    return "".join(template)

def simulate_reads(seq: str, coverage: float, read_len: int, insert_size: int) -> List[Tuple[str, str]]:
    cseq = "".join(COMPLEMENT[c] for c in seq)
  
    reads: List[Tuple[str, str]] = []
    n_reads = int(len(seq) * coverage / (2 * read_len))
    for _ in range(n_reads):
        start = random.randint(0, len(seq) - insert_size)
        r1 = seq[start:start+read_len]
        r2 = cseq[start+insert_size-read_len:start+insert_size]
        reads.append((r1, r2[::-1]))  # reverse second read
        
    return reads
  
@dataclass
class ReferenceToken:
    rid: int
    mutations: int
    coverage: float
    
def simulate_community(templates: List[str], members: List[ReferenceToken], seed: int, read_len: int = 100, insert_size: int = 250):
    all_reads = []
    
    for tkn in members:
        assert tkn.rid < len(templates)
        seq = mutate(templates[tkn.rid], tkn.mutations, seed)
        seed += 1
        all_reads.extend(simulate_reads(seq, tkn.coverage, read_len, insert_size))
    
    return all_reads
  