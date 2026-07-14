from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass
import math
from typing import List, Optional

import numpy as np


ALPHA_SET = {"A", "C", "G", "T"}


@dataclass(frozen=True)
class Allele:
    locus_tag: str
    allele_id: int
    sequence: str
    
    @property
    def id(self) -> str:
        return f"{self.locus_tag}-{self.allele_id}"
    

class AlleleGenerator(ABC):
    @abstractmethod
    def generate(self) -> Allele: ...
    

class MutationModel(ABC):
    @abstractmethod
    def mutate(self, template: str) -> str: ...


@dataclass(frozen=True)
class SimpleMutationModel(MutationModel):
    r: float
    
    def mutate(self, template: str) -> str:
        result = list(template)
        for i in np.random.randint(0, len(template), size=math.floor(self.r * len(template))):
            result[i] = np.random.choice([c for c in ALPHA_SET if c != result[i]])
        return "".join(result)
    

@dataclass(frozen=True)
class ConservedDomainMutationModel(SimpleMutationModel):
    domain_start: int
    domain_length: int
    
    def mutate(self, template: str) -> str:
        result = list(template)
        L = len(template) - self.domain_length
        for i in np.random.randint(0, range(L), size=min(math.floor(self.r * len(template)), L)):
            if i >= self.domain_start:
                i += self.domain_length
            result[i] = np.random.choice([c for c in ALPHA_SET if c != result[i]])
        return "".join(result)
    

class ProteinFamily(AlleleGenerator):
    def __init__(self, name: str, template: str, mmodel: MutationModel) -> None:
        self.name = name
        self.template = template
        self.model = mmodel
        self.i = 0
        
    def generate(self) -> Allele:
        self.i += 1
        return Allele(self.name, self.i, self.model.mutate(self.template))
    

class RepeatModel(ABC):
    @abstractmethod
    def num_repeats(self) -> int: ...


@dataclass(frozen=True)
class ConservedRepeat(RepeatModel):
    n: int
    
    def num_repeats(self) -> int:
        return self.n
    
    
@dataclass(frozen=True)
class VariableRepeat(RepeatModel):
    lam: int
    
    def num_repeats(self) -> int:
        return 1 + np.random.poisson(self.lam)
    

class RepeatFamily(AlleleGenerator):
    def __init__(self, name: str, motif: str, mmodel: RepeatModel) -> None:
        self.name = name
        self.motif = motif
        self.model = mmodel
        self.i = 0
    
    def generate(self) -> Allele:
        self.i += 1
        return Allele(self.name, self.i, self.motif * self.model.num_repeats())
    
    
class StructureGenerator(ABC):
    @abstractmethod
    def modify(self, ids: List[int]) -> List[int]: ...
    
    
class RearrangeFirstGene(StructureGenerator):
    def modify(self, ids: List[int]):
        if len(ids) <= 1:
            return ids
        
        src = 1 if len(ids) > 2 else 0
        tmp = ids[src]
        ids[src] = ids[src + 1]
        ids[src + 1] = tmp
        
        return ids

    
class DeleteFirstGene(StructureGenerator):
    def modify(self, ids: List[int]):
        return [x for i, x in enumerate(ids) if i != 1]
    

@dataclass(frozen=True)
class Genome:
    name: str
    elements: List[Allele]
    
    def __repr__(self):
        return f"Genome {self.name}, {len(self.elements)} elements: {', '.join(a.id for a in self.elements)}"
    
    def __len__(self):
        return sum(len(e.sequence) for e in self.elements)
    
    def allele_ends(self) -> np.ndarray:
        pos = 0
        result = []
        for elem in self.elements:
            pos += len(elem.sequence)
            result.append(pos)
        return np.array(result, dtype=np.uint64)
    
    def seq_from(self, elem: int, pos: int, length: int) -> str:
        chunks: List[str] = []
        while length:
            chunks.append(self.elements[elem].sequence[pos:pos+length])
            elem += 1
            pos = 0
            length -= len(chunks[-1])
        return "".join(chunks)
    
    
class Species:
    def __init__(self, name: str, alleles: List[AlleleGenerator], modifiers: Optional[List[StructureGenerator]] = None) -> None:
        self.name = name
        self.alleles = alleles
        self.modifiers = modifiers or []
        self.i = 0
        
    def generate(self) -> Genome:
        inds = list(range(len(self.alleles)))

        for mx in self.modifiers:
            inds = mx.modify(inds)
            
        self.i += 1
        return Genome(f"{self.name}-{self.i}", [self.alleles[i].generate() for i in inds])


### Tests


def basic_species(n: int, length: int, r: float):
    rng = np.random.default_rng(seed=1)
    return Species("Species1", [ProteinFamily(f"family{i+1}", "".join(rng.choice(list(ALPHA_SET), length)), SimpleMutationModel(r)) for i in range(n)])


def test_alleles():
    sp = basic_species(2, 20, 0.1)
    g1 = sp.generate()
    g2 = sp.generate()
    
    print(g1)
    for elem in g1.elements:
        print(elem.id, f"{elem.sequence[:10]}...", sep="\t")
    
    print(g2)
    for elem in g2.elements:
        print(elem.id, f"{elem.sequence[:10]}...", sep="\t")


if __name__ == "__main__":
    test_alleles()
