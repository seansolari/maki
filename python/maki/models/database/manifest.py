
from __future__ import annotations
from abc import ABC, abstractmethod
import csv
from dataclasses import dataclass, astuple, fields
from pathlib import Path
from typing import Dict, Iterable, Iterator, Optional, Tuple, overload

from maki.utils.io import open_maybe_gzip, write_maybe_gzip


@dataclass
class GenomeRecord:
    accession: str # unique sequence ID
    taxonomy: str # user-supplied taxonomic group
    taxid: Optional[str] # validated taxonomy ID calculated from `taxonomy`
    
    @classmethod
    def from_dict(cls, r: Dict[str, str]):
        return cls(accession=r["accession"], taxonomy=r["taxonomy"], taxid=r["taxid"] or None)
      
    def astuple(self) -> Tuple[str, ...]:
        return (self.accession, self.taxonomy, self.taxid or '')
      
    def __eq__(self, rhs) -> bool:
        if not isinstance(rhs, GenomeRecord):
            return False
        return (self.accession == rhs.accession) and (self.taxonomy == rhs.taxonomy)

    def __hash__(self) -> int:
        return (self.accession, self.taxonomy).__hash__()


@dataclass
class GenomeData(GenomeRecord):
    fasta: str
    gff: Optional[str]
    
    @classmethod
    def from_dict(cls, r: Dict[str, str]):
        return cls(accession=r["accession"], taxonomy=r["taxonomy"], taxid=r["taxid"] or None, fasta=r["fasta"], gff=r["gff"] or None)
      
    def astuple(self) -> Tuple[str, ...]:
        return (self.accession, self.taxonomy, self.taxid or '', self.fasta, self.gff or '')


class SequencePackage(ABC):
    @abstractmethod
    def genomes(self) -> Iterator[GenomeData]:
        pass
      
    def __enter__(self):
        return self
    
    @abstractmethod
    def __exit__(self, exc_type, exc_val, exc_tb):
        pass
      

class DatabasePackage(ABC):
    @abstractmethod
    def records(self) -> Iterator[GenomeRecord]:
        pass
  
    @abstractmethod
    def retrieve_data(self, accessions: Iterable[str]) -> SequencePackage:
        pass


@dataclass(frozen=True)
class ManifestSchema:
    accession: str
    taxonomy: str
    taxid: Optional[str] = None
    
    def parse_row(self, r: Dict[str, str]) -> GenomeRecord:
        return GenomeRecord(accession=r[self.accession], taxonomy=r[self.taxonomy], taxid=None if not self.taxid else r[self.taxid])
    

@dataclass(frozen=True)
class GenomeSchema:
    accession: str
    taxonomy: str
    fasta: str
    gff: Optional[str]

    def parse_row(self, r: Dict[str, str]) -> GenomeData:
        return GenomeData(accession=r[self.accession], taxonomy=r[self.taxonomy], taxid=None, fasta=r[self.taxonomy], gff=None if not self.gff else r[self.gff])


class Manifest[T: (GenomeData, GenomeRecord)](DatabasePackage, SequencePackage):
    def __init__(self, records: Optional[Iterable[T]] = None) -> None:
        self._records: Dict[str, T] = {}
        
        if records:
            for r in records:
                self.insert(r)

    def records(self):
        return self._records.values().__iter__()
      
    def retrieve_data(self, accessions: Iterable[str]):
        return Manifest(self._records[acc] for acc in accessions)
      
    def genomes(self):
        for v in self.records():
            if isinstance(v, GenomeData):
                yield v

    def __exit__(self, exc_type, exc_val, exc_tb):
        return
      
    def __len__(self):
        return self._records.__len__()

    def insert(self, rec: T):
        self._records[rec.accession] = rec
    
    def save(self, path: Path):
        it = self._records.values().__iter__()
        try:
            first_record = next(it)
        except StopIteration:
            print(f"[warning] no records, not saving to {path}")
            return
      
        with write_maybe_gzip(path) as f:
            fh = csv.writer(f)
            
            fh.writerow((f.name for f in fields(first_record)))
            fh.writerow(first_record.astuple())
            fh.writerows((r.astuple() for r in it))
            

@overload
def read_manifest(path: Path, schema: ManifestSchema) -> Manifest[GenomeRecord]: ...

@overload
def read_manifest(path: Path, schema: GenomeSchema) -> Manifest[GenomeData]: ...

def read_manifest[T: (ManifestSchema, GenomeSchema)](path: Path, schema: T) -> Manifest[GenomeRecord] | Manifest[GenomeData]:
    T = GenomeData if isinstance(schema, GenomeSchema) else GenomeRecord
    result = Manifest[T]()
    
    with open_maybe_gzip(path) as f:
        reader = csv.DictReader(f)
        
        if not reader.fieldnames:
            raise RuntimeError(f"Corrupted manifest: {path}")
        elif not all(c in reader.fieldnames for c in astuple(schema)):
            raise RuntimeError(f"Expected columns: {", ".join(astuple(schema))}")
          
        for row in reader:
            result.insert(schema.parse_row(row))

    return result
    