
from abc import ABC, abstractmethod
import csv
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from maki.utils.io import open_maybe_gzip


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


@dataclass
class GenomeData(GenomeRecord):
    fasta: str
    gff: Optional[str]


class SequencePackage(ABC):
    @abstractmethod
    def records(self) -> Iterable[GenomeData]:
        pass
      
    @abstractmethod
    def cleanup(self) -> None:
        pass
      
    def __enter__(self):
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cleanup()
        return False


class DatabasePackage(ABC):
    @abstractmethod
    def records(self) -> Iterable[GenomeRecord]:
        pass
  
    @abstractmethod
    def retrieve_data(self, accessions: Iterable[str]) -> SequencePackage:
        pass


class Manifest:
    REQUIRED_COLUMNS = {"accession"}

    def __init__(self, records):
        self.records: List[GenomeRecord] = records

    @classmethod
    def import_csv(cls, path: Path):
        """Create a new manifest by inspecting a csv.
        """
        with open_maybe_gzip(path) as f:
            reader = csv.DictReader(f)
            cols = reader.fieldnames
            
            if not cols:
                raise RuntimeError(f"Could not detect column names in manifest: {path}")

            # auto-detect taxonomy column
            tax_column = None
            for c in cols:
                if c.lower() in ["taxid", "tax_id", "taxonomy", "organism", "name"]:
                    tax_column = c
                    break

            if not tax_column:
                raise ValueError(
                    "Could not detect taxonomy column.\n"
                    "Expected one of: taxid, tax_id, taxonomy, organism, name"
                )

            records = []
            for row in reader:
                records.append(GenomeRecord(accession=row["accession"], taxonomy=row[tax_column], taxid=None))

        print(f"Detected taxonomy column: '{tax_column}'")

        return cls(records)
      
    @classmethod
    def load(cls, path: Path):
        with open_maybe_gzip(path) as f:
            reader = csv.DictReader(f)
            
            if not reader.fieldnames:
                raise RuntimeError(f"Corrupted manifest: {path}")
              
            records = []
            for row in reader:
                records.append(GenomeRecord.from_dict(row))
                
            return cls(records)
          
    def save(self, path: Path):
        with open_maybe_gzip(path) as f:
            fh = csv.writer(f)
            fh.writerow((f.name for f in fields(GenomeRecord)))
            fh.writerows((r.astuple() for r in self.records))
