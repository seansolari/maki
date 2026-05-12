import csv
from dataclasses import dataclass
from typing import List, Dict
from maki.utils.hashing import hash_strings


@dataclass
class GenomeRecord:
    accession: str
    tax_id: int
    fasta: str
    gff: str


class Manifest:
    def __init__(self, records: List[GenomeRecord]):
        self.records = records

    @classmethod
    def from_csv(cls, path: str) -> "Manifest":
        records = []
        with open(path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                records.append(
                    GenomeRecord(
                        accession=row["accession"],
                        tax_id=int(row["tax_id"]),
                        fasta=row["fasta"],
                        gff=row["gff"]
                    )
                )
        return cls(records)

    def accessions(self):
        return [r.accession for r in self.records]

    def compute_hash(self) -> str:
        return hash_strings(self.accessions())

    def group_by_rank(self, taxonomy, rank: str) -> Dict[int, List[GenomeRecord]]:
        grouped = {}

        for record in self.records:
            ancestor = taxonomy.get_ancestor_at_rank(record.tax_id, rank)
            grouped.setdefault(ancestor, []).append(record)

        return grouped