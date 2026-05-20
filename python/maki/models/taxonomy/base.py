
from abc import abstractmethod
from pathlib import Path
from typing import Dict, Iterable, Tuple


class BaseTaxonomy:
    name: str
  
    def __init__(self, db_root: Path):
        self.db_root = db_root
        self.tax_root = db_root / "taxonomy" / self.name
        self.tax_root.mkdir(parents=True, exist_ok=True)

    @abstractmethod
    def ensure_downloaded(self):
        """Ensure taxonomy database exists locally (no re-download if present)."""
        pass

    @abstractmethod
    def resolve_taxid(self, value: str) -> str:
        """
        Accepts taxid OR name.
        Returns taxid.
        """
        pass

    @abstractmethod
    def get_lineage(self, taxid: str) -> list[str]:
        pass

    @abstractmethod
    def get_ancestor_at_rank(self, taxid: str, rank: str) -> str:
        pass
      
      
def resolve_accession_taxids(records: Iterable[Tuple[str, str]], taxonomy: BaseTaxonomy):
    result: Dict[str, str] = {}

    unresolved = []

    for acc, taxid in records:
        try:
            result[acc] = taxonomy.resolve_taxid(taxid)
        except Exception as e:
            print(f"[warning] {acc}: {e}")
            unresolved.append(acc)

    if unresolved:
        print(f"\n[warning] {len(unresolved)} genomes could not be resolved:")
        for u in unresolved[:10]:
            print(f" - {u}")

    return result
