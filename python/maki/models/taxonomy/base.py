
from abc import abstractmethod
from dataclasses import dataclass
from enum import Enum
import logging
from pathlib import Path
from typing import Dict, Iterable, Tuple


logger = logging.getLogger(__name__)


class Ranks(str, Enum):
    Species = "species"
    Genus = "genus"
    Family = "family"
    Order = "order"
    Class = "class"
    Phylum = "phylum"
    Domain = "domain"


@dataclass(frozen=True)
class TaxidSearchResult:
    result: str
    found: bool


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
    def resolve_taxid(self, value: str) -> TaxidSearchResult:
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
      
      
class TaxonomyUnresolvedRankException(Exception):
    pass
      
      
def resolve_accession_taxids(records: Iterable[Tuple[str, str]], taxonomy: BaseTaxonomy):
    result: Dict[str, str] = {}

    for acc, taxid in records:
        qry = taxonomy.resolve_taxid(taxid)
        if not qry.found:
            logger.warning(f"Could not resolve {acc} taxid {taxid}, using {qry.result}")
        result[acc] = qry.result

    return result
