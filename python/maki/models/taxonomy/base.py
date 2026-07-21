
from abc import abstractmethod
from enum import Enum
import logging
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple


logger = logging.getLogger(__name__)


class Ranks(str, Enum):
    Species = "species"
    Genus = "genus"
    Family = "family"
    Order = "order"
    Class = "class"
    Phylum = "phylum"
    SuperKingdom = "superkingdom"
    
    
class TaxidNotFound(Exception):
    pass


class BaseTaxonomy:
    name: str
  
    def __init__(self, db_root: Path):
        self.db_root = db_root
        self.tax_root = db_root / "taxonomy" / self.name
        self.tax_root.mkdir(parents=True, exist_ok=True)

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
      
      
class TaxonomyUnresolvedRankException(Exception):
    pass
      
      
def resolve_accession_taxids(records: Iterable[Tuple[str, str]], taxonomy: BaseTaxonomy):
    result: Dict[str, str] = {}
    unresolved: Set[str] = set()

    for acc, taxid in records:
        try:
            result[acc] = taxonomy.resolve_taxid(taxid)
        except TaxidNotFound:
            logger.warning(f"Could not resolve {acc} taxid {taxid}")
            unresolved.add(taxid)
    
    if unresolved:
        unresolved_l = list(unresolved)
        raise TaxidNotFound(f"{';'.join(unresolved_l[:5])}{'...' if len(unresolved_l) > 5 else ''}")

    return result
