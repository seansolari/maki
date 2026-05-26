
from difflib import get_close_matches
from pathlib import Path
from typing import Optional

from ete4 import GTDBTaxa
from .base import BaseTaxonomy


class GTDBTaxonomy(BaseTaxonomy):
    name = "gtdb"

    def __init__(self, db_root: Path):
        super().__init__(db_root)
        self.gtdb: Optional[GTDBTaxa] = None

    def ensure_downloaded(self):
        print("Initializing GTDB taxonomy (ETE4)...")
        self.gtdb = GTDBTaxa(dbfile=str(self.tax_root / "gtdb.sqlite"))

    def resolve_taxid(self, value):
        assert self.gtdb, "Missing taxonomy data."
        value = value.strip()

        # Direct exact match
        res = self.gtdb.get_name_lineage([value])
        if res:
            return value
          
        # Check prefixes
        for prefixed in self._check_hierarchy(value):
            res = self.gtdb.get_name_lineage([prefixed])
            if res:
                return prefixed

        raise ValueError(f"Could not resolve GTDB taxonomy for: '{value}'")

    def get_lineage(self, taxid):
        assert self.gtdb, "Missing taxonomy data."
        
        l = self.gtdb.get_name_lineage([taxid])
        assert l, f"Could not identify taxid: {taxid}"
        
        return l[0][taxid]

    def get_ancestor_at_rank(self, taxid, rank):
        assert self.gtdb, "Missing taxonomy data."
        
        lineage = self.get_lineage(taxid)
        ranks = self.gtdb.get_rank(lineage)

        for t in lineage:
            if ranks.get(t) == rank:
                return t

        raise ValueError(f"No ancestor at rank '{rank}' for taxid {taxid}")
      
    ## Taxonomy helpers
    
    @staticmethod
    def _check_hierarchy(name: str):
        name = " ".join(name.split())
        words = name.count(" ") + 1
        
        if words > 1:
            genus, species, *_ = name.split(" ")
            yield f"s__{genus} {species}"
        
        uname = name.split(" ")[0]
        for prefix in ("g", "f", "o", "c", "p", "d"):
            yield f"{prefix}__{uname}"
