
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
        
        res = self.gtdb.get_name_lineage([value])

        if not res:
            raise ValueError(f"Could not resolve GTDB name: '{value}'")

        return value

    def get_lineage(self, taxid):
        assert self.gtdb, "Missing taxonomy data."
        
        l = self.gtdb.get_name_lineage([taxid])
        assert l, f"Could not identify taxid: {taxid}"
        print(l)
        
        return l

    def get_ancestor_at_rank(self, taxid, rank):
        assert self.gtdb, "Missing taxonomy data."
        
        lineage = self.get_lineage(taxid)
        ranks = self.gtdb.get_rank(lineage)

        for t in lineage:
            if ranks.get(t) == rank:
                return t

        raise ValueError(f"No ancestor at rank '{rank}' for taxid {taxid}")