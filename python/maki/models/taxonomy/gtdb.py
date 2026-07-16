
from pathlib import Path
import re
from typing import Optional

from ete4 import GTDBTaxa
from .base import BaseTaxonomy, TaxidSearchResult, TaxonomyUnresolvedRankException


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
            return TaxidSearchResult(value, True)
        else:
            return self._normalise_name(value)
        
    def _normalise_name(self, name: str):
        assert self.gtdb, "Missing taxonomy data."
        
        name = " ".join(name.split())
        words = name.count(" ") + 1
        
        if words > 1:
            genus, species, *_ = name.split(" ")
            revised = f"s__{genus} {species}"
            return TaxidSearchResult(revised, bool(self.gtdb.get_name_lineage([revised])))
        else:
            uname = name.split(" ")[0]
            for prefix in ("g", "f", "o", "c", "p", "d"):
                prefixed = f"{prefix}__{uname}"
                res = self.gtdb.get_name_lineage([prefixed])
                if res:
                    return TaxidSearchResult(prefixed, True)
            else:
                return TaxidSearchResult(f"u__{uname}", False)

    def get_lineage(self, taxid):
        assert self.gtdb, "Missing taxonomy data."
        
        l = self.gtdb.get_name_lineage([taxid])
        if l:
            return l[0][taxid]
        else:
            uname = re.findall(r"(?:\S__)?(\S+)", taxid)[0]
            for prefix in ("g", "f", "o", "c", "p", "d"):
                prefixed = f"{prefix}__{uname}"
                l = self.gtdb.get_name_lineage([prefixed])
                if l:
                    return l[0][prefixed]
            else:
                return []            

    def get_ancestor_at_rank(self, taxid, rank):
        assert self.gtdb, "Missing taxonomy data."
        
        lineage = self.get_lineage(taxid)
        ranks = self.gtdb.get_rank(lineage)

        for t in lineage:
            if ranks.get(t) == rank:
                return t

        raise TaxonomyUnresolvedRankException(f"No ancestor at rank '{rank}' for taxid {taxid}")
