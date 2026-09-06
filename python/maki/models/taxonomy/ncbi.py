
import logging
from pathlib import Path
from typing import Optional

from ete4 import NCBITaxa
from .base import BaseTaxonomy, TaxidNotFound


logger = logging.getLogger(__name__)


class NCBITaxonomy(BaseTaxonomy):
    name = "ncbi"

    def __init__(self, db_root: Path):
        super().__init__(db_root)
        
        sqlfile = self.root / "ncbi.sqlite"
        if not sqlfile.exists():
            logger.error("Uninitialised taxonomy: %s", db_root)
            raise RuntimeError("Uninitialised taxonomy: %s" % db_root)
        
        self.ncbi = NCBITaxa(dbfile=str(sqlfile))
        
    # Configuring taxonomy
    # --------------------
        
    @staticmethod
    def ensure_release(root: Path, release: Optional[str] = None) -> None:
        sqlfile = root / "ncbi.sqlite"
        
        if release:
            taxdmp = Path(release)
            
            if not taxdmp.exists():
                logger.error("Taxdump does not exist: %s", release)
                raise RuntimeError("Taxdump does not exist: %s" % release)
            
            NCBITaxa(dbfile=str(sqlfile), taxdump_file=taxdmp)
            
        else:
            NCBITaxa(dbfile=str(sqlfile))
    
    # Lookup
    # ------

    def resolve_taxid(self, value):
        if value.isdigit():
            name = self.ncbi.get_taxid_translator([value])
            if not name:
                raise TaxidNotFound(value)
            
            return value

        # resolve name → taxid
        res = self.ncbi.get_name_translator([value])
        if not res:
            raise TaxidNotFound(value)

        return res[value][0]

    def get_lineage(self, taxid):
        l = self.ncbi.get_lineage(taxid)
        assert l, f"Could not identify taxid: {taxid}"
        
        return [str(v) for v in l]

    def get_ancestor_at_rank(self, taxid, rank):
        lineage = self.get_lineage(taxid)
        ranks = self.ncbi.get_rank(lineage)

        for t in lineage:
            if ranks.get(t) == rank:
                return t

        raise ValueError(f"No ancestor at rank '{rank}' for taxid {taxid}")
