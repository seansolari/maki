
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
        self.dbfile = self.root / "ncbi.sqlite"
        self.ncbi: Optional[NCBITaxa] = None

    def ensure_downloaded(self):
        if self.dbfile.exists():
            logger.info("Using existing NCBI taxonomy.")
        else:
            logger.info("Downloading NCBI taxonomy (ETE4)...")

        self.ncbi = NCBITaxa(dbfile=str(self.dbfile))

    def resolve_taxid(self, value):
        assert self.ncbi, "Missing taxonomy data."
      
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
        assert self.ncbi, "Missing taxonomy data."
        
        l = self.ncbi.get_lineage(taxid)
        assert l, f"Could not identify taxid: {taxid}"
        
        return [str(v) for v in l]

    def get_ancestor_at_rank(self, taxid, rank):
        assert self.ncbi, "Missing taxonomy data."
        
        lineage = self.get_lineage(taxid)
        ranks = self.ncbi.get_rank(lineage)

        for t in lineage:
            if ranks.get(t) == rank:
                return t

        raise ValueError(f"No ancestor at rank '{rank}' for taxid {taxid}")
