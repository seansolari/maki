
from enum import Enum
import logging
from pathlib import Path
from typing import Optional

from ete4 import GTDBTaxa, update_ete_data
from .base import BaseTaxonomy, TaxidNotFound, TaxonomyUnresolvedRankException


logger = logging.getLogger(__name__)


class GTDBRelease(str, Enum):
    r95 = "gtdb95"
    r202 = "gtdb202"
    r207 = "gtdb207"
    r214 = "gtdb214"
    r220 = "gtdb220"
    r226 = "gtdb226"
    latest = "gtdblatest"
    
    def ete_dmp_url(self):
        return f"gtdb_taxonomy/{self.value}/{'gtdb_latest_' if self.name == "latest" else self.value}dump.tar.gz"


class GTDBTaxonomy(BaseTaxonomy):
    name = "gtdb"

    def __init__(self, db_root: Path):
        super().__init__(db_root)
        
        sqlfile = db_root / "gtdb.sqlite"
        if not sqlfile.exists():
            logger.error("Uninitialised taxonomy: %s", db_root)
            raise RuntimeError("Uninitialised taxonomy: %s" % db_root)
        
        self.gtdb = GTDBTaxa(dbfile=str(sqlfile))
        
    # Configuring taxonomy
    # --------------------
    
    @staticmethod
    def ensure_release(root: Path, release: Optional[str] = None) -> None:
        if release:
            if release.endswith(".tar.gz"):
                taxdmp = Path(release)
                
                if not taxdmp.exists():
                    logger.error("Taxdump does not exist: %s", release)
                    raise RuntimeError("Taxdump does not exist: %s" % release)
                
                GTDBTaxonomy.ensure_local_release(root, taxdmp)
                
            else:
                try:
                    gtdb_version = GTDBRelease[release]
                except KeyError:
                    logger.error("Unrecognised GTDB release: %s", release)
                    raise RuntimeError("Unrecognised GTDB release: %s" % release)
                else:
                    GTDBTaxonomy.ensure_external_release(root, gtdb_version)

        else:
            GTDBTaxonomy.ensure_external_release(root, GTDBRelease.latest)
    
    @staticmethod
    def ensure_local_release(root: Path, taxdump: Path):
        rel_name = taxdump.name.removesuffix(".tar.gz")
        rel_base = root / rel_name
        rel_base.mkdir(parents=True, exist_ok=True)
        
        sql = rel_base / f"{rel_name}.sqlite"
        if not sql.exists():
            GTDBTaxa(dbfile=str(sql), taxdump_file=str(taxdump))

    @staticmethod
    def ensure_external_release(root: Path, rel: GTDBRelease):
        rel_base = root / rel.name
        rel_base.mkdir(parents=True, exist_ok=True)
        
        sql = rel_base / f"{rel.value}.sqlite"
        if not sql.exists():
            local_dmp = (rel_base / f"{rel.value}dump.tar.gz").absolute()
            remote_dmp = rel.ete_dmp_url()
            update_ete_data(str(local_dmp), str(remote_dmp))
            GTDBTaxa(dbfile=str(sql), taxdump_file=str(local_dmp))
            local_dmp.unlink()
    
    # Lookup
    # ------
    
    def search_taxid(self, value: str) -> Optional[str]:
        res = self.gtdb.get_name_lineage([value])
        return value if res else None

    def resolve_taxid(self, value):
        value = value.strip()

        if self.search_taxid(value):
            return value
        else:
            for fixed_value in self._normalise_name(value):
                if self.search_taxid(fixed_value):
                    return fixed_value
        
        raise TaxidNotFound(value)
        
    def _normalise_name(self, name: str):
        name = " ".join(name.split())
        words = name.count(" ") + 1
        
        if words > 1:
            genus, species, *_ = name.split(" ")
            yield f"s__{genus} {species}"
        else:
            uname = name.split(" ")[0]
            for prefix in ("g", "f", "o", "c", "p", "d"):
                yield f"{prefix}__{uname}"

    def get_lineage(self, taxid):
        l: list[str] = self.gtdb.get_name_lineage([taxid])
        if l:
            return l
        else:
            raise TaxidNotFound(taxid)

    def get_ancestor_at_rank(self, taxid, rank):
        lineage = self.get_lineage(taxid)
        ranks = self.gtdb.get_rank(lineage)
    
        for t in lineage:
            if ranks.get(t) == rank:
                return t
        else:
            raise TaxonomyUnresolvedRankException(f"No ancestor at rank '{rank}' for taxid {taxid}")
