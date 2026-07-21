
from enum import Enum
from pathlib import Path
from typing import Iterable, List, Mapping, Optional

from ete4 import GTDBTaxa, update_ete_data
from .base import BaseTaxonomy, TaxidNotFound, TaxonomyUnresolvedRankException


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

    def __init__(self, db_root: Path, releases: Optional[GTDBRelease | Iterable[GTDBRelease]] = None):
        super().__init__(db_root)
        
        if isinstance(releases, GTDBRelease):
            releases = [releases]
        elif releases is None:
            releases = []
        else:
            releases = list(releases)
        
        releases += self.list_releases()
        if not releases:
            releases.append(GTDBRelease.latest)
        
        self.gtdb: Mapping[GTDBRelease, GTDBTaxa] = {rel: self._ensure_release(rel) for rel in releases}
        
    def list_releases(self):
        for p in self.tax_root.iterdir():
            try:
                yield GTDBRelease[p.name]
            except KeyError:
                continue
            
    def _ensure_release(self, rel: GTDBRelease) -> GTDBTaxa:
        rel_base = self.tax_root / rel.name
        rel_base.mkdir(parents=True, exist_ok=True)
        
        sql = rel_base / f"{rel.value}.sqlite"
        if not sql.exists():
            local_dmp = (rel_base / f"{rel.value}dump.tar.gz").absolute()
            remote_dmp = rel.ete_dmp_url()
            update_ete_data(str(local_dmp), str(remote_dmp))
            GTDBTaxa(dbfile=str(sql), taxdump_file=str(local_dmp))
            local_dmp.unlink()
            
        return GTDBTaxa(dbfile=str(sql))
    
    def search_taxid(self, value: str):
        for rel, db in self.gtdb.items():
            res = db.get_name_lineage([value])
            if res:
                yield rel, value

    def resolve_taxid(self, value):
        assert self.gtdb, "Missing taxonomy data."
        value = value.strip()

        if any(1 for _ in self.search_taxid(value)):
            return value
        else:
            for fixed_value in self._normalise_name(value):
                if any(1 for _ in self.search_taxid(fixed_value)):
                    return fixed_value
                
        raise TaxidNotFound(value)        
        
    def _normalise_name(self, name: str):
        assert self.gtdb, "Missing taxonomy data."
        
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
        assert self.gtdb, "Missing taxonomy data."
        
        lineages: List[List[str]] = []
        for db in self.gtdb.values():
            l = db.get_name_lineage([taxid])
            if l:
                lineages.append(l[0][taxid])
        
        if lineages:
            l = lineages.pop(max(range(len(lineages)), key = lambda i: len(lineages[i])))
            
            for other_lineage in lineages:
                if not all(j in l for j in other_lineage):
                    raise ValueError(f"Inconsistent lineages for taxid={taxid}: {",".join(l)} != {','.join(other_lineage)}")
            else:
                return l
        else:
            raise TaxidNotFound(taxid)

    def get_ancestor_at_rank(self, taxid, rank):
        assert self.gtdb, "Missing taxonomy data."
        
        lineage = self.get_lineage(taxid)
        lineage_ranks = {rel: db.get_rank(lineage) for rel, db in self.gtdb.items()}
        
        result: Mapping[GTDBRelease, str] = {}
        for rel, ranks in lineage_ranks.items():
            for t in lineage:
                if ranks.get(t) == rank:
                    result[rel] = t
                    break
        
        if not result:
            raise TaxonomyUnresolvedRankException(f"No ancestor at rank '{rank}' for taxid {taxid}")
        else:
            ancestors = set(result.values())
            if len(ancestors) > 1:
                raise TaxonomyUnresolvedRankException(f"Inconsistent ancestors between taxonomy releases: {", ".join(f'{rel.name}={v}' for rel, v in result.items())}")
            else:
                return ancestors.pop()
