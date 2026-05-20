
from __future__ import annotations
from enum import Enum
import json
from pathlib import Path

from maki.models.taxonomy import BaseTaxonomy, GTDBTaxonomy, NCBITaxonomy, resolve_accession_taxids
from .archive import XzArchive, RawArchive
from .manifest import Manifest


class TaxonomySource(Enum):
    gtdb = 0
    ncbi = 1
    
    @classmethod
    def from_str(cls, val: str) -> TaxonomySource:
        try:
            return cls[val.lower()]
        except KeyError:
            raise ValueError(f"Unrecognised taxonomy source: {val}")


class StaticDatabase:
    def __init__(self, root: Path, kmer_size: int, rank: str, clusters: XzArchive, taxonomy: BaseTaxonomy, manifest: Manifest) -> None:
        self.root = root
        self.metadata_file = root / "metadata.json"
        self.kmer_size = kmer_size
        self.rank = rank
        self.clusters = clusters
        self.taxonomy = taxonomy
        self.manifest = manifest
  
    @classmethod
    def load(cls, root: Path):
        assert root.exists()
        
        meta = json.loads((root / "metadata.json").read_text())
        
        return StaticDatabase(root, meta["kmer_size"], meta["rank"], XzArchive(root / "clusters.tar.xz"), cls._init_taxonomy(root, TaxonomySource.from_str(meta["taxonomy_source"])), Manifest.load(root / "manifest.csv.gz"))
  
    @classmethod
    def create(cls, root: Path, kmer_size: int, rank: str, tax: TaxonomySource, manifest_src: Path):
        root = root
        root.mkdir(parents=True, exist_ok=True)
        
        clusters = XzArchive(root / "clusters.tar.xz")
        taxonomy = cls._init_taxonomy(root, tax)
        manifest = Manifest.import_csv(manifest_src)
        
        result = cls(root, kmer_size, rank, clusters, taxonomy, manifest)
        result.update_taxids()
        result.preserve()
        
        return result
    
    @staticmethod
    def _init_taxonomy(root: Path, source: TaxonomySource) -> BaseTaxonomy:
        if source == TaxonomySource.ncbi:
            taxonomy = NCBITaxonomy(root)
        else:
            taxonomy = GTDBTaxonomy(root)

        taxonomy.ensure_downloaded()
        return taxonomy
      
    def update_taxids(self):
        print("Resolving genome taxonomy...")
        mapping = resolve_accession_taxids(((r.accession, r.taxonomy) for r in self.manifest.records if not r.taxid), self.taxonomy)

        # attach taxids to records
        for rec in self.manifest.records:
            if not rec.taxid:
                rec.taxid = mapping.get(rec.accession)

        # remove unresolved
        valid_records = sum(1 for r in self.manifest.records if r.taxid)
        missing_records = len(self.manifest.records) - valid_records

        if missing_records:
            print(f"[warning] {missing_records} genomes could not be assigned taxonomy IDs")

        print(f"{valid_records} genomes with valid taxonomy.")
      
    def preserve(self):
        meta = {
            "kmer_size": self.kmer_size,
            "rank": self.rank,
            "taxonomy_source": "gtdb" if isinstance(self.taxonomy, GTDBTaxonomy) else "ncbi",
        }

        with self.metadata_file.open("w") as f:
            json.dump(meta, f, indent=2)
            
        self.manifest.save(self.root / "manifest.csv.gz")
            
    def decompress(self) -> WriteableDatabase:
        return WriteableDatabase(self.root, self.kmer_size, self.rank, self.clusters.decompress(), self.taxonomy, self.manifest)
    

class WriteableDatabase:
    def __init__(self, root: Path, kmer_size: int, rank: str, clusters: RawArchive, taxonomy: BaseTaxonomy, manifest: Manifest) -> None:
        self.root = root
        self.kmer_size = kmer_size
        self.rank = rank
        self.clusters = clusters
        self.taxonomy = taxonomy
        self.manifest = manifest
        
    def compress(self) -> StaticDatabase:
        return StaticDatabase(self.root, self.kmer_size, self.rank, self.clusters.compress(), self.taxonomy, self.manifest)
