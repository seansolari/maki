
from __future__ import annotations
from dataclasses import dataclass
from enum import Enum
import json
from pathlib import Path
from typing import Dict, List

from maki.models.taxonomy import BaseTaxonomy, GTDBTaxonomy, NCBITaxonomy
from maki.models.taxonomy.base import resolve_accession_taxids
from .archive import XzArchive, RawArchive
from .manifest import DatabasePackage, GenomeRecord, Manifest, ManifestSchema, read_manifest


class UpdateMode(Enum):
    fixed = 0
    updateable = 1
    

class TaxonomySource(Enum):
    gtdb = 0
    ncbi = 1


@dataclass(frozen=True)
class DatabaseOptions:
    root: Path
    kmer_size: int
    rank: str
    update_mode: UpdateMode
    tax: TaxonomySource


class _Database:
    _SCHEMA = ManifestSchema("accession", "taxonomy", "taxid")
    
    def __init__(self, root: Path, kmer_size: int, rank: str, update_mode: UpdateMode, taxonomy: BaseTaxonomy, manifest: Manifest[GenomeRecord]):
        self.root = root
        self.kmer_size = kmer_size
        self.rank = rank
        self.update_mode = update_mode
        self.taxonomy = taxonomy
        self.manifest = manifest
        
    def save_manifest(self):
        self.manifest.save(self.root / "manifest.csv.gz")
        
    @classmethod
    def read_manifest(cls, root: Path):
        return read_manifest(root / "manifest.csv.gz", cls._SCHEMA)


class StaticDatabase(_Database):
    def __init__(self, root: Path, kmer_size: int, rank: str, update_mode: UpdateMode, clusters: XzArchive, taxonomy: BaseTaxonomy, manifest: Manifest[GenomeRecord]) -> None:
        super().__init__(root, kmer_size, rank, update_mode, taxonomy, manifest)
        self.metadata_file = root / "metadata.json"
        self.clusters = clusters

    @classmethod
    def load(cls, root: Path):
        assert root.exists()
        
        meta = json.loads((root / "metadata.json").read_text())
        
        return StaticDatabase(root, meta["kmer_size"], meta["rank"], UpdateMode[meta["update_mode"]], XzArchive(root / "clusters.tar.xz"), cls._init_taxonomy(root, TaxonomySource[meta["taxonomy_source"]]), cls.read_manifest(root))
  
    @classmethod
    def create(cls, db_otps: DatabaseOptions):
        db_otps.root.mkdir(parents=True, exist_ok=True)
        
        clusters = XzArchive(db_otps.root / "clusters.tar.xz")
        taxonomy = cls._init_taxonomy(db_otps.root, db_otps.tax)
        manifest = Manifest[GenomeRecord]()
        
        result = cls(db_otps.root, db_otps.kmer_size, db_otps.rank, db_otps.update_mode, clusters, taxonomy, manifest)
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
      
    def preserve(self):
        meta = {
            "kmer_size": self.kmer_size,
            "rank": self.rank,
            "update_mode": self.update_mode.name,
            "taxonomy_source": "gtdb" if isinstance(self.taxonomy, GTDBTaxonomy) else "ncbi",
        }

        with self.metadata_file.open("w") as f:
            json.dump(meta, f, indent=2)
            
    def decompress(self) -> WriteableDatabase:
        return WriteableDatabase(self.root, self.kmer_size, self.rank, self.update_mode, self.clusters.decompress(), self.taxonomy, self.manifest)
    

class WriteableDatabase(_Database):
    def __init__(self, root: Path, kmer_size: int, rank: str, update_mode: UpdateMode, clusters: RawArchive, taxonomy: BaseTaxonomy, manifest: Manifest[GenomeRecord]) -> None:
        super().__init__(root, kmer_size, rank, update_mode, taxonomy, manifest)
        self.clusters = clusters
        
    def compress(self) -> StaticDatabase:
        return StaticDatabase(self.root, self.kmer_size, self.rank, self.update_mode, self.clusters.compress(), self.taxonomy, self.manifest)
      
    def insert(self, package: DatabasePackage, concurrency: int):
        taxids = self._get_taxids(package)
        groups = self._group_by_rank(package, taxids)
        
        for cluster_id, accessions in groups:
            # fetch cluster
            new_cluster = cluster_id not in self.clusters
            
            print(f"[insert] retrieving cluster {cluster_id}")
            cluster = self.clusters.get_or_create(cluster_id)
            
            if (not new_cluster) and (not cluster.updateable):
                raise RuntimeError(f"Attempted insertion to non-updateable cluster: {cluster_id}")
            
            # build cluster
            with package.retrieve_data(accessions) as data:
                print(f"[insert] inserting sequences into cluster {cluster.root}")
                cluster.insert(data.genomes())
                
                print(f"[insert] constructing cluster index {cluster.root}")
                cluster.build(self.kmer_size, concurrency)
                
                if cluster.updateable or self.update_mode == UpdateMode.updateable:
                    print(f"[insert] persisting sources at {cluster.source_dir}")
                    cluster.persist_sources()
                else:
                    print(f"[insert] clearing sources at {cluster.source_dir}")
                    cluster.remove_sources()
                
                # update manifest
                for rec in data.genomes():
                    if not rec.taxid:
                        rec.taxid = taxids.get(rec.accession)
                    
                    self.manifest.insert(rec)
        
        # preserve manifest
        self.save_manifest()
        
    def _get_taxids(self, package: DatabasePackage):
        return resolve_accession_taxids(((r.accession, r.taxonomy) for r in package.records() if not r.taxid), self.taxonomy)
        
    def _group_by_rank(self, package: DatabasePackage, taxids: Dict[str, str]):
        grouped: Dict[str, List[str]] = {}
        
        for record in package.records():
            taxid = taxids.get(record.accession)
            ancestor = record.taxonomy if not taxid else self.taxonomy.get_ancestor_at_rank(taxid, self.rank)
            
            key = "/".join(self.taxonomy.get_lineage(ancestor))
            grouped.setdefault(key, []).append(record.accession)

        return list(grouped.items())
