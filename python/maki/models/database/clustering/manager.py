from __future__ import annotations
import bz2
import json
from dataclasses import asdict, dataclass
from pathlib import Path
import random
from typing import Iterable, List, Optional
from .pairwise import ClusteringResult


@dataclass
class ClusterNode:
    cluster_id: str
    accessions: List[str]
    
    @property
    def is_singleton(self) -> bool:
        return len(self.accessions) == 1
    
    def add_member(self, accession: str):
        self.accessions.append(accession)


class ClusterManager:
    def __init__(
        self,
        storage_path: str | Path,
    ):
        self.storage_path = Path(storage_path)

        if self.storage_path.exists():
            self._load()
        
        else:
            self._data: dict[str, dict[str, ClusterNode]] = {}

    @property
    def tags(self) -> Iterable[str]:
        return self._data.keys()
    
    def has_tag(self, qry: str) -> bool:
        return qry in self._data

    def save(self) -> None:        
        obj = {
            tag: {
                cluster_id: asdict(cluster)
                for cluster_id, cluster in tag_data.items()
            }
            for tag, tag_data in self._data.items()
        }

        with bz2.open(self.storage_path, "wt", encoding="utf-8") as fh:
            json.dump(obj, fh, indent=2)

    def _load(self) -> None:
        with bz2.open(self.storage_path, "rt", encoding="utf-8") as fh:
            data = json.load(fh)

        self._data = {
            tag: {
                cluster_id: ClusterNode(**cluster_data)
                for cluster_id, cluster_data in tag_data.items()
            }
            for tag, tag_data in data.items()
        }

    def iter_tag(self, tag: str):
        yield from self._data[tag].values()
        
    def get_cluster(self, tag: str, cluster_id: str):
        return self._data[tag][cluster_id]

    def add_anon_clustering(
        self,
        obj: ClusteringResult,
        tag: str
    ):
        if tag in self._data:
            raise ValueError(f"Tag {tag} already exists")
        
        self._data[tag] = {}
        
        next_id = 0
        
        for cluster in obj.clusters:
            cluster_id = str(next_id)
            next_id += 1

            self._data[tag][cluster_id] = ClusterNode(cluster_id, cluster)

        for singleton in obj.singletons:
            cluster_id = str(next_id)
            next_id += 1
            
            self._data[tag][cluster_id] = ClusterNode(cluster_id, [singleton])
    
    def add_named_clustering(
        self,
        obj: dict[str, str],
        *,
        tag: str
    ):
        if tag in self._data:
            raise ValueError(f"Tag {tag} already exists")
        
        self._data[tag] = {}
        
        for cluster_member, cluster_id in obj.items():
            self._data[tag].setdefault(cluster_id, ClusterNode(cluster_id, [])).add_member(cluster_member)
            
    def group_by(
        self,
        *,
        tag: str,
        derep_tag: Optional[str] = None,
        seed: int = 42
    ):
        rng = random.Random()
        rng.seed(seed)
        
        groups = {
            cx.cluster_id: cx.accessions
            for cx in self.iter_tag(tag)
        }
        
        if derep_tag:
            derep_map = {
                accn: cx.cluster_id
                for cx in self.iter_tag(derep_tag)
                for accn in cx.accessions
            }
            
            for cid in sorted(groups.keys()):
                # Group by de-replication cluster
                subgroups: dict[str, list[str]] = {}
                
                for mem in groups[cid]:
                    subgroups.setdefault(derep_map[mem], []).append(mem)
                
                # Sort to enforce reproducible selection per seed
                for subgroup in subgroups:
                    subgroups[subgroup].sort()
                
                # Choose one per de-replication cluster
                groups[cid] = [
                    rng.choice(subgroups[subgroup])
                    for subgroup in sorted(subgroups.keys())
                ]

        return groups
