from __future__ import annotations
import bz2
import csv
import concurrent.futures
from dataclasses import dataclass
import json
import logging
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Iterable

from maki.sketch.core import PairwiseResults, cluster_from_pairwise


logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------
# Results Containers
# ---------------------------------------------------------------------

@dataclass(slots=True, frozen=True)
class EdgeSet:
    fieldnames: list[str]
    data: list[dict[str, str]]
    
    @property
    def empty(self) -> bool:
        return len(self.data) == 0
    

class ClusteringResult:
    def __init__(self) -> None:
        self.clusters: list[list[str]] = []
        self.singletons: set[str] = set()
    
    def add_cluster(self, cluster: list[str]):
        self.clusters.append(cluster)
    
    def add_singleton(self, name: str):
        self.singletons.add(name)
    
    def add_singletons(self, values: set[str]):
        self.singletons |= values
    
    @classmethod
    def combine_disjoint(cls, objs: Iterable[ClusteringResult]):
        result = cls()
        
        for obj in objs:
            for cluster in obj.clusters:
                result.add_cluster(cluster)
            
            result.add_singletons(obj.singletons)
            
        return result
    
    @property
    def num_clusters(self) -> int:
        return len(self.clusters) + len(self.singletons)
    
    @property
    def num_singletons(self) -> int:
        return len(self.singletons)
    
    def stats(self) -> tuple[int, int, int]:
        csizes = [len(cluster) for cluster in self.clusters]
        csizes.sort()
        
        mi = 1 if self.singletons else csizes[0] if csizes else 0
        me = 1 if (self.singletons and len(self.singletons) > len(csizes)) else csizes[max((self.num_clusters//2)-self.num_singletons, 0)] if csizes else 0
        ma = csizes[-1] if csizes else 1 if self.singletons else 0
        
        return mi, me, ma


# ---------------------------------------------------------------------
# Thread Workers
# ---------------------------------------------------------------------

def _cluster_worker(param_set: tuple[EdgeSet, float, int, Path]):
    edges, ani, cores, tmp_dir = param_set

    clusters = ClusteringResult()

    # Filter edges to ANI
    filtered_edges, new_singletons = _edge_filter(edges, ani)
    
    # Singletons are nodes that have no edges meeting ANI filter
    clusters.add_singletons(set(new_singletons))
    
    # Cluster connected components
    if not filtered_edges.empty:
        for refined_cluster in _cluster_recs(filtered_edges, ani, cores, tmp_dir):
            clusters.add_cluster(refined_cluster)

    return clusters


def _edge_filter(edges: EdgeSet, min_ani: float):
    population: set[str] = set()
    observed: set[str] = set()
    
    records: list[dict[str, str]] = []
    
    for rec in edges.data:
        population.add(rec["query_name"])
        population.add(rec["match_name"])
        
        if float(rec["average_containment_ani"]) >= min_ani:
            records.append(rec)
            
            observed.add(rec["query_name"])
            observed.add(rec["match_name"])
    
    return EdgeSet(edges.fieldnames, records), [member for member in population if member not in observed]


def _cluster_recs(edges: EdgeSet, ani: float, cores: int, tmp_dir: Path) -> list[list[str]]:
    to_be_clustered = {m for rec in edges.data for m in (rec["query_name"], rec["match_name"])}
    
    with TemporaryDirectory(dir=tmp_dir) as tmp:
        # Export recs to csv
        csv_file = Path(tmp) / "pairwise.csv"
        _export_recs(edges, csv_file)
        
        # Use sourmash pairwise clustering
        clustered = cluster_from_pairwise(csv_file, ani, cores)
        
        # Check for members in the input edge set but not in the output clusters
        missing = to_be_clustered - {member for cluster in clustered for member in cluster}
        
        if missing:
            raise RuntimeError(
                f"{len(missing)} accessions were included in "
                "input edge set but not in final cluster list."
            )
        
        else:
            return clustered


def _export_recs(edges: EdgeSet, file: Path):
    with file.open(mode="wt", newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=edges.fieldnames)
        writer.writeheader()
        writer.writerows(edges.data)


# ---------------------------------------------------------------------
# Data Management
# ---------------------------------------------------------------------

class PairwiseManager:
    """
    Manages a hierarchical clustering structure.

    The current clustering state is defined by all leaf nodes.
    Non-leaf nodes represent previous clustering levels.
    """
    
    THREADS_PER_TASK = 5
    
    def __init__(
        self,
        storage_path: str | Path,
    ):
        self.storage_path = Path(storage_path)
        
        self.fieldnames: list[str] | None = None
        self.data: dict[str, list[dict[str, str]]] = {}
            
    # ---------------------------------------------------------------------
    # Persistence
    # ---------------------------------------------------------------------

    def save(self, empty_ok: bool = True) -> None:
        if not self.fieldnames:
            if empty_ok:
                return
            else:
                raise ValueError("No records to save.")
        
        data = {
            "fieldnames": self.fieldnames,
            "results": self.data
        }

        with bz2.open(self.storage_path, "wt", encoding="utf-8") as fh:
            json.dump(data, fh)
            
    def try_load(self):
        if self.storage_path.exists():
            self._load()

    def _load(self) -> None:
        with bz2.open(self.storage_path, "rt", encoding="utf-8") as fh:
            data = json.load(fh)

        self.fieldnames = data["fieldnames"]        
        self.data: dict[str, list[dict[str, str]]] = data["results"]
        
    # ---------------------------------------------------------------------
    # Insert Methods
    # ---------------------------------------------------------------------

    def group_ids(self):
        yield from self.data

    def insert_results(self, group_id: str, results: PairwiseResults, force: bool = False):
        if group_id in self.data and not force:
            raise ValueError(f"Pairwise results already exist for {group_id}.")
        
        if not self.fieldnames:
            self.fieldnames = results.fieldnames
        
        self.data[group_id] = results.data
        
    # ---------------------------------------------------------------------
    # Clustering Methods
    # ---------------------------------------------------------------------
    
    def refine_clusters(self, cluster_ids: Iterable[str], ani_thresholds: list[float], parallel: int):
        assert self.fieldnames
        
        results: list[ClusteringResult] = []
        
        parallel_workers = (parallel + self.THREADS_PER_TASK - 1) // self.THREADS_PER_TASK
        with concurrent.futures.ProcessPoolExecutor(max_workers=parallel_workers) as executor:
            
            for ani in ani_thresholds:
                logger.info("Clustering at ANI threshold %.4f.", ani)

                tasks = (
                    (EdgeSet(self.fieldnames, self.data[cid]), ani, self.THREADS_PER_TASK, self.storage_path.parent)
                    for cid in cluster_ids
                )

                results.append(
                    ClusteringResult.combine_disjoint(
                        executor.map(
                            _cluster_worker,
                            tasks
                        )
                    )
                )
        
        return results
