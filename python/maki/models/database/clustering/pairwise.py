
import bz2
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Iterable, Optional

from maki.models.database.sketch.core import PairwiseResults, cluster_from_pairwise


class PairwiseManager:
    """
    Manages a hierarchical clustering structure.

    The current clustering state is defined by all leaf nodes.
    Non-leaf nodes represent previous clustering levels.
    """

    def __init__(
        self,
        storage_path: str | Path,
    ):
        self.storage_path = Path(storage_path)
        
        self.fieldnames: list[str] | None = None
        self.data: dict[str, list[dict[str, str]]] = {}

        if self.storage_path.exists():
            self._load()
            
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
    
    def cluster_recursive(self, group_id: str, ani_thresholds: list[float], cores: int):
        if not self.fieldnames:
            raise RuntimeError("No pairwise data for recursive clustering.")
        
        data = self.data[group_id]
        
        ani_thresholds.sort()
        sub_clusters: list[list[list[str]]] = []
        
        for ani in ani_thresholds:
            if sub_clusters:
                new_clustering: list[list[str]] = []
                
                for cluster in sub_clusters[-1]:
                    if len(cluster) > 1:
                        filtered_data = filter(
                            lambda rec: rec["query_name"] in cluster and rec["match_name"] in cluster,
                            data
                        )
                        
                        for sub_cluster in self._cluster_recs(
                            filtered_data,
                            ani,
                            cores
                        ):
                            
                            new_clustering.append(sub_cluster)
                    
                    else:
                        new_clustering.append(cluster)
                
                sub_clusters.append(new_clustering)
            
            else:
                sub_clusters.append(self._cluster_recs(data, ani, cores))
        
        return sub_clusters
            
    def _cluster_recs(self, data: Iterable[dict[str, str]], ani: float, cores: int) -> list[list[str]]:
        with TemporaryDirectory(dir=self.storage_path.parent) as tmp:
            # Export recs to csv
            csv_file = Path(tmp) / "pairwise.csv"
            self._export_recs(data, csv_file)
            
            # Use sourmash pairwise clustering
            return cluster_from_pairwise(csv_file, ani, cores)
    
    def _export_recs(self, data: Iterable[dict[str, str]], file: Path):
        assert self.fieldnames
        
        with file.open(mode="wt", newline='', encoding='utf-8') as fh:
            writer = csv.DictWriter(fh, fieldnames=self.fieldnames)
            writer.writeheader()
            writer.writerows(data)
        