from __future__ import annotations

from pathlib import Path
import subprocess
import sys
from typing import List, Sequence

from maki.models.database.sketch.core import zip_signatures
from maki.models.database.sketch.store import SourmashSketchStore
from maki.utils.hashing import hash_strings

from .manager import ClusterRefinementOperator


class SourmashClusterRefiner(ClusterRefinementOperator):
    def __init__(
        self,
        ani: float,
        sketch_db: SourmashSketchStore,
        working_dir: Path
    ) -> None:
        super().__init__(f"ani-{ani:.3f}")
        
        self.ani = ani
        self.sketch_db = sketch_db
        
        self.working_dir = working_dir
    
    def refine(self, accessions: Sequence[str]) -> List[List[str]]:
        grp_name = hash_strings(accessions)[:8]
        sigs = self.sketch_db.load_many(list(accessions))

        # Combine signatures into ZIP file
        zar_file = self.working_dir / f"{grp_name}.zip"
        
        zip_signatures(
            sigs.values(),
            zar_file
        )
        
        # Pairwise comparisons within group
        stat_file = self.working_dir / f"{grp_name}.csv"
        
        cmd = [
            sys.executable, "-m",
            "sourmash", "scripts", "pairwise",
            zar_file,
            "-o", stat_file
        ]
        subprocess.run(cmd, check=True)
        
        # Cluster results
        cluster_file = self.working_dir / f"{grp_name}-clusters.csv"
        cluster_sizes_file = self.working_dir / f"{grp_name}-cluster-sizes.csv"
        cmd = [
            sys.executable, "-m",
            "sourmash", "scripts", "cluster",
            stat_file,
            "--similarity_column", "average_containment_ani",
            "-t", str(self.ani),
            "-o", cluster_file,
            "--cluster-sizes", cluster_sizes_file
        ]
        subprocess.run(cmd, check=True)
        