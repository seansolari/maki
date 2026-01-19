
"""Cluster file generation helpers.

Builds a representative-to-member mapping for annotated genes across genomes.
"""
from __future__ import annotations
from pathlib import Path
from typing import Mapping, Sequence

def write_cluster_map(
    path: Path,
    representative_genome: str,
    accession_map: Mapping[str, str],
    gene_ids: Sequence[str],
) -> None:
    """Write the cluster TSV used by the Random Gene Distance pipeline."""
    path.parent.mkdir(parents=True, exist_ok=True)
    rep_acc = accession_map[representative_genome]
    with open(path, "w") as fh:
        for gid in gene_ids:
            rep = f"{rep_acc}-gene-{gid}"
            for g, acc in accession_map.items():
                mem = f"{acc}-gene-{gid}"
                fh.write(f"{rep}\t{mem}\n")
