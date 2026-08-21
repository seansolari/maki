from __future__ import annotations
from collections import defaultdict
from dataclasses import asdict
from typing import Mapping

from sourmash import SourmashSignature

from .core import FORMAT_VERSION, ClusterAssignment, SketchParameters, downsample_signature
from .store import SourmashSketchStore


class GenomeClusterBuilder:
    """
    Streaming, approximate genome clustering.

    Algorithm
    ---------
    1. Sketch genomes in parallel.
    2. Process sketches in input order.
    3. Use banded hashes to retrieve likely cluster representatives.
    4. Compare only those representatives using sourmash Jaccard.
    5. Assign to the best representative above the threshold.
    6. Otherwise create a new cluster.

    Cluster representatives are fixed: the first genome assigned to a cluster.
    This prevents progressive union sketches from making clusters increasingly
    broad and reduces memory consumption.
    """

    def __init__(
        self,
        *,
        ani_threshold: float = 0.95,
        downsample_factor: int = 10,
        max_candidates: int = 128,
    ):
        if not 0.0 < ani_threshold <= 1.0:
            raise ValueError("ani_threshold must be in (0, 1]")
        if downsample_factor < 1:
            raise ValueError("downsample_factor must be >=1")
        if max_candidates <= 0:
            raise ValueError("max_candidates must be positive")

        self.ani_threshold = ani_threshold
        self.downsample_factor = downsample_factor
        self.max_candidates = max_candidates

    def cluster(
        self,
        sketches: SourmashSketchStore,
    ) -> tuple[list[ClusterAssignment], dict]:
        
        # cluster_id -> representative signature
        representatives: dict[int, SourmashSignature] = {}

        # key -> cluster IDs
        reverse_index: dict[int, int] = {}

        records: list[ClusterAssignment] = []
        next_cluster_id = 0

        for genome_id in sketches.accessions():
            signature = sketches.load_signature(genome_id)

            candidates, signature_lsh = self._candidate_clusters(
                signature,
                reverse_index
            )

            best_cluster: int | None = None
            best_ani = -1.0

            for cluster_id in candidates:
                representative = representatives[cluster_id]

                # Use the sourmash comparison API rather than constructing an
                # explicit pairwise matrix.
                jaccard = signature.minhash.jaccard(
                    representative.minhash,
                )
                ani = self.jaccard_to_ani(
                    jaccard,
                    sketches.params.ksize,
                )

                if ani >= self.ani_threshold and ani > best_ani:
                    best_cluster = cluster_id
                    best_ani = ani

            if best_cluster is None:
                best_cluster = next_cluster_id
                next_cluster_id += 1

                representatives[best_cluster] = signature
                self._index_representative(
                    best_cluster,
                    signature_lsh,
                    reverse_index,
                )

            records.append(
                ClusterAssignment(
                    genome_id=genome_id,
                    cluster_id=best_cluster,
                )
            )

        return records, self._summarise_metadata(records, sketches.params)

    def _candidate_clusters(
        self,
        signature: SourmashSignature,
        reverse_index: Mapping[int, int],
    ) -> tuple[list[int], list[int]]:
        subsig = downsample_signature(signature, self.downsample_factor * signature.minhash.scaled)
        subhashes = list(subsig.minhash.hashes)
        
        counts: dict[int, int] = defaultdict(int)

        for hash_value in subhashes:
            try:
                counts[reverse_index[hash_value]] += 1
            except KeyError:
                continue

        # Highest number of matching bands first. Capping this bound keeps an
        # exceptionally repetitive genome from causing a very large comparison
        # burst.
        top_matches = sorted(
            filter(lambda cluster_id: cluster_id > -1, counts),
            key=lambda cluster_id: (-counts[cluster_id], cluster_id),
        )
        return top_matches[: self.max_candidates], subhashes

    def _index_representative(
        self,
        cluster_id: int,
        hashes: list[int],
        reverse_index: dict[int, int],
    ) -> None:
        for hash_value in hashes:
            if hash_value in reverse_index:
                reverse_index[hash_value] = -1
            else:
                reverse_index[hash_value] = cluster_id

    @staticmethod
    def jaccard_to_ani(jaccard: float, ksize: int) -> float:
        """
        Mash-style ANI approximation from Jaccard.

        J = x / (2 - x), where x = ANI**k.
        Therefore x = 2J / (1 + J).
        """
        if jaccard <= 0:
            return 0.0
        if jaccard >= 1:
            return 1.0

        x = (2.0 * jaccard) / (1.0 + jaccard)
        return x ** (1.0 / ksize)

    def _summarise_metadata(
        self,
        assignments: list[ClusterAssignment],
        params: SketchParameters
    ) -> dict:
        n_clusters = 1 + max(r.cluster_id for r in assignments)

        metadata = {
            "format_version": FORMAT_VERSION,
            "phase": "genome_clustering",
            "sketch": asdict(params),
            "clustering": {
                "method": "streaming_greedy_lsh",
                "ani_threshold": self.ani_threshold,
                "downsample_factor": self.downsample_factor,
                "max_candidates": self.max_candidates,
                "representative": "first_genome",
            },
            "genome_count": len(assignments),
            "cluster_count": n_clusters,
        }

        return metadata
            