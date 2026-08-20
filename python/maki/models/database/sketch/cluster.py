from __future__ import annotations
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import csv
from dataclasses import asdict
import json
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence

from sourmash import MinHash, SourmashSignature

from .core import (
    FORMAT_VERSION,
    SketchParameters,
    load_one_signature,
    new_minhash,
    safe_filename,
    save_one_signature,
    validate_signature
)
from maki.utils.io import read_fasta


def _sketch_worker(args: tuple[GenomeInput, SketchParameters, str]):
    genome, params, output_dir = args

    sketch_name = safe_filename(genome.genome_id) + ".sig"
    sketch_path = Path(output_dir) / sketch_name

    mh = new_minhash(params)

    for _, sequence in read_fasta(genome.fasta_path):
        # force=True skips k-mers containing ambiguous characters rather than
        # failing the complete genome.
        mh.add_sequence(sequence, force=True)

    if len(mh) == 0:
        raise ValueError(
            f"Genome {genome.genome_id!r} produced an empty sketch"
        )

    signature = SourmashSignature(
        mh,
        name=genome.genome_id,
        filename=str(genome.fasta_path),
    )
    save_one_signature(signature, sketch_path)

    return genome.genome_id, genome.fasta_path, str(sketch_path)


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
        output_dir: str | Path,
        *,
        sketch_params: SketchParameters = SketchParameters(),
        ani_threshold: float = 0.95,
        lsh_hashes: int = 64,
        band_size: int = 4,
        max_candidates: int = 128,
        processes: int | None = None,
    ):
        sketch_params.validate()

        if not 0.0 < ani_threshold <= 1.0:
            raise ValueError("ani_threshold must be in (0, 1]")
        if lsh_hashes <= 0:
            raise ValueError("lsh_hashes must be positive")
        if band_size <= 0:
            raise ValueError("band_size must be positive")
        if lsh_hashes < band_size:
            raise ValueError("lsh_hashes must be >= band_size")
        if max_candidates <= 0:
            raise ValueError("max_candidates must be positive")

        self.output_dir = Path(output_dir)
        self.genome_sketch_dir = self.output_dir / "genome_sketches"
        self.params = sketch_params
        self.ani_threshold = ani_threshold
        self.lsh_hashes = lsh_hashes
        self.band_size = band_size
        self.max_candidates = max_candidates
        self.processes = processes or os.cpu_count() or 1

    @property
    def metadata_path(self) -> Path:
        return self.output_dir / "build_metadata.json"

    @property
    def assignments_path(self) -> Path:
        return self.output_dir / "genome_clusters.tsv"

    def build(
        self,
        genomes: Iterable[GenomeInput | tuple[str, str]],
    ) -> list:
        genome_inputs = [
            item if isinstance(item, GenomeInput)
            else GenomeInput(genome_id=item[0], fasta_path=item[1])
            for item in genomes
        ]

        if not genome_inputs:
            raise ValueError("No genomes were supplied")

        genome_ids = [g.genome_id for g in genome_inputs]
        if len(genome_ids) != len(set(genome_ids)):
            raise ValueError("genome_id values must be unique")

        self.genome_sketch_dir.mkdir(parents=True, exist_ok=True)

        sketch_jobs = [
            (genome, self.params, str(self.genome_sketch_dir))
            for genome in genome_inputs
        ]

        with ProcessPoolExecutor(max_workers=self.processes) as executor:
            sketched = list(executor.map(_sketch_worker, sketch_jobs))

        assignments = self._cluster_sketches(sketched)
        self._write_outputs(assignments)

        return assignments

    def _cluster_sketches(
        self,
        sketched: Sequence[tuple[str, str, str]],
    ) -> list:
        # cluster_id -> representative signature
        representatives: dict[int, SourmashSignature] = {}

        # band key -> cluster IDs
        band_index: dict[tuple[int, ...], list[int]] = defaultdict(list)

        records: list[GenomeRecord] = []
        next_cluster_id = 0

        for genome_id, fasta_path, sketch_path in sketched:
            signature = load_one_signature(sketch_path)
            validate_signature(
                signature,
                self.params,
                allow_finer_scaled=False,
            )

            candidates = self._candidate_clusters(
                signature.minhash,
                band_index,
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
                    self.params.ksize,
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
                    signature.minhash,
                    band_index,
                )

            records.append(
                GenomeRecord(
                    genome_id=genome_id,
                    fasta_path=fasta_path,
                    sketch_path=sketch_path,
                    cluster_id=best_cluster,
                )
            )

        return records

    def _candidate_clusters(
        self,
        mh: MinHash,
        band_index: Mapping[tuple[int, ...], list[int]],
    ) -> list:
        counts: dict[int, int] = defaultdict(int)

        for key in self._band_keys(mh):
            for cluster_id in band_index.get(key, []):
                counts[cluster_id] += 1

        # Highest number of matching bands first. Capping this bound keeps an
        # exceptionally repetitive genome from causing a very large comparison
        # burst.
        ordered = sorted(
            counts,
            key=lambda cluster_id: (-counts[cluster_id], cluster_id),
        )
        return ordered[: self.max_candidates]

    def _index_representative(
        self,
        cluster_id: int,
        mh: MinHash,
        band_index: dict[tuple[int, ...], list[int]],
    ) -> None:
        for key in self._band_keys(mh):
            band_index[key].append(cluster_id)

    def _band_keys(self, mh: MinHash) -> Iterator[tuple[int, ...]]:
        selected = sorted(mh.hashes)[: self.lsh_hashes]

        # Very small or highly repetitive assemblies may have fewer available
        # hashes. They still receive a key, although discrimination is lower.
        if len(selected) < self.band_size:
            if selected:
                yield tuple(selected)
            return

        limit = len(selected) - self.band_size + 1

        # Non-overlapping bands reduce index size. Offset is included to avoid
        # treating an identical tuple at distinct positions as the same band.
        for start in range(0, limit, self.band_size):
            band = selected[start:start + self.band_size]
            yield (start, *band)

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

    def _write_outputs(self, assignments: Sequence[GenomeRecord]) -> None:
        self.output_dir.mkdir(parents=True, exist_ok=True)

        with self.assignments_path.open("wt", newline="") as fp:
            writer = csv.writer(fp, delimiter="\t")
            writer.writerow(
                ["genome_id", "fasta_path", "sketch_path", "cluster_id"]
            )

            for record in assignments:
                writer.writerow(
                    [
                        record.genome_id,
                        record.fasta_path,
                        record.sketch_path,
                        record.cluster_id,
                    ]
                )

        n_clusters = 1 + max(r.cluster_id for r in assignments)

        metadata = {
            "format_version": FORMAT_VERSION,
            "phase": "genome_clustering",
            "sketch": asdict(self.params),
            "clustering": {
                "method": "streaming_greedy_lsh",
                "ani_threshold": self.ani_threshold,
                "lsh_hashes": self.lsh_hashes,
                "band_size": self.band_size,
                "max_candidates": self.max_candidates,
                "representative": "first_genome",
            },
            "genome_count": len(assignments),
            "cluster_count": n_clusters,
        }

        with self.metadata_path.open("wt") as fp:
            json.dump(metadata, fp, indent=2, sort_keys=True)
            