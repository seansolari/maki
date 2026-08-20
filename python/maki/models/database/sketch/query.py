from __future__ import annotations
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import csv
from dataclasses import asdict
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import sourmash

from .core import (
    FORMAT_VERSION,
    ContainmentHit,
    SketchParameters,
    downsample_signature,
    load_one_signature,
    new_minhash,
    save_one_signature,
    validate_signature
)


def _merge_cluster_worker(
    args: tuple[int, list[str], SketchParameters, str],
) -> tuple[int, str, int]:
    cluster_id, sketch_paths, params, output_dir = args

    union = new_minhash(params)

    for sketch_path in sketch_paths:
        signature = load_one_signature(sketch_path)
        validate_signature(signature, params, allow_finer_scaled=False)
        union.add_many(signature.minhash)

    signature = sourmash.SourmashSignature(
        union,
        name=f"cluster-{cluster_id:09d}",
        filename="",
    )

    output_path = (
        Path(output_dir) /
        f"cluster-{cluster_id:09d}.sig"
    )
    save_one_signature(signature, output_path)

    return cluster_id, str(output_path), len(union)
  

class ClusterContainmentDatabase:
    """
    Build and query an SBT containing cluster-union sketches.

    The SBT is constructed through sourmash's own module entry point rather
    than private SBT implementation classes. This is intentional: the
    command/index format is substantially more stable than the internal
    GraphFactory, SigLeaf, storage, and node APIs.

    sourmash remains the only non-standard runtime dependency.
    """

    def __init__(self, database_dir: str | Path):
        self.database_dir = Path(database_dir)
        self.cluster_sketch_dir = self.database_dir / "cluster_sketches"
        self.metadata_path = self.database_dir / "database_metadata.json"
        self.cluster_table_path = self.database_dir / "clusters.tsv"
        self.sbt_path = self.database_dir / "clusters.sbt.zip"

        self._metadata: dict | None = None

        if self.metadata_path.exists():
            with self.metadata_path.open("rt") as fp:
                self._metadata = json.load(fp)

    @classmethod
    def build(
        cls,
        phase1_dir: str | Path,
        database_dir: str | Path,
        *,
        processes: int | None = None,
        delete_genome_sketches: bool = False,
        sbt_threshold: float = 0.01,
        overwrite: bool = False,
    ) -> "ClusterContainmentDatabase":
        phase1_dir = Path(phase1_dir)
        database_dir = Path(database_dir)

        phase1_metadata_path = phase1_dir / "build_metadata.json"
        assignments_path = phase1_dir / "genome_clusters.tsv"

        if not phase1_metadata_path.exists():
            raise FileNotFoundError(phase1_metadata_path)
        if not assignments_path.exists():
            raise FileNotFoundError(assignments_path)

        if database_dir.exists() and any(database_dir.iterdir()):
            if not overwrite:
                raise FileExistsError(
                    f"{database_dir} is non-empty; use overwrite=True"
                )
            shutil.rmtree(database_dir)

        database_dir.mkdir(parents=True, exist_ok=True)
        db = cls(database_dir)
        db.cluster_sketch_dir.mkdir(parents=True, exist_ok=True)

        with phase1_metadata_path.open("rt") as fp:
            phase1_metadata = json.load(fp)

        params = SketchParameters(**phase1_metadata["sketch"])
        params.validate()

        clusters: dict[int, list[str]] = defaultdict(list)
        genome_ids: dict[int, list[str]] = defaultdict(list)

        with assignments_path.open("rt", newline="") as fp:
            reader = csv.DictReader(fp, delimiter="\t")

            for row in reader:
                cluster_id = int(row["cluster_id"])
                sketch_path = row["sketch_path"]

                if not Path(sketch_path).exists():
                    raise FileNotFoundError(sketch_path)

                clusters[cluster_id].append(sketch_path)
                genome_ids[cluster_id].append(row["genome_id"])

        processes = processes or os.cpu_count() or 1

        merge_jobs = [
            (
                cluster_id,
                sketch_paths,
                params,
                str(db.cluster_sketch_dir),
            )
            for cluster_id, sketch_paths in sorted(clusters.items())
        ]

        with ProcessPoolExecutor(max_workers=processes) as executor:
            merged = list(executor.map(_merge_cluster_worker, merge_jobs))

        cluster_paths = {
            cluster_id: sketch_path
            for cluster_id, sketch_path, _ in merged
        }
        cluster_hash_counts = {
            cluster_id: hash_count
            for cluster_id, _, hash_count in merged
        }

        with db.cluster_table_path.open("wt", newline="") as fp:
            writer = csv.writer(fp, delimiter="\t")
            writer.writerow(
                [
                    "cluster_id",
                    "cluster_name",
                    "cluster_sketch_path",
                    "genome_count",
                    "hash_count",
                    "genome_ids_json",
                ]
            )

            for cluster_id in sorted(clusters):
                writer.writerow(
                    [
                        cluster_id,
                        f"cluster-{cluster_id:09d}",
                        cluster_paths[cluster_id],
                        len(genome_ids[cluster_id]),
                        cluster_hash_counts[cluster_id],
                        json.dumps(genome_ids[cluster_id]),
                    ]
                )

        # Build the SBT using sourmash's supported index entry point.
        #
        # Calling Python via sys.executable ensures the sourmash module comes
        # from the same environment as this package.
        cmd = [
            sys.executable,
            "-m",
            "sourmash",
            "index",
            str(db.sbt_path),
            *[
                cluster_paths[cluster_id]
                for cluster_id in sorted(cluster_paths)
            ],
            "-k",
            str(params.ksize),
        ]

        subprocess.run(cmd, check=True)

        metadata = {
            "format_version": FORMAT_VERSION,
            "phase": "cluster_containment_database",
            "sourmash_version": getattr(
                sourmash,
                "__version__",
                "unknown",
            ),
            "sketch": asdict(params),
            "cluster_count": len(clusters),
            "source_phase1_dir": str(phase1_dir.resolve()),
            "sbt_path": str(db.sbt_path.resolve()),
            "cluster_table_path": str(db.cluster_table_path.resolve()),
            "default_query_threshold": sbt_threshold,
        }

        with db.metadata_path.open("wt") as fp:
            json.dump(metadata, fp, indent=2, sort_keys=True)

        if delete_genome_sketches:
            # Delete only after all cluster signatures, the SBT, the cluster
            # table, and metadata have been written successfully.
            for sketch_paths in clusters.values():
                for sketch_path in sketch_paths:
                    Path(sketch_path).unlink(missing_ok=True)

            genome_sketch_dir = phase1_dir / "genome_sketches"
            try:
                genome_sketch_dir.rmdir()
            except OSError:
                # Leave it alone if it contains unrelated files.
                pass

        db._metadata = metadata
        return db

    @property
    def metadata(self) -> dict:
        if self._metadata is None:
            if not self.metadata_path.exists():
                raise FileNotFoundError(
                    f"No database metadata at {self.metadata_path}"
                )

            with self.metadata_path.open("rt") as fp:
                self._metadata = json.load(fp)

        return self._metadata # type: ignore

    @property
    def sketch_params(self) -> SketchParameters:
        return SketchParameters(**self.metadata["sketch"])

    def validate_query(self, query: sourmash.SourmashSignature) -> None:
        validate_signature(
            query,
            self.sketch_params,
            allow_finer_scaled=True,
        )

    def query(
        self,
        query: sourmash.SourmashSignature | str | Path,
        *,
        threshold: float | None = None,
        max_results: int | None = None,
    ) -> list:
        """Return cluster containment hits at or above threshold.

        Containment direction is:

            fraction of cluster sketch contained in query sketch

        This is the useful direction for asking whether the sample contains
        evidence for a reference cluster.
        """
        if isinstance(query, (str, Path)):
            query_signature = load_one_signature(query)
        else:
            query_signature = query

        self.validate_query(query_signature)

        params = self.sketch_params
        query_signature = downsample_signature(
            query_signature,
            params.scaled,
        )

        if threshold is None:
            threshold = float(
                self.metadata.get("default_query_threshold", 0.01)
            )

        if not 0.0 <= threshold <= 1.0:
            raise ValueError("threshold must be between 0 and 1")

        with tempfile.TemporaryDirectory(
            prefix="cluster-prefilter-"
        ) as tmp_dir:
            query_path = Path(tmp_dir) / "query.sig"
            csv_path = Path(tmp_dir) / "results.csv"

            save_one_signature(query_signature, query_path)

            cmd = [
                sys.executable,
                "-m",
                "sourmash",
                "search",
                str(query_path),
                str(self.sbt_path),
                "--containment",
                "--threshold",
                str(threshold),
                "--csv",
                str(csv_path),
                "-k",
                str(params.ksize),
                "--quiet",
            ]

            subprocess.run(cmd, check=True)

            hits = self._read_search_results(csv_path)

        hits.sort(
            key=lambda hit: (-hit.containment, hit.cluster_id)
        )

        if max_results is not None:
            if max_results < 0:
                raise ValueError("max_results must be non-negative")
            hits = hits[:max_results]

        return hits

    def _read_search_results(
        self,
        csv_path: str | Path,
    ) -> list:
        hits: list[ContainmentHit] = []

        if not Path(csv_path).exists():
            return hits

        with Path(csv_path).open("rt", newline="") as fp:
            reader = csv.DictReader(fp)

            for row in reader:
                name = row.get("name") or row.get("match_name") or ""
                cluster_id = self._parse_cluster_id(name)

                # sourmash versions have used slightly different CSV column
                # labels. Try the known containment/similarity alternatives.
                containment_text = (
                    row.get("containment")
                    or row.get("similarity")
                    or row.get("score")
                )

                if containment_text is None:
                    raise ValueError(
                        "Could not find a containment/similarity column in "
                        f"sourmash output: {reader.fieldnames}"
                    )

                hits.append(
                    ContainmentHit(
                        cluster_id=cluster_id,
                        containment=float(containment_text),
                    )
                )

        return hits

    @staticmethod
    def _parse_cluster_id(name: str) -> int:
        prefix = "cluster-"

        if not name.startswith(prefix):
            raise ValueError(
                f"Unexpected cluster signature name: {name!r}"
            )

        value = name[len(prefix):].split()[0]
        return int(value)

    def cluster_metadata(self) -> dict[int, dict]:
        result: dict[int, dict] = {}

        with self.cluster_table_path.open("rt", newline="") as fp:
            reader = csv.DictReader(fp, delimiter="\t")

            for row in reader:
                cluster_id = int(row["cluster_id"])
                result[cluster_id] = {
                    "cluster_name": row["cluster_name"],
                    "cluster_sketch_path": row["cluster_sketch_path"],
                    "genome_count": int(row["genome_count"]),
                    "hash_count": int(row["hash_count"]),
                    "genome_ids": json.loads(row["genome_ids_json"]),
                }

        return result
