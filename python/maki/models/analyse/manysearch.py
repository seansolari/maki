from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Iterable, Iterator

from maki.sketch.core import build_standalone_manifest

from .sketch import MetagenomeSketch


class ManySearchStore:
    """
    Incrementally manage sourmash branchwater manysearch containment results.

    Layout:

        root/
        ├── SRR1/
        │   └── genomes_abcd1234.json.gz
        ├── SRR2/
        │   └── genomes_abcd1234.json.gz
        └── ...

    Each file contains a JSON list of manysearch records for a
    single metagenome accession.
    """

    def __init__(
        self,
        root: str | Path,
        rocksdb: str | Path,
    ) -> None:
        self.root = Path(root)
        self.rocksdb = Path(rocksdb)

        self.root.mkdir(parents=True, exist_ok=True)

    @property
    def db_id(self) -> str:
        """
        Stable identifier for the database.

        Uses path + modification time so results are invalidated
        when the database is rebuilt.
        """
        stat = self.rocksdb.stat()

        h = hashlib.sha256()
        h.update(str(self.rocksdb.resolve()).encode())
        h.update(str(stat.st_mtime_ns).encode())

        return h.hexdigest()[:8]

    @property
    def result_filename(self) -> str:
        return f"{self.rocksdb.stem}_{self.db_id}.json.gz"

    def accession_dir(self, accession: str) -> Path:
        return self.root / accession

    def result_path(self, accession: str) -> Path:
        return self.accession_dir(accession) / self.result_filename

    def has_result(self, accession: str) -> bool:
        return self.result_path(accession).exists()

    def list_accessions(self) -> list[str]:
        accessions = []

        for d in self.root.iterdir():
            if not d.is_dir():
                continue

            if self.result_path(d.name).exists():
                accessions.append(d.name)

        return sorted(accessions)

    def load_results(self, accession: str) -> list:
        with gzip.open(self.result_path(accession), "rt") as fp:
            return json.load(fp)

    def iter_results(self) -> Iterator[tuple[str, list[dict]]]:
        for accession in self.list_accessions():
            yield accession, self.load_results(accession)

    def remove(self, accession: str) -> None:
        self.result_path(accession).unlink(missing_ok=True)

    def missing(
        self,
        metagenomes: Iterable[MetagenomeSketch],
    ) -> list[MetagenomeSketch]:
        return [
            mg
            for mg in metagenomes
            if not self.has_result(mg.accession)
        ]

    # Perform search
    # --------------

    def update(
        self,
        metagenomes: Iterable[MetagenomeSketch],
        *,
        threshold_bp: int | None = None,
        threads: int | None = None,
    ) -> list[str]:
        """
        Run manysearch only for metagenomes lacking results.

        Returns the list of newly processed accessions.
        """

        pending = self.missing(metagenomes)

        if not pending:
            return []

        self._run_manysearch(
            pending,
            threshold_bp=threshold_bp,
            threads=threads,
        )

        return [x.accession for x in pending]

    def _run_manysearch(
        self,
        metagenomes: list[MetagenomeSketch],
        *,
        threshold_bp: int | None,
        threads: int | None,
    ) -> None:

        with tempfile.TemporaryDirectory() as tmpdir:

            tmpdir = Path(tmpdir)

            # Build manifest
            manifest = tmpdir / "metagenomes.csv"
            
            build_standalone_manifest(
                (m.accession for m in metagenomes),
                manifest
            )
            
            # Perform search
            output = tmpdir / "manysearch.csv"

            cmd = [
                sys.executable, "-m",
                "sourmash",
                "scripts",
                "manysearch",
                str(self.rocksdb),
                str(manifest),
                "-o", str(output),
            ]

            if threads is not None:
                cmd.extend(
                    [
                        "--threads",
                        str(threads),
                    ]
                )

            if threshold_bp is not None:
                cmd.extend(
                    [
                        "--threshold-bp",
                        str(threshold_bp),
                    ]
                )

            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
            )

            if proc.returncode:
                raise RuntimeError(
                    "manysearch failed\n\n"
                    f"Command:\n{' '.join(cmd)}\n\n"
                    f"STDOUT:\n{proc.stdout}\n\n"
                    f"STDERR:\n{proc.stderr}"
                )

            self._split_results(
                output_csv=output,
                metagenomes=metagenomes,
            )

    def _split_results(
        self,
        output_csv: Path,
        metagenomes: list[MetagenomeSketch],
    ) -> None:
        """
        Split manysearch output into per-accession result files.

        Also creates empty result files for metagenomes
        with no matches.
        """

        per_accession: dict[str, list[dict]] = {
            mg.accession: []
            for mg in metagenomes
        }

        with output_csv.open(newline="") as fp:
            reader = csv.DictReader(fp)

            for row in reader:
                query = row["match_name"]
                per_accession.setdefault(query, []).append(row)

        for accession, records in per_accession.items():
            self._store_accession_result(
                accession,
                records,
            )

    def _store_accession_result(
        self,
        accession: str,
        records: list[dict],
    ) -> None:

        acc_dir = self.accession_dir(accession)
        acc_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        final_path = self.result_path(accession)

        with tempfile.NamedTemporaryFile(
            dir=acc_dir,
            delete=False,
            suffix=".tmp",
        ) as tmp:

            tmp_path = Path(tmp.name)

        try:
            with gzip.open(tmp_path, "wt") as fp:
                json.dump(records, fp)

            shutil.move(tmp_path, final_path)

        finally:
            tmp_path.unlink(missing_ok=True)

    def refresh(
        self,
        metagenome: MetagenomeSketch,
        *,
        threshold_bp: int | None = None,
        threads: int | None = None,
    ) -> None:
        self.remove(metagenome.accession)

        self.update(
            [metagenome],
            threshold_bp=threshold_bp,
            threads=threads,
        )
