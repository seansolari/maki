from __future__ import annotations
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import os
from pathlib import Path
import sqlite3
from typing import Iterable, Iterator, Literal, Sequence
import uuid


ManifestState = Literal["reserved", "complete", "failed"]


def utc_now() -> str:
    """Return an ISO 8601 UTC timestamp suitable for storing in SQLite."""
    return datetime.now(timezone.utc).isoformat()


@dataclass(frozen=True, slots=True)
class ManifestRecord:
    accession: str
    taxid: int
    state: ManifestState
    package_id: str
    reserved_at: str
    completed_at: str | None
    error: str | None


class ManifestError(RuntimeError):
    """Base exception for manifest operations."""


class AccessionsAlreadyExistError(ManifestError):
    def __init__(self, accessions: Iterable[str]):
        self.accessions = tuple(sorted(accessions))
        joined = ", ".join(self.accessions)
        super().__init__(f"Accessions already exist in the manifest: {joined}")


class SQLiteManifest:
    """
    Transactional accession manifest backed by a standalone SQLite file.

    Connections are deliberately short-lived so this object remains safe and
    lightweight when hooks are instantiated in many independent processes.

    Do not share an individual sqlite3.Connection between processes.
    """

    SCHEMA_VERSION = 1

    def __init__(
        self,
        path: str | Path,
        *,
        timeout: float = 60.0,
    ) -> None:
        self.path = Path(path)
        self.timeout = timeout

        if not self.path.is_file():
            raise FileNotFoundError(
                f"Manifest database does not exist: {self.path}"
            )

        self._validate_schema()

    @classmethod
    def create(
        cls,
        path: str | Path,
        *,
        timeout: float = 60.0,
    ) -> SQLiteManifest:
        path = Path(path)

        if path.exists():
            raise FileExistsError(f"Manifest already exists: {path}")

        path.parent.mkdir(parents=True, exist_ok=True)

        connection = sqlite3.connect(
            path,
            timeout=timeout,
            isolation_level=None,
        )

        try:
            # WAL allows readers and a writer to operate concurrently.
            connection.execute("PRAGMA journal_mode = WAL")
            connection.execute("PRAGMA synchronous = FULL")
            connection.execute("PRAGMA foreign_keys = ON")

            connection.executescript(
                """
                BEGIN IMMEDIATE;

                CREATE TABLE manifest_metadata (
                    key   TEXT PRIMARY KEY,
                    value TEXT NOT NULL
                );

                CREATE TABLE accessions (
                    accession   TEXT PRIMARY KEY,
                    taxid       TEXT NOT NULL,
                    state       TEXT NOT NULL
                                CHECK (state IN (
                                    'reserved',
                                    'complete',
                                    'failed'
                                )),
                    package_id  TEXT NOT NULL,
                    reserved_at TEXT NOT NULL,
                    completed_at TEXT,
                    error       TEXT
                );

                CREATE INDEX accessions_taxid_idx
                    ON accessions(taxid);

                CREATE INDEX accessions_state_idx
                    ON accessions(state);

                CREATE INDEX accessions_package_idx
                    ON accessions(package_id);

                INSERT INTO manifest_metadata(key, value)
                VALUES ('schema_version', '1');

                COMMIT;
                """
            )
        except Exception:
            try:
                connection.execute("ROLLBACK")
            except sqlite3.Error:
                pass

            connection.close()

            # Avoid leaving a partially initialized database behind.
            path.unlink(missing_ok=True)
            raise
        else:
            connection.close()

        return cls(path, timeout=timeout)

    def _connect(self) -> sqlite3.Connection:
        connection = sqlite3.connect(
            self.path,
            timeout=self.timeout,
            isolation_level=None,
        )
        connection.row_factory = sqlite3.Row

        connection.execute("PRAGMA foreign_keys = ON")
        connection.execute(f"PRAGMA busy_timeout = {int(self.timeout * 1000)}")

        return connection

    def _validate_schema(self) -> None:
        with self._connect() as connection:
            row = connection.execute(
                """
                SELECT value
                FROM manifest_metadata
                WHERE key = 'schema_version'
                """
            ).fetchone()

        if row is None:
            raise ManifestError(
                f"Manifest has no schema version: {self.path}"
            )

        try:
            version = int(row["value"])
        except (TypeError, ValueError) as exc:
            raise ManifestError("Invalid manifest schema version") from exc

        if version != self.SCHEMA_VERSION:
            raise ManifestError(
                "Unsupported manifest schema version "
                f"{version}; expected {self.SCHEMA_VERSION}"
            )

    def reserve(
        self,
        records: Sequence[tuple[str, str]],
        *,
        package_id: str | None = None,
        fail_if_exists: bool = True,
    ) -> tuple[str, tuple[str, ...]]:
        """
        Atomically reserve accessions.

        Returns:
            A tuple containing:
              - the package ID
              - accessions newly reserved by this call

        If fail_if_exists is True, the entire operation fails when any
        accession is already present.

        If fail_if_exists is False, existing accessions are skipped and only
        previously unseen accessions are reserved.
        """
        normalized = self._validate_records(records)
        package_id = package_id or str(uuid.uuid4())

        if not normalized:
            return package_id, ()

        requested = {accession for accession, _ in normalized}
        placeholders = ",".join("?" for _ in requested)

        connection = self._connect()

        try:
            # BEGIN IMMEDIATE obtains the write reservation before checking
            # uniqueness. Another writer waits instead of racing between
            # SELECT and INSERT.
            connection.execute("BEGIN IMMEDIATE")

            rows = connection.execute(
                f"""
                SELECT accession
                FROM accessions
                WHERE accession IN ({placeholders})
                """,
                tuple(requested),
            ).fetchall()

            existing = {row["accession"] for row in rows}

            if existing and fail_if_exists:
                raise AccessionsAlreadyExistError(existing)

            new_records = [
                (accession, taxid)
                for accession, taxid in normalized
                if accession not in existing
            ]

            reserved_at = utc_now()

            connection.executemany(
                """
                INSERT INTO accessions (
                    accession,
                    taxid,
                    state,
                    package_id,
                    reserved_at,
                    completed_at,
                    error
                )
                VALUES (?, ?, 'reserved', ?, ?, NULL, NULL)
                """,
                [
                    (accession, taxid, package_id, reserved_at)
                    for accession, taxid in new_records
                ],
            )

            connection.execute("COMMIT")

            return (
                package_id,
                tuple(accession for accession, _ in new_records),
            )

        except Exception:
            try:
                connection.execute("ROLLBACK")
            except sqlite3.Error:
                pass
            raise
        finally:
            connection.close()

    def mark_complete(
        self,
        accessions: Iterable[str],
        *,
        package_id: str,
    ) -> None:
        self._set_state(
            accessions,
            package_id=package_id,
            state="complete",
        )

    def mark_failed(
        self,
        accessions: Iterable[str],
        *,
        package_id: str,
        error: str,
    ) -> None:
        self._set_state(
            accessions,
            package_id=package_id,
            state="failed",
            error=error,
        )

    def _set_state(
        self,
        accessions: Iterable[str],
        *,
        package_id: str,
        state: ManifestState,
        error: str | None = None,
    ) -> None:
        accession_list = tuple(dict.fromkeys(accessions))

        if not accession_list:
            return

        connection = self._connect()

        try:
            connection.execute("BEGIN IMMEDIATE")

            completed_at = utc_now() if state == "complete" else None

            cursor = connection.executemany(
                """
                UPDATE accessions
                SET state = ?,
                    completed_at = ?,
                    error = ?
                WHERE accession = ?
                  AND package_id = ?
                  AND state = 'reserved'
                """,
                [
                    (
                        state,
                        completed_at,
                        error,
                        accession,
                        package_id,
                    )
                    for accession in accession_list
                ],
            )

            # sqlite3's executemany rowcount is supported for UPDATE, although
            # validating the final state explicitly gives clearer diagnostics.
            placeholders = ",".join("?" for _ in accession_list)

            rows = connection.execute(
                f"""
                SELECT accession, state, package_id
                FROM accessions
                WHERE accession IN ({placeholders})
                """,
                accession_list,
            ).fetchall()

            by_accession = {row["accession"]: row for row in rows}

            invalid = []

            for accession in accession_list:
                row = by_accession.get(accession)

                if (
                    row is None
                    or row["package_id"] != package_id
                    or row["state"] != state
                ):
                    invalid.append(accession)

            if invalid:
                raise ManifestError(
                    "Could not transition these accessions to "
                    f"{state!r}: {', '.join(sorted(invalid))}"
                )

            connection.execute("COMMIT")

        except Exception:
            try:
                connection.execute("ROLLBACK")
            except sqlite3.Error:
                pass
            raise
        finally:
            connection.close()

    def get(self, accession: str) -> ManifestRecord | None:
        with self._connect() as connection:
            row = connection.execute(
                """
                SELECT
                    accession,
                    taxid,
                    state,
                    package_id,
                    reserved_at,
                    completed_at,
                    error
                FROM accessions
                WHERE accession = ?
                """,
                (accession,),
            ).fetchone()

        return self._row_to_record(row) if row is not None else None

    def iter_records(
        self,
        *,
        state: ManifestState | None = None,
    ) -> Iterator:
        query = """
            SELECT
                accession,
                taxid,
                state,
                package_id,
                reserved_at,
                completed_at,
                error
            FROM accessions
        """
        params: tuple[object, ...] = ()

        if state is not None:
            query += " WHERE state = ?"
            params = (state,)

        query += " ORDER BY accession"

        # Keep the connection alive for the duration of iteration.
        connection = self._connect()

        try:
            cursor = connection.execute(query, params)

            for row in cursor:
                yield self._row_to_record(row)
        finally:
            connection.close()

    def count(self, *, state: ManifestState | None = None) -> int:
        with self._connect() as connection:
            if state is None:
                row = connection.execute(
                    "SELECT COUNT(*) AS count FROM accessions"
                ).fetchone()
            else:
                row = connection.execute(
                    """
                    SELECT COUNT(*) AS count
                    FROM accessions
                    WHERE state = ?
                    """,
                    (state,),
                ).fetchone()

        return int(row["count"])

    def export_taxonomy_tsv(
        self,
        destination: str | Path,
        *,
        complete_only: bool = True,
        include_header: bool = False,
        overwrite: bool = False,
    ) -> Path:
        """
        Export accession-to-taxid mappings.

        Output columns:
            accession <TAB> taxid

        The file is written to a temporary sibling and atomically renamed into
        place, so readers never observe a partially written export.
        """
        destination = Path(destination)

        if destination.exists() and not overwrite:
            raise FileExistsError(destination)

        destination.parent.mkdir(parents=True, exist_ok=True)

        query = "SELECT accession, taxid FROM accessions"
        params: tuple[object, ...] = ()

        if complete_only:
            query += " WHERE state = ?"
            params = ("complete",)

        query += " ORDER BY accession"

        temporary = self._temporary_sibling(destination)

        try:
            with self._connect() as connection, temporary.open(
                "w",
                encoding="utf-8",
                newline="",
            ) as output:
                writer = csv.writer(
                    output,
                    delimiter="\t",
                    lineterminator="\n",
                )

                if include_header:
                    writer.writerow(("accession", "taxid"))

                cursor = connection.execute(query, params)

                for row in cursor:
                    writer.writerow((row["accession"], row["taxid"]))

                output.flush()
                os.fsync(output.fileno())

            os.replace(temporary, destination)
        except Exception:
            temporary.unlink(missing_ok=True)
            raise

        return destination

    def export_manifest_tsv(
        self,
        destination: str | Path,
        *,
        complete_only: bool = False,
        include_header: bool = True,
        overwrite: bool = False,
    ) -> Path:
        """
        Export the complete audit manifest, including state information.
        """
        destination = Path(destination)

        if destination.exists() and not overwrite:
            raise FileExistsError(destination)

        destination.parent.mkdir(parents=True, exist_ok=True)

        query = """
            SELECT
                accession,
                taxid,
                state,
                package_id,
                reserved_at,
                completed_at,
                error
            FROM accessions
        """
        params: tuple[object, ...] = ()

        if complete_only:
            query += " WHERE state = ?"
            params = ("complete",)

        query += " ORDER BY accession"

        temporary = self._temporary_sibling(destination)

        try:
            with self._connect() as connection, temporary.open(
                "w",
                encoding="utf-8",
                newline="",
            ) as output:
                writer = csv.writer(
                    output,
                    delimiter="\t",
                    lineterminator="\n",
                )

                if include_header:
                    writer.writerow(
                        (
                            "accession",
                            "taxid",
                            "state",
                            "package_id",
                            "reserved_at",
                            "completed_at",
                            "error",
                        )
                    )

                for row in connection.execute(query, params):
                    writer.writerow(
                        (
                            row["accession"],
                            row["taxid"],
                            row["state"],
                            row["package_id"],
                            row["reserved_at"],
                            row["completed_at"] or "",
                            row["error"] or "",
                        )
                    )

                output.flush()
                os.fsync(output.fileno())

            os.replace(temporary, destination)
        except Exception:
            temporary.unlink(missing_ok=True)
            raise

        return destination

    def integrity_check(self) -> None:
        with self._connect() as connection:
            row = connection.execute("PRAGMA integrity_check").fetchone()

        if row[0] != "ok":
            raise ManifestError(
                f"SQLite integrity check failed: {row[0]}"
            )

    @staticmethod
    def _validate_records(
        records: Sequence[tuple[str, str]],
    ) -> tuple[tuple[str, str], ...]:
        normalized: list[tuple[str, str]] = []
        seen: set[str] = set()

        for accession, taxid in records:
            accession = accession.strip()

            if not accession:
                raise ValueError("Accession must not be empty")

            if accession in seen:
                raise ValueError(
                    f"Duplicate accession within package: {accession}"
                )

            if isinstance(taxid, bool) or not isinstance(taxid, str):
                raise TypeError(
                    f"Taxid for {accession} must be an integer"
                )

            seen.add(accession)
            normalized.append((accession, taxid))

        return tuple(normalized)

    @staticmethod
    def _row_to_record(row: sqlite3.Row) -> ManifestRecord:
        return ManifestRecord(
            accession=row["accession"],
            taxid=row["taxid"],
            state=row["state"],
            package_id=row["package_id"],
            reserved_at=row["reserved_at"],
            completed_at=row["completed_at"],
            error=row["error"],
        )

    @staticmethod
    def _temporary_sibling(destination: Path) -> Path:
        return destination.with_name(
            f".{destination.name}.{uuid.uuid4().hex}.tmp"
        )
  