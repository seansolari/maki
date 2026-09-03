from __future__ import annotations

import os
import sqlite3
from collections import defaultdict
from dataclasses import dataclass
from itertools import islice
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence, TypeVar


@dataclass(frozen=True, slots=True)
class AnnotationRecord:
    accession: str
    seed_name: str
    database_name: str | None
    label: str | None


@dataclass(frozen=True, slots=True)
class SeedTaxonomy:
    accession: str
    taxid: str


@dataclass(frozen=True, slots=True)
class SeedInfo:
    seed_name: str
    accession: str
    taxid: str | None


T = TypeVar("T")

def _batched(iterable: Iterable[T], size: int) -> Iterator[list[T]]:
    """
    Yield lists containing at most `size` items.

    This is provided instead of itertools.batched so that the module works
    with Python versions earlier than 3.12.
    """
    if size < 1:
        raise ValueError("Batch size must be at least 1")

    iterator = iter(iterable)

    while batch := list(islice(iterator, size)):
        yield batch


def _chunks(values: Sequence, size: int = 900) -> Iterator[Sequence]:
    """
    Split values into chunks small enough for SQLite parameter limits.
    """
    for start in range(0, len(values), size):
        yield values[start:start + size]


def _quote_identifier(identifier: str) -> str:
    """
    Quote an SQLite identifier.

    Database names supplied by users are never used as table names directly,
    but quoting generated identifiers makes the boundary explicit.
    """
    return '"' + identifier.replace('"', '""') + '"'


class AnnotationIndexBuilder:
    """
    Build or append to an SQLite annotation index.

    The builder owns a read-write SQLite connection and is not thread-safe.
    Distinct builder instances should not write to the same database
    concurrently.

    Parameters
    ----------
    path:
        SQLite database path.
    """

    SCHEMA_VERSION = 1

    def __init__(self, path: str | os.PathLike[str]) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)

        self._connection = sqlite3.connect(self.path)
        self._connection.execute("PRAGMA foreign_keys = ON")
        self._connection.execute("PRAGMA journal_mode = WAL")
        self._connection.execute("PRAGMA synchronous = NORMAL")
        self._connection.execute("PRAGMA temp_store = MEMORY")

        self._database_table_cache: dict[str, str] = {}

        self._create_schema()
        self._check_schema_version()
        self._assert_not_finalized()

    def __enter__(self) -> AnnotationIndexBuilder:
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()

    def close(self) -> None:
        if self._connection is not None:
            self._connection.close()
            self._connection = None  # type: ignore[assignment]

    def _create_schema(self) -> None:
        assert self._connection
        
        with self._connection:
            self._connection.executescript(
                """
                CREATE TABLE IF NOT EXISTS metadata (
                    key   TEXT PRIMARY KEY,
                    value TEXT NOT NULL
                ) WITHOUT ROWID;

                CREATE TABLE IF NOT EXISTS seed (
                    seed_id   INTEGER PRIMARY KEY,
                    seed_name TEXT NOT NULL UNIQUE,
                    accession TEXT NOT NULL
                );

                CREATE INDEX IF NOT EXISTS seed_accession_idx
                    ON seed(accession);

                CREATE TABLE IF NOT EXISTS annotation_database (
                    database_id INTEGER PRIMARY KEY,
                    name        TEXT NOT NULL UNIQUE,
                    table_name  TEXT NOT NULL UNIQUE
                );

                CREATE TABLE IF NOT EXISTS accession_taxonomy (
                    accession TEXT PRIMARY KEY,
                    taxid     TEXT NOT NULL
                ) WITHOUT ROWID;
                """
            )

            self._connection.execute(
                """
                INSERT OR IGNORE INTO metadata(key, value)
                VALUES ('schema_version', ?)
                """,
                (str(self.SCHEMA_VERSION),),
            )

            self._connection.execute(
                """
                INSERT OR IGNORE INTO metadata(key, value)
                VALUES ('finalized', '0')
                """
            )

    def _check_schema_version(self) -> None:
        assert self._connection
        
        row = self._connection.execute(
            """
            SELECT value
            FROM metadata
            WHERE key = 'schema_version'
            """
        ).fetchone()

        if row is None:
            raise RuntimeError("Database has no schema version")

        version = int(row[0])

        if version != self.SCHEMA_VERSION:
            raise RuntimeError(
                f"Unsupported annotation index schema version {version}; "
                f"expected {self.SCHEMA_VERSION}"
            )

    def _assert_not_finalized(self) -> None:
        assert self._connection
        
        row = self._connection.execute(
            """
            SELECT value
            FROM metadata
            WHERE key = 'finalized'
            """
        ).fetchone()

        if row is not None and row[0] == "1":
            raise RuntimeError(
                f"Annotation index has already been finalized: {self.path}"
            )

    def _get_or_create_database_table(
        self,
        database_name: str,
    ) -> str:
        """
        Return the internal table name for an annotation database.

        External names are represented by generated table names such as
        annotation_1. This means names containing spaces, punctuation, or SQL
        keywords remain safe.
        """
        assert self._connection
        
        if not database_name:
            raise ValueError("database_name cannot be empty")

        cached = self._database_table_cache.get(database_name)
        if cached is not None:
            return cached

        row = self._connection.execute(
            """
            SELECT table_name
            FROM annotation_database
            WHERE name = ?
            """,
            (database_name,),
        ).fetchone()

        if row is not None:
            table_name = row[0]
            self._database_table_cache[database_name] = table_name
            return table_name

        cursor = self._connection.execute(
            """
            INSERT INTO annotation_database(name, table_name)
            VALUES (?, '')
            """,
            (database_name,),
        )
        database_id = cursor.lastrowid

        if database_id is None:
            raise RuntimeError(
                f"Could not allocate database ID for {database_name!r}"
            )

        table_name = f"annotation_{database_id}"
        quoted_table_name = _quote_identifier(table_name)

        self._connection.execute(
            """
            UPDATE annotation_database
            SET table_name = ?
            WHERE database_id = ?
            """,
            (table_name, database_id),
        )

        self._connection.execute(
            f"""
            CREATE TABLE {quoted_table_name} (
                seed_id INTEGER NOT NULL,
                label   TEXT NOT NULL,

                PRIMARY KEY (seed_id, label),
                FOREIGN KEY (seed_id)
                    REFERENCES seed(seed_id)
                    ON DELETE CASCADE
            ) WITHOUT ROWID
            """
        )

        self._database_table_cache[database_name] = table_name
        return table_name

    def add_annotations(
        self,
        records: Iterable[AnnotationRecord],
        *,
        batch_size: int = 10_000,
    ) -> int:
        """
        Insert AnnotationRecord objects into the index.

        Existing seed-label pairs are ignored, making repeated insertion
        idempotent. Each batch is committed independently so that very large
        iterables do not need to be materialized in memory.

        A seed name already associated with a different accession causes a
        ValueError rather than silently changing the association.

        Returns
        -------
        int
            Number of input records processed. This includes records whose
            annotation pair was already present.
        """
        assert self._connection
        
        self._assert_not_finalized()

        processed = 0

        for batch in _batched(records, batch_size):
            self._validate_records(batch)

            with self._connection:
                seed_accessions = {
                    record.seed_name: record.accession
                    for record in batch
                }

                self._validate_batch_seed_accessions(batch)

                self._connection.executemany(
                    """
                    INSERT OR IGNORE INTO seed(seed_name, accession)
                    VALUES (?, ?)
                    """,
                    seed_accessions.items(),
                )

                self._validate_stored_seed_accessions(seed_accessions)

                seed_ids = self._lookup_seed_ids(
                    list(seed_accessions)
                )

                annotations_by_database: dict[
                    str,
                    set[tuple[int, str]],
                ] = defaultdict(set)

                for record in batch:
                    if record.database_name and record.label:
                        annotations_by_database[record.database_name].add(
                            (seed_ids[record.seed_name], record.label)
                        )

                for database_name, annotation_rows in (
                    annotations_by_database.items()
                ):
                    table_name = self._get_or_create_database_table(
                        database_name
                    )
                    quoted_table_name = _quote_identifier(table_name)

                    self._connection.executemany(
                        f"""
                        INSERT OR IGNORE INTO {quoted_table_name}(
                            seed_id,
                            label
                        )
                        VALUES (?, ?)
                        """,
                        annotation_rows,
                    )

            processed += len(batch)

        return processed

    @staticmethod
    def _validate_records(
        records: Sequence[AnnotationRecord],
    ) -> None:
        for record in records:
            if not record.accession:
                raise ValueError("accession cannot be empty")
            if not record.seed_name:
                raise ValueError("seed_name cannot be empty")

    @staticmethod
    def _validate_batch_seed_accessions(
        records: Sequence[AnnotationRecord],
    ) -> None:
        seed_accessions: dict[str, str] = {}

        for record in records:
            previous = seed_accessions.setdefault(
                record.seed_name,
                record.accession,
            )

            if previous != record.accession:
                raise ValueError(
                    f"Seed {record.seed_name!r} occurs with multiple "
                    f"accessions in the same batch: "
                    f"{previous!r} and {record.accession!r}"
                )

    def _validate_stored_seed_accessions(
        self,
        expected: Mapping[str, str],
    ) -> None:
        stored = self._lookup_seed_rows(list(expected))

        for seed_name, accession in stored:
            expected_accession = expected[seed_name]

            if accession != expected_accession:
                raise ValueError(
                    f"Seed {seed_name!r} is already associated with "
                    f"accession {accession!r}, not "
                    f"{expected_accession!r}"
                )

    def _lookup_seed_rows(
        self,
        seed_names: Sequence[str],
    ) -> list[tuple[str, str]]:
        assert self._connection
        
        rows: list[tuple[str, str]] = []

        for chunk in _chunks(seed_names):
            placeholders = ",".join("?" for _ in chunk)
            rows.extend(
                self._connection.execute(
                    f"""
                    SELECT seed_name, accession
                    FROM seed
                    WHERE seed_name IN ({placeholders})
                    """,
                    chunk,
                )
            )

        return rows

    def _lookup_seed_ids(
        self,
        seed_names: Sequence[str],
    ) -> dict[str, int]:
        assert self._connection
        
        result: dict[str, int] = {}

        for chunk in _chunks(seed_names):
            placeholders = ",".join("?" for _ in chunk)

            for seed_id, seed_name in self._connection.execute(
                f"""
                SELECT seed_id, seed_name
                FROM seed
                WHERE seed_name IN ({placeholders})
                """,
                chunk,
            ):
                result[seed_name] = seed_id

        missing = set(seed_names).difference(result)
        if missing:
            raise RuntimeError(
                f"Failed to resolve inserted seeds: {sorted(missing)!r}"
            )

        return result

    def update_taxonomy(
        self,
        items: Iterable[SeedTaxonomy],
        *,
        batch_size: int = 10_000,
    ) -> int:
        """
        Insert or update accession-to-taxonomy mappings.

        Parameters
        ----------
        taxonomy:
            One of:

            * a mapping of accession to taxonomy ID;
            * an iterable of SeedTaxonomy objects; or
            * an iterable of ``(accession, taxid)`` tuples.

        batch_size:
            Number of mappings committed per transaction.

        Returns
        -------
        int
            Number of input mappings processed.
        """
        assert self._connection
        
        self._assert_not_finalized()

        processed = 0

        for batch in _batched(items, batch_size):
            normalized = [
                (str(item.accession), str(item.taxid))
                for item in batch
            ]

            for accession, taxid in normalized:
                if not accession:
                    raise ValueError("accession cannot be empty")
                if not taxid:
                    raise ValueError("taxid cannot be empty")

            with self._connection:
                self._connection.executemany(
                    """
                    INSERT INTO accession_taxonomy(accession, taxid)
                    VALUES (?, ?)
                    ON CONFLICT(accession) DO UPDATE SET
                        taxid = excluded.taxid
                    """,
                    normalized,
                )

            processed += len(normalized)

        return processed

    def finalize(
        self,
        *,
        compact: bool = True,
        make_file_read_only: bool = False,
    ) -> None:
        """
        Finalize the index for read-only use.

        Finalization:

        1. marks the database as finalized;
        2. updates SQLite query-planner statistics;
        3. checkpoints and removes the WAL journal;
        4. optionally compacts the database with VACUUM; and
        5. optionally removes filesystem write permissions.

        Once finalized, AnnotationIndexBuilder will reject subsequent writes.

        Setting ``make_file_read_only=True`` changes filesystem permissions.
        This is useful for a distributed static index, but is optional because
        permission semantics vary between operating systems and deployment
        environments.
        """
        assert self._connection
        
        self._assert_not_finalized()

        with self._connection:
            self._connection.execute(
                """
                UPDATE metadata
                SET value = '1'
                WHERE key = 'finalized'
                """
            )
            self._connection.execute("ANALYZE")
            self._connection.execute("PRAGMA optimize")

        if compact:
            self._connection.execute("VACUUM")

        # Move out of WAL mode so the finalized index is a single file.
        self._connection.execute("PRAGMA wal_checkpoint(TRUNCATE)")
        self._connection.execute("PRAGMA journal_mode = DELETE")

        self.close()

        if make_file_read_only:
            current_mode = self.path.stat().st_mode
            self.path.chmod(current_mode & ~0o222)


class AnnotationIndex:
    """
    Read-only interface to a finalized annotation index.

    The SQLite URI ``mode=ro`` prevents accidental modification even if the
    underlying file remains writable at the filesystem level.
    """

    SUPPORTED_SCHEMA_VERSION = AnnotationIndexBuilder.SCHEMA_VERSION

    def __init__(self, path: str | os.PathLike[str]) -> None:
        self.path = Path(path)

        if not self.path.is_file():
            raise FileNotFoundError(self.path)

        uri = f"{self.path.resolve().as_uri()}?mode=ro"
        self._connection = sqlite3.connect(uri, uri=True)
        self._connection.execute("PRAGMA query_only = ON")
        self._connection.execute("PRAGMA foreign_keys = ON")

        self._database_tables: dict[str, str] | None = None

        self._check_schema()

    def __enter__(self) -> AnnotationIndex:
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()

    def close(self) -> None:
        if self._connection is not None:
            self._connection.close()
            self._connection = None  # type: ignore[assignment]

    def _check_schema(self) -> None:
        assert self._connection
        
        metadata = dict(
            self._connection.execute(
                """
                SELECT key, value
                FROM metadata
                WHERE key IN ('schema_version', 'finalized')
                """
            )
        )

        version = int(metadata.get("schema_version", -1))

        if version != self.SUPPORTED_SCHEMA_VERSION:
            raise RuntimeError(
                f"Unsupported annotation index schema version {version}; "
                f"expected {self.SUPPORTED_SCHEMA_VERSION}"
            )

        if metadata.get("finalized") != "1":
            raise RuntimeError(
                "The annotation index has not been finalized"
            )

    def _load_database_tables(self) -> dict[str, str]:
        assert self._connection
        
        if self._database_tables is None:
            self._database_tables = dict(
                self._connection.execute(
                    """
                    SELECT name, table_name
                    FROM annotation_database
                    ORDER BY name
                    """
                )
            )

        return self._database_tables

    def list_databases(self) -> list:
        """
        Return external annotation database names in alphabetical order.
        """
        return list(self._load_database_tables())

    def _resolve_seed_ids(
        self,
        seed_names: Sequence[str],
    ) -> dict[str, int]:
        assert self._connection
        
        seed_names = list(dict.fromkeys(seed_names))
        result: dict[str, int] = {}

        for chunk in _chunks(seed_names):
            placeholders = ",".join("?" for _ in chunk)

            for seed_id, seed_name in self._connection.execute(
                f"""
                SELECT seed_id, seed_name
                FROM seed
                WHERE seed_name IN ({placeholders})
                """,
                chunk,
            ):
                result[seed_name] = seed_id

        return result

    def get_annotations(
        self,
        seed_names: Iterable[str],
        *,
        database_name: str | None = None,
    ) -> dict[str, dict[str, list[str]]]:
        """
        Retrieve annotations for a collection of seed names.

        Parameters
        ----------
        seed_names:
            Seed names to query.
        database_name:
            Restrict the query to one annotation database. If omitted,
            annotations from every registered database are returned.

        Returns
        -------
        dict
            Nested mapping in the form::

                {
                    seed_name: {
                        database_name: [label, ...]
                    }
                }

            Requested seeds with no matching annotations are returned with an
            empty dictionary. Unknown seed names are also represented by an
            empty dictionary.
        """
        assert self._connection
        
        requested = list(dict.fromkeys(map(str, seed_names)))
        result: dict[str, dict[str, list[str]]] = {
            seed_name: {}
            for seed_name in requested
        }

        if not requested:
            return result

        seed_ids = self._resolve_seed_ids(requested)
        if not seed_ids:
            return result

        seed_names_by_id = {
            seed_id: seed_name
            for seed_name, seed_id in seed_ids.items()
        }

        database_tables = self._load_database_tables()

        if database_name is None:
            selected_databases = database_tables.items()
        else:
            try:
                table_name = database_tables[database_name]
            except KeyError as error:
                raise KeyError(
                    f"Unknown annotation database {database_name!r}"
                ) from error

            selected_databases = [(database_name, table_name)]

        seed_id_values = list(seed_names_by_id)

        for current_database, table_name in selected_databases:
            quoted_table_name = _quote_identifier(table_name)

            for chunk in _chunks(seed_id_values):
                placeholders = ",".join("?" for _ in chunk)

                for seed_id, label in self._connection.execute(
                    f"""
                    SELECT seed_id, label
                    FROM {quoted_table_name}
                    WHERE seed_id IN ({placeholders})
                    ORDER BY seed_id, label
                    """,
                    chunk,
                ):
                    seed_name = seed_names_by_id[seed_id]

                    result[seed_name].setdefault(
                        current_database,
                        [],
                    ).append(label)

        return result

    def get_seed_info(
        self,
        seed_names: Iterable[str],
    ) -> dict[str, SeedInfo | None]:
        """
        Retrieve accession and taxonomy information for seed names.

        Unknown seeds map to ``None``. Known seeds without an accession
        taxonomy mapping have ``taxid=None``.

        Returns
        -------
        dict
            Mapping from each requested seed name to SeedInfo or None.
        """
        assert self._connection
        
        requested = list(dict.fromkeys(map(str, seed_names)))
        result: dict[str, SeedInfo | None] = {
            seed_name: None
            for seed_name in requested
        }

        for chunk in _chunks(requested):
            placeholders = ",".join("?" for _ in chunk)

            rows = self._connection.execute(
                f"""
                SELECT
                    s.seed_name,
                    s.accession,
                    t.taxid
                FROM seed AS s
                LEFT JOIN accession_taxonomy AS t
                    ON t.accession = s.accession
                WHERE s.seed_name IN ({placeholders})
                """,
                chunk,
            )

            for seed_name, accession, taxid in rows:
                result[seed_name] = SeedInfo(
                    seed_name=seed_name,
                    accession=accession,
                    taxid=taxid,
                )

        return result
