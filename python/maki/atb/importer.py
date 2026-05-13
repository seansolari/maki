import re
import gzip
import csv
from pathlib import Path
import sqlite3
from typing import List


def infer_delimiter(file_path: Path):
    ptn = re.compile(r"([^\.]+)(?:\.gz)?$")
    sfx = ptn.findall(file_path.name)
    if not sfx or sfx[0] not in {"csv", "tsv"}:
      raise ValueError(f"Could not identify file type for table {file_path}")
    else:
      return "\t" if sfx[0] == "tsv" else ","


def open_maybe_gzip(file_path: Path):
    if file_path.suffix == ".gz":
        return gzip.open(file_path, "rt")
    return open(file_path, "r")


def create_table(conn: sqlite3.Connection, table_name: str, columns: List[str], unique_col: str = "sample"):
    cur = conn.cursor()

    col_defs = []

    for col in columns:
        if col == unique_col:
            col_defs.append(f"{col} TEXT PRIMARY KEY")
        else:
            col_defs.append(f"{col} TEXT")

    schema = ", ".join(col_defs)

    query = f"CREATE TABLE IF NOT EXISTS {table_name} ({schema})"
    cur.execute(query)


def insert_rows(conn: sqlite3.Connection, table_name: str, columns: List[str], rows: List[List[str]]):
    cur = conn.cursor()

    placeholders = ",".join(["?"] * len(columns))
    col_str = ",".join(columns)

    query = f"""
    INSERT OR IGNORE INTO {table_name} ({col_str})
    VALUES ({placeholders})
    """

    cur.executemany(query, rows)
    conn.commit()


def import_csv_to_sqlite(conn: sqlite3.Connection, file_path: Path, table_name: str):
    delimiter = infer_delimiter(file_path)

    with open_maybe_gzip(file_path) as f:
        reader = csv.reader(f, delimiter=delimiter)

        # header
        columns = next(reader)

        # sanitize column names
        columns = [c.strip().replace(" ", "_") for c in columns]

        # ensure "sample" exists for deduplication
        if "sample" not in columns:
            raise ValueError("Input file must contain a 'sample' column")

        create_table(conn, table_name, columns)

        batch = []
        batch_size = 10000

        for row in reader:
            batch.append(row)

            if len(batch) >= batch_size:
                insert_rows(conn, table_name, columns, batch)
                batch = []

        if batch:
            insert_rows(conn, table_name, columns, batch)

    return table_name
  