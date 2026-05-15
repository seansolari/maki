import csv
from pathlib import Path
import sqlite3

from .db import create_table, insert_rows
from .lists import AssemblyFileLists
from .utils import infer_delimiter, open_maybe_gzip


def import_assembly_batch_files(conn: sqlite3.Connection):
    axl = AssemblyFileLists(Path.cwd() / "atb.assembly.list.csv.gz")
    
    # column definition and sanitize
    columns = 


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
  