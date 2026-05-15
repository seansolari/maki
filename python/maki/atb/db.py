import sqlite3
import lzma
import shutil
import subprocess
from pathlib import Path
from typing import List


def connect(db_path: Path):
    if not db_path.exists():
        raise FileNotFoundError(f"Database not found: {db_path}")
    return sqlite3.connect(db_path)


def ensure_db(remote_src: str, db_path: Path):
    db_path = db_path.expanduser()
    db_path.parent.mkdir(parents=True, exist_ok=True)

    if db_path.exists():
        return db_path

    gz_path = db_path.with_suffix(".sqlite.xz")

    print("Downloading metadata database...")
    subprocess.check_call(["wget", "-O", gz_path, remote_src])

    print("Extracting...")
    with lzma.open(gz_path, "rb") as f_in:
        with open(db_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    gz_path.unlink()
    print("Done.")

    return db_path

  
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
