import sqlite3
from typing import Dict, List


def inspect_schema(conn: sqlite3.Connection) -> Dict[str, List[str]]:
    cur = conn.cursor()

    cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = [r[0] for r in cur.fetchall()]

    schema: Dict[str, List[str]] = {}

    for t in tables:
        cur.execute(f"PRAGMA table_info({t});")
        schema[t] = [row[1] for row in cur.fetchall()]

    return schema
