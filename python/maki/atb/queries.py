import sqlite3

from .columns import SampleColumns


def build_taxon_filter(mapping: SampleColumns, taxon):
    return f"{mapping.species} = ?", [taxon]


def count_taxon(conn: sqlite3.Connection, mapping: SampleColumns, taxon: str) -> int:
    where, params = build_taxon_filter(mapping, taxon)

    query = f"""
    SELECT COUNT(*)
    FROM {mapping.table}
    WHERE {where}
    """

    cur = conn.cursor()
    cur.execute(query, params)
    return int(cur.fetchone()[0])


def count_all(conn: sqlite3.Connection, mapping: SampleColumns):
    query = f"""
    SELECT {mapping.species}, COUNT(*)
    FROM {mapping.table}
    GROUP BY {mapping.species}
    """

    cur = conn.cursor()
    cur.execute(query)
    return cur.fetchall()


def list_sample_rows(conn: sqlite3.Connection, mapping: SampleColumns, high_quality: bool, assemblies: bool, limit: int):
    conditions = []
    params = []

    if high_quality and mapping.hq:
        conditions.append(f"{mapping.hq} = 'PASS'")

    if assemblies and mapping.assembly:
        conditions.append(f"{mapping.assembly} = 1")

    query = f"""
    SELECT *
    FROM {mapping.table}
    {("WHERE " + " AND ".join(conditions)) if conditions else ""}
    LIMIT ?
    """

    params.append(limit)

    cur = conn.cursor()
    cur.execute(query, params)

    return cur.fetchall()
  