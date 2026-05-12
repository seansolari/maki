def inspect_schema(conn):
    cur = conn.cursor()

    cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = [r[0] for r in cur.fetchall()]

    schema = {}

    for t in tables:
        cur.execute(f"PRAGMA table_info({t});")
        schema[t] = [row[1] for row in cur.fetchall()]

    return schema


def find_sample_table(schema):
    for table, cols in schema.items():
        if "sample" in cols:
            return table
    raise RuntimeError("No sample table found")
