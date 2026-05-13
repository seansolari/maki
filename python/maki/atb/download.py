
def plan_download(conn, schema, mapping, taxon, high_quality, assemblies):
    file_table = find_file_table(schema)

    if not file_table:
        raise RuntimeError("No file list table found in DB")

    conditions = []
    params = []

    from .queries import build_taxon_filter

    where, p = build_taxon_filter(mapping, taxon)
    conditions.append(where)
    params.extend(p)

    if high_quality and mapping["hq"]:
        conditions.append(f"s.{mapping['hq']} = 'PASS'")

    if assemblies and mapping["assembly"]:
        conditions.append(f"s.{mapping['assembly']} = 1")

    query = f"""
    SELECT f.archive_path, COUNT(*)
    FROM {mapping['table']} s
    JOIN {file_table} f ON s.sample = f.sample
    WHERE {" AND ".join(conditions)}
    GROUP BY f.archive_path
    ORDER BY COUNT(*) DESC
    """

    cur = conn.cursor()
    cur.execute(query, params)

    return cur.fetchall()
  