def build_taxon_filter(mapping, taxon, rank):
    if rank == "species" and mapping["species"]:
        return f"{mapping['species']} = ?", [taxon]

    if rank in mapping["ranks"]:
        return f"{mapping['ranks'][rank]} = ?", [taxon]

    if mapping["lineage"]:
        prefix = rank[0] + "__"
        return f"{mapping['lineage']} LIKE ?", [f"%{prefix}{taxon}%"]

    # fallback
    return f"{mapping['species']} LIKE ?", [f"%{taxon}%"]


def count(conn, mapping, taxon, rank):
    where, params = build_taxon_filter(mapping, taxon, rank)

    query = f"""
    SELECT COUNT(*) FROM {mapping['table']}
    WHERE {where}
    """

    cur = conn.cursor()
    cur.execute(query, params)
    return cur.fetchone()[0]


def list_samples(conn, mapping, taxon, rank, high_quality, assemblies, limit):
    conditions = []
    params = []

    where, p = build_taxon_filter(mapping, taxon, rank)
    conditions.append(where)
    params.extend(p)

    if high_quality and mapping["hq"]:
        conditions.append(f"{mapping['hq']} = 'PASS'")

    if assemblies and mapping["assembly"]:
        conditions.append(f"{mapping['assembly']} = 1")

    query = f"""
    SELECT {mapping['sample']}
    FROM {mapping['table']}
    WHERE {" AND ".join(conditions)}
    LIMIT ?
    """

    params.append(limit)

    cur = conn.cursor()
    cur.execute(query, params)

    return [r[0] for r in cur.fetchall()]
  