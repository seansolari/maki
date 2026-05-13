from .query_builder import QueryBuilder
from .models import CORE, STATUS, FILE_LIST


def apply_taxon_filter(qb: QueryBuilder, core, taxon, rank):
    if rank == "species":
        qb.where(f"{core.name}.{core.species_col} = ?", taxon)

    elif rank == "genus" and core.genus_col:
        qb.where(f"{core.name}.{core.genus_col} = ?", taxon)

    elif rank == "family" and core.family_col:
        qb.where(f"{core.name}.{core.family_col} = ?", taxon)

    elif core.lineage_col:
        prefix = rank[0] + "__"
        qb.where(f"{core.name}.{core.lineage_col} LIKE ?", f"%{prefix}{taxon}%")

    else:
        # fallback
        qb.where(f"{core.name}.{core.species_col} LIKE ?", f"%{taxon}%")

    return qb


def count(conn, taxon, rank, high_quality=False, assemblies=False):
    qb = QueryBuilder(CORE)

    qb = apply_taxon_filter(qb, CORE, taxon, rank)

    if high_quality or assemblies:
        qb.join(
            STATUS,
            f"{CORE.name}.{CORE.sample_col}",
            f"{STATUS.name}.{STATUS.sample_col}",
        )

        if high_quality:
            qb.where(f"{STATUS.name}.{STATUS.hq_col} = 'PASS'")

        if assemblies:
            qb.where(f"{STATUS.name}.{STATUS.assembly_col} = 1")

    query, params = qb.build("COUNT(*)")

    cur = conn.cursor()
    cur.execute(query, params)
    return cur.fetchone()[0]


def list_samples(conn, taxon, rank, high_quality=False, assemblies=False, limit=20):
    qb = QueryBuilder(CORE)

    qb = apply_taxon_filter(qb, CORE, taxon, rank)

    if high_quality or assemblies:
        qb.join(
            STATUS,
            f"{CORE.name}.{CORE.sample_col}",
            f"{STATUS.name}.{STATUS.sample_col}",
        )

        if high_quality:
            qb.where(f"{STATUS.name}.{STATUS.hq_col} = 'PASS'")

        if assemblies:
            qb.where(f"{STATUS.name}.{STATUS.assembly_col} = 1")

    query, params = qb.build(f"{CORE.name}.{CORE.sample_col}")

    query += " LIMIT ?"
    params.append(limit)

    cur = conn.cursor()
    cur.execute(query, params)
    return [r[0] for r in cur.fetchall()]


def plan_download(conn, taxon, rank, high_quality=True, assemblies=True):
    qb = QueryBuilder(CORE)

    qb = apply_taxon_filter(qb, CORE, taxon, rank)

    # join status
    qb.join(
        STATUS,
        f"{CORE.name}.{CORE.sample_col}",
        f"{STATUS.name}.{STATUS.sample_col}",
    )

    if high_quality:
        qb.where(f"{STATUS.name}.{STATUS.hq_col} = 'PASS'")

    if assemblies:
        qb.where(f"{STATUS.name}.{STATUS.assembly_col} = 1")

    # join file list
    qb.join(
        FILE_LIST,
        f"{CORE.name}.{CORE.sample_col}",
        f"{FILE_LIST.name}.{FILE_LIST.sample_col}",
    )

    query, params = qb.build(
        f"{FILE_LIST.name}.{FILE_LIST.archive_col}, COUNT(*)"
    )

    query += f"""
    GROUP BY {FILE_LIST.name}.{FILE_LIST.archive_col}
    ORDER BY COUNT(*) DESC
    """

    cur = conn.cursor()
    cur.execute(query, params)

    return cur.fetchall()
