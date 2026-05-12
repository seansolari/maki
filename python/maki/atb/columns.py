def resolve_columns(schema):
    sample_table = None

    for t, cols in schema.items():
        if "sample" in cols:
            sample_table = t
            break

    if not sample_table:
        raise RuntimeError("Could not find sample table")

    cols = schema[sample_table]

    def find(name):
        for c in cols:
            if name in c.lower():
                return c
        return None

    mapping = {
        "table": sample_table,
        "sample": "sample",
        "species": find("species"),
        "hq": find("hq"),
        "assembly": find("assembly_on_osf"),
        "lineage": find("lineage"),
    }

    # detect rank columns if present
    mapping["ranks"] = {
        r: r for r in ["phylum", "class", "order", "family", "genus"]
        if r in cols
    }

    return mapping
