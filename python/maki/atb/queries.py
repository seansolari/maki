from dataclasses import dataclass
import sqlite3
from typing import Optional

from .query_builder import QueryBuilder
from .models import ASSEMBLY_SCHEMA, ASSEMBLY_STATS_SCHEMA, CHECKM2_SCHEMA, BAKTA_SCHEMA


def join_table(qb: QueryBuilder, tbl):
    return qb.join(
        tbl,
        f"{ASSEMBLY_SCHEMA.name}.{ASSEMBLY_SCHEMA.sample}",
        f"{tbl.name}.{tbl.sample}",
    )


@dataclass
class QueryOptions:
    taxon: Optional[str] = None
    high_quality: bool = False
    has_assembly: bool = False
    max_contamination: Optional[float] = None
    min_completeness: Optional[float] = None
    has_annotation: bool = False
    
    def enforce(self, qb: QueryBuilder):
        if self.taxon:
            qb.where(f"{ASSEMBLY_SCHEMA.name}.{ASSEMBLY_SCHEMA.species} = ?", self.taxon)
            
        if self.high_quality:
            qb.where(f"{ASSEMBLY_SCHEMA.name}.{ASSEMBLY_SCHEMA.hq} = 'PASS'")
          
        if self.has_assembly:
            qb.where(f"{ASSEMBLY_SCHEMA.name}.{ASSEMBLY_SCHEMA.assembly_exists} = 1")
        
        if self.max_contamination or self.min_completeness:
            join_table(qb, CHECKM2_SCHEMA)
            
            if self.max_contamination:
                qb.where(f"{CHECKM2_SCHEMA.name}.{CHECKM2_SCHEMA.contamination} <= ?", self.max_contamination)
                
            if self.min_completeness:
                qb.where(f"{CHECKM2_SCHEMA.name}.{CHECKM2_SCHEMA.completeness} >= ?", self.min_completeness)
        
        if self.has_annotation:
            join_table(qb, BAKTA_SCHEMA)
            
        return qb


def count_taxon(conn: sqlite3.Connection, opts: QueryOptions):
    qb = QueryBuilder(ASSEMBLY_SCHEMA)
    
    opts.enforce(qb)
    
    query, params = qb.build("COUNT(*)")

    cur = conn.cursor()
    cur.execute(query, params)
    return cur.fetchone()[0]


def count_all(conn: sqlite3.Connection, opts: QueryOptions):
    qb = QueryBuilder(ASSEMBLY_SCHEMA)
    
    opts.enforce(qb)

    qb.group_by(ASSEMBLY_SCHEMA.species)
    qb.order_by("COUNT(*) DESC")

    query, params = qb.build(f"{ASSEMBLY_SCHEMA.species}, COUNT(*)")
    
    cur = conn.cursor()
    cur.execute(query, params)
    return cur.fetchall()


def list_sample_rows(conn: sqlite3.Connection, opts: QueryOptions, limit: int = 20):
    qb = QueryBuilder(ASSEMBLY_SCHEMA)

    opts.enforce(qb)
    qb.set_limit(limit)

    query, params = qb.build(f"*")

    cur = conn.cursor()
    cur.execute(query, params)
    return cur.fetchall()


# def plan_download(conn: sqlite3.Connection, opts: QueryOptions):
#     qb = QueryBuilder(ASSEMBLY_SCHEMA)
# 
#     opts.enforce(qb)
# 
#     query, params = qb.build(
#         f"{FILE_LIST.name}.{FILE_LIST.archive_col}, COUNT(*)"
#     )
# 
#     query += f"""
#     GROUP BY {FILE_LIST.name}.{FILE_LIST.archive_col}
#     ORDER BY COUNT(*) DESC
#     """
# 
#     cur = conn.cursor()
#     cur.execute(query, params)
# 
#     return cur.fetchall()
