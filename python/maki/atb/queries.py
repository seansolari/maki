from dataclasses import dataclass
import sqlite3
from typing import Optional

from .query_builder import QueryBuilder
from .models import ASSEMBLY_BATCH_SCHEMA, ASSEMBLY_SCHEMA, ASSEMBLY_STATS_SCHEMA, BAKTA_BATCH_SCHEMA, CHECKM2_SCHEMA, BAKTA_SCHEMA


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
            qb.join(CHECKM2_SCHEMA,
                    f"{ASSEMBLY_SCHEMA.name}.{ASSEMBLY_SCHEMA.sample}",
                    f"{CHECKM2_SCHEMA.name}.{CHECKM2_SCHEMA.sample}")
            
            if self.max_contamination:
                qb.where(f"{CHECKM2_SCHEMA.name}.{CHECKM2_SCHEMA.contamination} <= ?", self.max_contamination)
                
            if self.min_completeness:
                qb.where(f"{CHECKM2_SCHEMA.name}.{CHECKM2_SCHEMA.completeness} >= ?", self.min_completeness)
        
        if self.has_annotation:
            qb.join(BAKTA_SCHEMA,
                    f"{ASSEMBLY_SCHEMA.name}.{ASSEMBLY_SCHEMA.sample}",
                    f"{BAKTA_SCHEMA.name}.{BAKTA_SCHEMA.sample}")
          
            qb.where(f"{BAKTA_SCHEMA.name}.{BAKTA_SCHEMA.status} = 'PASS'")
            
        return qb


def count_taxon(conn: sqlite3.Connection, opts: QueryOptions):
    assert opts.taxon

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

    columns = qb.colnames()
    query, params = qb.build(qb.autocols())

    cur = conn.cursor()
    cur.execute(query, params)
    return columns, cur.fetchall()


def plan_download(conn: sqlite3.Connection, opts: QueryOptions):
    qb = QueryBuilder(ASSEMBLY_SCHEMA)

    opts.enforce(qb)
    qb.join(ASSEMBLY_STATS_SCHEMA,
            f"{ASSEMBLY_SCHEMA.name}.{ASSEMBLY_SCHEMA.sample}",
            f"{ASSEMBLY_STATS_SCHEMA.name}.{ASSEMBLY_STATS_SCHEMA.sample}")
    qb.join(ASSEMBLY_BATCH_SCHEMA,
            f"{ASSEMBLY_SCHEMA.name}.{ASSEMBLY_SCHEMA.tar_xz}",
            f"{ASSEMBLY_BATCH_SCHEMA.name}.{ASSEMBLY_BATCH_SCHEMA.tar_xz}")
    if opts.has_annotation:
        qb.join(BAKTA_BATCH_SCHEMA,
                f"{BAKTA_SCHEMA.name}.{BAKTA_SCHEMA.tar_xz}",
                f"{BAKTA_BATCH_SCHEMA.name}.{BAKTA_BATCH_SCHEMA.tar_xz}")

    columns = qb.colnames()
    query, params = qb.build(qb.autocols())

    cur = conn.cursor()
    cur.execute(query, params)
    return columns, cur.fetchall()
