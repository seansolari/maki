import sqlite3
from typing import Dict, Iterable, List, Tuple

from .lists import _StringDataClass, AssemblyFileLists, BaktaFileLists
from .models import ASSEMBLY_BATCH_SCHEMA, BAKTA_BATCH_SCHEMA, BAKTA_SCHEMA, Table
from .sql import DbHandle, create_table, delete_table, insert_rows


def insert_file_list(conn: sqlite3.Connection, schema: Dict[str, List[str]], file_list: Iterable[_StringDataClass], table_def: Table):
    columns = table_def.columns()
    
    if table_def.name not in schema:
        create_table(conn, table_def.name, columns)
    
    batch: List[Tuple[str, ...]] = []
    batch_size = 10000

    for row in file_list:
        batch.append(row.astuple())

        if len(batch) >= batch_size:
            insert_rows(conn, table_def.name, columns, batch)
            batch = []

    if batch:
        insert_rows(conn, table_def.name, columns, batch)
        
    return table_def.name


def import_assembly_batches(db: DbHandle, schema: Dict[str, List[str]], force: bool = False):
    mani = AssemblyFileLists(db.path.parent / "atb.assembly.list.csv.gz")
    
    if ASSEMBLY_BATCH_SCHEMA.name in schema:
        if force:
            delete_table(db.conn, ASSEMBLY_BATCH_SCHEMA.name)
        else:
            return
    
    insert_file_list(db.conn, schema, mani.batch_files, ASSEMBLY_BATCH_SCHEMA)


def import_annotation_batches(db: DbHandle, schema: Dict[str, List[str]], force: bool = False):
    mani = BaktaFileLists(db.path.parent / "atb.bakta.list.csv.gz")
    
    if BAKTA_BATCH_SCHEMA.name in schema:
        if force:
            delete_table(db.conn, BAKTA_BATCH_SCHEMA.name)
        else:
            return _import_annotation_file_lists(db.conn, schema, mani, force=force)
    
    insert_file_list(db.conn, schema, mani.batch_files, BAKTA_BATCH_SCHEMA)
    _import_annotation_file_lists(db.conn, schema, mani, force=force)


def _import_annotation_file_lists(conn: sqlite3.Connection, schema: Dict[str, List[str]], mani: BaktaFileLists, force: bool = False):
    if BAKTA_SCHEMA.name in schema:
        if force:
            delete_table(conn, BAKTA_SCHEMA.name)
        else:
            return

    for lst in mani.lists():
        insert_file_list(conn, schema, lst.data, BAKTA_SCHEMA)
