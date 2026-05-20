
from pathlib import Path
from typing import Dict, Iterable, Tuple, TypeVar

from .base import BaseTaxonomy


def load_mapping(path: Path) -> Dict[str, str]:
    if not path.exists():
        return {}

    mapping: Dict[str, str] = {}
    with open(path) as f:
        for line in f:
            acc, taxid = line.strip().split("\t")
            mapping[acc] = taxid

    return mapping


def save_mapping(path: Path, mapping: Dict[str, str]):
    with open(path, "w") as f:
        for k, v in mapping.items():
            f.write(f"{k}\t{v}\n")


T = TypeVar("T")


def resolve_accession_taxids(records: Iterable[Tuple[str, str]], taxonomy: BaseTaxonomy, db_root: Path):
    mapping_file = db_root / "taxonomy" / "accession_to_taxid.tsv"

    existing = load_mapping(mapping_file)
    updated = dict(existing)

    unresolved = []

    for acc, taxid in records:
        if acc in existing:
            continue

        try:
            taxid = taxonomy.resolve_taxid(taxid)
            updated[acc] = taxid
        except Exception as e:
            print(f"[warning] {acc}: {e}")
            unresolved.append(acc)

    if unresolved:
        print(f"\n[warning] {len(unresolved)} genomes could not be resolved:")
        for u in unresolved[:10]:
            print(f" - {u}")

    save_mapping(mapping_file, updated)

    return updated
