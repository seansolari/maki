
"""Phylogeny generation helpers."""
from __future__ import annotations
from pathlib import Path
from typing import Iterable

def make_fork_newick(genome_names: Iterable[str], root_label: str = "d__Bacteria", branch_length: float = 1.0) -> str:
    """Construct a fork (star) phylogeny in Newick format."""
    names = list(genome_names)
    if not names:
        return f"{root_label};"
    bl = f":{branch_length}"
    leaves = ",".join(f"{n}{bl}" for n in names)
    return f"({leaves}){root_label};"

def write_newick(path: Path, newick: str) -> None:
    """Write `newick` to `path` as plain text."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write(newick + "\n")
