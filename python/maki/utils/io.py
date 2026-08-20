import gzip
import hashlib
import io
import re
from pathlib import Path
from typing import Iterator, Literal, Tuple, overload


def infer_delimiter(file_path: Path):
    ptn = re.compile(r"([^\.]+)(?:\.gz)?$")
    sfx = ptn.findall(file_path.name)
    if not sfx or sfx[0] not in {"csv", "tsv"}:
      raise ValueError(f"Could not identify file type for table {file_path}")
    else:
      return "\t" if sfx[0] == "tsv" else ","


@overload
def open_maybe_gzip(path, mode: Literal["wt"]) -> io.TextIOWrapper:
    ...
    
@overload
def open_maybe_gzip(path, mode: Literal["rt"]) -> io.TextIOWrapper:
    ...

def open_maybe_gzip(path, mode="rt"):
    """
    Open plain text or gzip-compressed files based on filename suffix.
    """
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def read_fasta(path: str | Path):
    """Minimal FASTA reader to avoid adding screed or Biopython as a direct
    dependency of this module.

    sourmash itself may install screed, but this keeps the package boundary
    small and explicit.
    """
    with open_maybe_gzip(path, "rt") as fp:
        yield from _yield_fasta_data(fp)


def read_fasta_from_gff(path: str | Path):
    with open_maybe_gzip(path, "rt") as fp:
        for line in fp:
            if line.startswith("##FASTA"):
                break

        yield from _yield_fasta_data(fp)


def _yield_fasta_data(src: io.TextIOWrapper) -> Iterator[Tuple[str, str]]:
    header: str | None = None
    sequence_parts: list[str] = []

    for raw_line in src:
        line = raw_line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if header and sequence_parts:
                yield header, "".join(sequence_parts)
                sequence_parts.clear()
            
            header = line[1:].split()[0].strip()
        else:
            sequence_parts.append(line)

    if header and sequence_parts:
        yield header, "".join(sequence_parts)


def md5sum(path: Path, chunk_size: int = 8192):
    h = hashlib.md5()
    with open(path, "rb") as f:
        while chunk := f.read(chunk_size):
            h.update(chunk)
    return h.hexdigest()
