import gzip
import hashlib
import io
import re
from pathlib import Path
from typing import Literal, overload


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


def md5sum(path: Path, chunk_size: int = 8192):
    h = hashlib.md5()
    with open(path, "rb") as f:
        while chunk := f.read(chunk_size):
            h.update(chunk)
    return h.hexdigest()
