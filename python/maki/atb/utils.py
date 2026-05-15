import gzip
import re
from pathlib import Path


def infer_delimiter(file_path: Path):
    ptn = re.compile(r"([^\.]+)(?:\.gz)?$")
    sfx = ptn.findall(file_path.name)
    if not sfx or sfx[0] not in {"csv", "tsv"}:
      raise ValueError(f"Could not identify file type for table {file_path}")
    else:
      return "\t" if sfx[0] == "tsv" else ","


def open_maybe_gzip(file_path: Path):
    if file_path.suffix == ".gz":
        return gzip.open(file_path, "rt")
    return open(file_path, "r")
