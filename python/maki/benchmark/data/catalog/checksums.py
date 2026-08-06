# dataset_manager/checksums.py

from pathlib import Path
import hashlib


def file_checksum(path: str | Path, algorithm: str = "sha256") -> str:
    path = Path(path)
    h = hashlib.new(algorithm)

    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)

    return h.hexdigest()


def verify_checksum(path: str | Path, algorithm: str, expected: str) -> bool:
    actual = file_checksum(path, algorithm)
    return actual.lower() == expected.lower()