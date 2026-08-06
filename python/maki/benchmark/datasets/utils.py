from __future__ import annotations

from pathlib import Path
import hashlib
import shutil


def compute_checksum(
    path: str | Path,
    algorithm: str = "sha256",
) -> str:
    path = Path(path)
    hasher = hashlib.new(algorithm)

    if path.is_dir():
        for child in sorted(path.rglob("*")):
            if child.is_file():
                hasher.update(str(child.relative_to(path)).encode())
                with child.open("rb") as handle:
                    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                        hasher.update(chunk)
    else:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                hasher.update(chunk)

    return hasher.hexdigest()


def copy_path(source: str | Path, destination: str | Path) -> None:
    source = Path(source)
    destination = Path(destination)

    if destination.exists():
        if destination.is_dir():
            shutil.rmtree(destination)
        else:
            destination.unlink()

    if source.is_dir():
        shutil.copytree(source, destination)
    else:
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
