import sqlite3
import lzma
import shutil
import subprocess
from pathlib import Path
from .config import SQLITE_URL


def connect(db_path: Path):
    if not db_path.exists():
        raise FileNotFoundError(f"Database not found: {db_path}")
    return sqlite3.connect(db_path)


def ensure_db(db_path: Path):
    db_path = db_path.expanduser()
    db_path.parent.mkdir(parents=True, exist_ok=True)

    if db_path.exists():
        return db_path

    gz_path = db_path.with_suffix(".sqlite.xz")

    print("Downloading metadata database...")
    subprocess.check_call(["wget", "-O", gz_path, SQLITE_URL])

    print("Extracting...")
    with lzma.open(gz_path, "rb") as f_in:
        with open(db_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    gz_path.unlink()
    print("Done.")

    return db_path
