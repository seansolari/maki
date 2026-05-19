
from itertools import repeat
from pathlib import Path
import shutil
import subprocess
from tempfile import TemporaryDirectory
from typing import List

from tqdm.contrib.concurrent import process_map


def _json_to_gff(json_path: str, output_base: Path, tmp: Path) -> None:
    gff_trg = (output_base / "gff" / json_path).with_suffix(".gff")
    
    if not gff_trg.exists():
        prefix = gff_trg.name.removesuffix(".gff")
        subprocess.check_call(["bakta_io", "--output", tmp / prefix, "--prefix", prefix, output_base / "bakta" / json_path],
                              stdout=subprocess.DEVNULL,
                              stderr=subprocess.DEVNULL)
        
        gff_trg.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(tmp / prefix / f"{prefix}.gff3", gff_trg)


def pjson_to_gff(json_files: List[str], output_base: Path, concurrency: int):
    if not shutil.which("bakta_io"):
        raise RuntimeError("Cannot find `bakta_io` command.")
    
    (output_base / "gff").mkdir(parents=True, exist_ok=True)
    
    with TemporaryDirectory(dir=output_base) as tmp:
        tmp = Path(tmp)
        process_map(_json_to_gff, json_files, repeat(output_base), repeat(tmp), max_workers=concurrency, chunksize=1)
