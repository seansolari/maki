from __future__ import annotations
import csv
import gzip
from dataclasses import astuple, dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable

from .osf_get_files_for_project import get_assembly_files, get_bakta_files, OSFFile


@dataclass(frozen=True)
class RemoteBatchFile:
    download_url: str
    tar_name: str
    tar_md5: str
    
    def __repr__(self) -> str:
        return f"RemoteBatchFile({self.tar_name}, {self.download_url}, {self.tar_md5})"


class _FileLists:
    def __init__(self, path: Path, getter: Callable[[], Iterable[OSFFile]], override: bool = False):
        if not path.exists() or override:
            self.download(path, getter)
            
        self._data: Dict[str, RemoteBatchFile] = {}
            
        with gzip.open(path, "rt") as f:
            reader = csv.DictReader(f)
            
            for r in reader:
                self._data[r["filename"]] = RemoteBatchFile(download_url=r["url"], tar_name=r["filename"], tar_md5=r["md5"])

    @staticmethod
    def download(path: Path, getter: Callable[[], Iterable[OSFFile]]):        
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.suffix != ".gz":
            path = path.with_suffix(".csv.gz")
            
        with gzip.open(path, "wt") as f:
            writer = csv.writer(f)
            writer.writerow(["project", "project_id", "filename", "url", "md5", "size"])
            for r in getter():
                writer.writerow(astuple(r))

    def __getitem__(self, key: str):
        return self._data.__getitem__(key)
      
    def __contains__(self, key: str):
        return self._data.__contains__(key)


class AssemblyFileLists(_FileLists):
    def __init__(self, path: Path, override: bool = False) -> None:
        super().__init__(path, get_assembly_files, override)


class BaktaFileLists(_FileLists):
    def __init__(self, path: Path, override: bool = False) -> None:
        super().__init__(path, get_bakta_files, override)
