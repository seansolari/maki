from __future__ import annotations
import csv
import gzip
from dataclasses import astuple, dataclass
import os
from pathlib import Path
import requests
from typing import List, Tuple

from .utils import infer_delimiter
from .extern.osf_get_files_for_project import to_table, get_all_data


@dataclass(frozen=True)
class _StringDataClass:
    def astuple(self) -> Tuple[str, ...]:
      return astuple(self)


@dataclass(frozen=True)
class RemoteBatchFile(_StringDataClass):
    tar_xz: str
    tar_xz_url: str
    tar_xz_md5: str
    tar_xz_size_MB: str
    
    def __repr__(self) -> str:
        return f"RemoteBatchFile({self.tar_xz}, {self.tar_xz_url}, {self.tar_xz_md5}, {self.tar_xz_size_MB})"


class _FileLists:
    def __init__(self, path: Path, project_id: str, override: bool = False):
        self.path = path
      
        if not path.exists() or override:
            self.download(path, project_id)
        
        self.batch_files: List[RemoteBatchFile] = []
        self.file_lists: List[RemoteBatchFile] = []
        
        with gzip.open(path, "rt") as f:
            reader = csv.DictReader(f)
            
            for r in reader:
                if ".batch." in r["filename"]:
                    self.batch_files.append(RemoteBatchFile(os.path.basename(r["filename"]), r["url"], r["md5"], r["size"]))
                elif r["filename"].startswith("File_Lists/"):
                    self.file_lists.append(RemoteBatchFile(r["filename"], r["url"], r["md5"], r["size"]))

    @staticmethod
    def download(path: Path, project_id: str):        
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.suffix != ".gz":
            path = path.with_suffix(".csv.gz")
        
        with gzip.open(path, "wt") as f:
            writer = csv.writer(f)
            writer.writerow(["project", "project_id", "filename", "url", "md5", "size"])
            for r in to_table(get_all_data(project_id)):
                writer.writerow(astuple(r))


class AssemblyFileLists(_FileLists):
    def __init__(self, path: Path, override: bool = False) -> None:
        super().__init__(path, "zxfmy", override)


@dataclass(frozen=True)
class BaktaFile(_StringDataClass):
    sample: str
    status: str
    file_name: str
    file_md5: str
    tar_xz: str
    
    def __repr__(self) -> str:
        return f"BaktaFile({self.sample}, {self.status}, {self.file_name}, {self.file_md5}, {self.tar_xz})"


class BaktaManifest:
    def __init__(self, local: Path, remote: RemoteBatchFile):
        local.mkdir(parents=True, exist_ok=True)
        
        self.local = local / Path(remote.tar_xz).name
        
        if not self.local.exists():
            resp = requests.get(remote.tar_xz_url, stream=True)
            with self.local.open(mode="wb") as f:
                for chunk in resp.iter_content(chunk_size=8192):
                    f.write(chunk)

        self.data: List[BaktaFile] = []

        with gzip.open(self.local, "rt") as f:
            reader = csv.DictReader(f, delimiter=infer_delimiter(self.local))
            
            for r in reader:
                self.data.append(BaktaFile(r["sample"], r["status"], r["file_name"], r["file_md5"], r["tar_xz"]))


class BaktaFileLists(_FileLists):
    def __init__(self, path: Path, override: bool = False) -> None:
        super().__init__(path, "zt57s", override)
        
    def lists(self):
        for fl in self.file_lists:
            yield BaktaManifest(self.path.parent / "bakta_lists", fl)

