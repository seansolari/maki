
from enum import Enum
from pathlib import Path
import typing


class cdbg_files:
    def __init__(self, base: Path) -> None: ...
    
    
class cdbg:
    def __init__(self) -> None: ...
    
    @staticmethod
    def disk_load(graph: cdbg, path: Path) -> None: ...


class wdbg_files:
    def __init__(self, base: Path) -> None: ...
    

class wdbg:
    def __init__(self) -> None: ...
    
    @staticmethod
    def disk_load(graph: wdbg, path: Path) -> None: ...


class FileType(Enum):
    FNA: FileType
    GFF3: FileType
    FQ: FileType
    

class genome_manifest:
    def __init__(self) -> None: ...
    
    @property
    def files(self) -> list[str]: ...
    
    @property
    def type(self) -> FileType: ...
    
    
def read_manifest(file: str, filter: FileType) -> genome_manifest: ...


class read_pair:    
    @typing.overload
    def __init__(self) -> None: ...
    
    @typing.overload
    def __init__(self, forward: str, reverse: str) -> None: ...
    
    @property
    def forward(self) -> str: ...
    
    @property
    def reverse(self) -> str: ...


class build_opts:
    @typing.overload
    def __init__(self) -> None: ...
    
    @typing.overload
    def __init__(self, k: int, s: int, out: Path, threads: int) -> None: ...
    
    @property
    def kmer_size(self) -> int: ...
    
    @property
    def suffix_size(self) -> int: ...
    
    @property
    def out(self) -> Path: ...
    
    @property
    def pool_size(self) -> int: ...
    
    @property
    def reserve_per_chunk(self) -> int: ...
    
    @property
    def threads(self) -> int: ...
    

def construct_cdbg(manifest: genome_manifest, opts: build_opts) -> cdbg_files: ...

def construct_wdbg(reads: read_pair, opts: build_opts) -> wdbg_files: ...
