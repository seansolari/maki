from dataclasses import dataclass
from typing import Tuple


@dataclass(frozen=True)
class Table:
    name: str

    @staticmethod
    def table_name() -> str: ...

    @staticmethod
    def columns() -> Tuple[str, ...]: ...
    
    @classmethod
    def default(cls):
        return cls(cls.table_name(), *cls.columns())
      
    def take_cols(self):
        for c in self.columns():
          yield (self.name, c)


@dataclass(frozen=True)
class SampleTable(Table):
    sample: str


@dataclass(frozen=True)
class AssemblyTable(SampleTable):
    species: str
    hq: str
    assembly_exists: str
    
    @staticmethod
    def table_name() -> str:
        return "assembly"
    
    @staticmethod
    def columns() -> Tuple[str, ...]:
        return ("assembly_accession", "asm_fasta_on_osf", "dataset", "scientific_name", "sylph_species", "hq_filter", "osf_tarball_filename", "osf_tarball_url", "comments")


@dataclass(frozen=True)
class AssemblyBatchTable(Table):
    tar_xz: str
    tar_xz_url: str
    tar_xz_md5: str
    tar_xz_size_MB: str
    
    @staticmethod
    def table_name() -> str:
        return "asm_batches"
    
    @staticmethod
    def columns() -> Tuple[str, ...]:
        return ("asm_tar_xz", "asm_tar_xz_url", "asm_tar_xz_md5", "asm_tar_xz_size_MB")


@dataclass(frozen=True)
class AssemblyStatsTable(SampleTable):
    total_length: str
    contigs: str
    mean_length: str
    longest: str
    shortest: str
    N_count: str
    Gaps: str
    N50: str
    N70: str
    N90: str
    
    @staticmethod
    def table_name() -> str:
        return "assembly_stats"
      
    @staticmethod
    def columns() -> Tuple[str, ...]:
        return ("sample_accession", "total_length", "number", "mean_length", "longest", "shortest", "N_count", "Gaps", "N50", "N70", "N90")


@dataclass(frozen=True)
class CheckM2Table(SampleTable):
    completeness: str
    contamination: str
    
    @staticmethod
    def table_name() -> str:
        return "checkm2"
      
    @staticmethod
    def columns() -> Tuple[str, ...]:
        return ("sample_accession", "Completeness_General", "Contamination")


@dataclass(frozen=True)
class BaktaTable(SampleTable):
    status: str
    file_name: str
    file_md5: str
    tar_xz: str
    tar_xz_md5: str
    tar_xz_size_MB: str
    
    @staticmethod
    def table_name() -> str:
        return "bakta"
      
    @staticmethod
    def columns() -> Tuple[str, ...]:
        return ("sample", "status", "file_name", "file_md5", "tar_xz", "tar_xz_md5", "tar_xz_size_MB")


@dataclass(frozen=True)
class BaktaBatchTable(Table):
    tar_xz: str
    tar_xz_url: str
    tar_xz_md5: str
    tar_xz_size_MB: str
    
    @staticmethod
    def table_name() -> str:
        return "ann_batches"
    
    @staticmethod
    def columns() -> Tuple[str, ...]:
        return ("ann_tar_xz", "ann_tar_xz_url", "ann_tar_xz_md5", "ann_tar_xz_size_MB")


ASSEMBLY_SCHEMA = AssemblyTable.default()
ASSEMBLY_BATCH_SCHEMA = AssemblyBatchTable.default()
ASSEMBLY_STATS_SCHEMA = AssemblyStatsTable.default()
CHECKM2_SCHEMA = CheckM2Table.default()
BAKTA_SCHEMA = BaktaTable.default()
BAKTA_BATCH_SCHEMA = BaktaBatchTable.default()
