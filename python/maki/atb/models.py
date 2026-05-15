from dataclasses import dataclass, fields
from typing import Tuple


@dataclass(frozen=True)
class Table:
    name: str

    @classmethod
    def columns(cls) -> Tuple[str, ...]: ...
    
    @staticmethod
    def table_name() -> str: ...
    
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
    
    @classmethod
    def columns(cls) -> Tuple[str, ...]:
        return ("assembly_accession", "asm_fasta_on_osf", "dataset", "scientific_name", "sylph_species", "hq_filter", "osf_tarball_filename", "osf_tarball_url", "comments")
    
    @staticmethod
    def table_name() -> str:
        return "assembly"


@dataclass(frozen=True)
class AssemblyBatchTable(Table):
    name: str
    tar_xz: str
    tar_xz_url: str
    tar_xz_md5: str
    tar_xz_size_MB: str
    
    @classmethod
    def columns(cls) -> Tuple[str, ...]:
        return ("asm_tar_xz", "asm_tar_xz_url", "asm_tar_xz_md5", "asm_tar_xz_size_MB")
    
    @staticmethod
    def table_name() -> str:
        return "asm_batches"
    

@dataclass(frozen=True)
class AssemblyStatsTable(Table):
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
    
    @classmethod
    def default(cls):
        return cls(name="assembly_stats", sample="sample_accession", total_length="total_length", contigs="number", mean_length="mean_length", longest="longest", shortest="shortest", N_count="N_count", Gaps="Gaps", N50="N50", N70="N70", N90="N90")

    def take_cols(self):
        for c in ("sample_accession", "total_length", "number", "mean_length", "longest", "shortest", "N_count", "Gaps", "N50", "N70", "N90"):
            yield (self.name, c)


@dataclass(frozen=True)
class CheckM2Table(Table):
    completeness: str
    contamination: str
    
    @classmethod
    def default(cls):
        return cls(name="checkm2", sample="sample_accession", completeness="Completeness_General", contamination="Contamination")

    def take_cols(self):
        for c in ("sample_accession", "Completeness_General", "Contamination"):
            yield (self.name, c)


@dataclass(frozen=True)
class BaktaTable(Table):
    status: str
    file_name: str
    file_md5: str
    tar_xz: str
    tar_xz_md5: str
    tar_xz_size_MB: str

    @classmethod
    def default(cls):
        return cls(name="bakta", sample="sample", status="status", file_name="file_name", file_md5="file_md5", tar_xz="tar_xz", tar_xz_md5="tar_xz_md5", tar_xz_size_MB="tar_xz_size_MB")

    def take_cols(self):
        for c in ("sample", "status", "file_name", "file_md5", "tar_xz", "tar_xz_md5", "tar_xz_size_MB"):
            yield (self.name, c)
  

ASSEMBLY_SCHEMA = AssemblyTable.default()
ASSEMBLY_STATS_SCHEMA = AssemblyStatsTable.default()
CHECKM2_SCHEMA = CheckM2Table.default()
BAKTA_SCHEMA = BaktaTable.default()
