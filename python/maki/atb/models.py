from dataclasses import dataclass


@dataclass(frozen=True)
class Table:
    name: str
    sample: str


@dataclass(frozen=True)
class AssemblyTable(Table):
    species: str
    hq: str
    assembly_exists: str
    
    @classmethod
    def default(cls):
        return cls(name="assembly", sample="sample_accession", species="sylph_species", hq="hq_filter", assembly_exists="asm_fasta_on_osf")


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
        return cls(name="assembly_stats", sample="sample_accession", total_length="total_length", contigs="contigs", mean_length="mean_length", longest="longest", shortest="shortest", N_count="N_count", Gaps="Gaps", N50="N50", N70="N70", N90="N90")


@dataclass(frozen=True)
class CheckM2Table(Table):
    completeness: str
    contamination: str
    
    @classmethod
    def default(cls):
        return cls(name="checkm2", sample="sample_accession", completeness="Completeness_General", contamination="Contamination")


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


ASSEMBLY_SCHEMA = AssemblyTable.default()
ASSEMBLY_STATS_SCHEMA = AssemblyStatsTable.default()
CHECKM2_SCHEMA = CheckM2Table.default()
BAKTA_SCHEMA = BaktaTable.default()
