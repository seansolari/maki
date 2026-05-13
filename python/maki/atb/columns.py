
from dataclasses import dataclass


@dataclass
class SampleColumns:
    table: str
    sample: str
    species: str
    hq: str
    assembly: str


@dataclass
class BaktaColumns:
    table: str
    sample: str
    status: str
    file_name: str
    file_md5: str
    tar_xz: str
    tar_xz_md5: str
    tar_xz_size_MB: str


def resolve_columns(schema) -> SampleColumns:
    return SampleColumns("assembly", "sample_accession", "sylph_species", "hq_filter", "asm_fasta_on_osf")


def resolve_bakta_columns(schema) -> BaktaColumns:
    return BaktaColumns("bakta", "sample", "status", "file_name", "file_md5", "tar_xz", "tar_xz_md5", "tar_xz_size_MB")
