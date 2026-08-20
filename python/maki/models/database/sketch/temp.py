



@dataclass(frozen=True)
class GenomeInput:
    genome_id: str
    fasta_path: str


@dataclass(frozen=True)
class GenomeRecord:
    genome_id: str
    fasta_path: str
    sketch_path: str
    cluster_id: int

