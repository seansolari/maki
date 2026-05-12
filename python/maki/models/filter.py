class FilterCluster:
    def __init__(self, root, kmer_size: int):
        self.path = root / "filter"
        self.index_path = self.path / "index.dat"
        self.kmer_size = kmer_size

        self.path.mkdir(parents=True, exist_ok=True)

    def build(self, fasta_files, threads: int):
        # PLUGIN HOOK
        self._build_index_backend(fasta_files, threads)

    def classify(self, reads):
        # PLUGIN HOOK
        return self._classify_backend(reads)

    def _build_index_backend(self, fasta_files, threads):
        raise NotImplementedError

    def _classify_backend(self, reads):
        raise NotImplementedError