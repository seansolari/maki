import tarfile
from pathlib import Path


class GenomeCluster:
    def __init__(self, cluster_id, tax_id, rank, kmer_size, root, mode="fixed"):
        self.cluster_id = cluster_id
        self.tax_id = tax_id
        self.rank = rank
        self.kmer_size = kmer_size
        self.mode = mode

        self.path = root / cluster_id
        self.path.mkdir(parents=True, exist_ok=True)

        self.index_path = self.path / "index.dat"
        self.accession_file = self.path / "accessions.txt"
        self.source_archive = self.path / "source.tar.gz"

    # ======================
    # BUILD
    # ======================
    def build(self, records, threads):
        fasta_files = [r.fasta for r in records]

        self._build_index_backend(fasta_files, threads)

        # Store accessions
        self.save_accessions([r.accession for r in records])

        # Store compressed sources if updatable
        if self.mode == "updatable":
            self._compress_sources(records)

    # ======================
    # UPDATE
    # ======================
    def update(self, new_records, threads):
        if self.mode != "updatable":
            raise RuntimeError(
                f"Cluster '{self.cluster_id}' is not updatable.\n"
                "This database was built in 'fixed' mode.\n"
                "Rebuild the database with '--mode updatable' to enable updates."
            )

        # Step 1: extract existing sources
        extracted_dir = self._extract_sources()

        # Step 2: merge existing + new
        all_records = self._merge_sources(extracted_dir, new_records)

        # Step 3: rebuild index
        fasta_files = [r.fasta for r in all_records]
        self._build_index_backend(fasta_files, threads)

        # Step 4: recompress
        self._compress_sources(all_records)

        # cleanup
        self._cleanup_temp(extracted_dir)

    # ======================
    # SOURCE HANDLING
    # ======================
    def _compress_sources(self, records):
        with tarfile.open(self.source_archive, "w:gz") as tar:
            for r in records:
                tar.add(r.fasta, arcname=Path(r.fasta).name)
                tar.add(r.gff, arcname=Path(r.gff).name)

    def _extract_sources(self):
        extract_dir = self.path / "source_tmp"
        extract_dir.mkdir(exist_ok=True)

        with tarfile.open(self.source_archive, "r:gz") as tar:
            tar.extractall(extract_dir)

        return extract_dir

    def _merge_sources(self, extracted_dir, new_records):
        """
        Combine existing extracted files with new records.
        Returns combined GenomeRecord list.
        """
        combined = []

        # Existing
        for f in extracted_dir.glob("*"):
            if f.suffix in (".fa", ".fasta"):
                combined.append(f)

        # New
        combined.extend(new_records)

        return new_records  # simplified placeholder

    def _cleanup_temp(self, path):
        import shutil
        shutil.rmtree(path, ignore_errors=True)

    # ======================
    # ACCESSIONS
    # ======================
    def save_accessions(self, accessions):
        with open(self.accession_file, "w") as f:
            f.write("\n".join(accessions))

    def load_accessions(self):
        if not self.accession_file.exists():
            return set()
        return set(self.accession_file.read_text().splitlines())

    # ======================
    # PLUGIN HOOK
    # ======================
    def _build_index_backend(self, fasta_files, threads):
        raise NotImplementedError
      
      
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
