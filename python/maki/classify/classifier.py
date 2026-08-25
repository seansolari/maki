
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from maki.classify.sample_manifest import Sample, SampleManifest
from maki.models.database.clustering import ClusterHandle
from maki.models.database.database import StaticDatabase
import maki.core as mx


class ClassificationResult:
    def __init__(self, assignments):
        self.assignments = assignments

    def to_json(self):
        import json
        return json.dumps(self.assignments, indent=2)


class ClassifyThread:
    def __init__(self, root: Path, sample: Sample) -> None:
        self.root = root
        self.sample = sample
        
        self.graph_dir = root / "index"
        self.classify_dir = root / "classify"
        
    @property
    def indexed(self) -> bool:
        return (self.graph_dir / "index" / ".ready").exists()
        
    def index(self, kmer_size: int, threads: int):
        """Create sample graph at `self.graph_dir`.
        """
        if not self.indexed:
            rp = mx.read_pair(str(self.sample.forward), str(self.sample.reverse))
            opts = mx.build_opts(kmer_size, 7, self.graph_dir, threads)
            
            mx.construct_wdbg(rp, opts)
            
            (self.graph_dir / "index" / ".ready").touch()
        
    def classify(self, h: ClusterHandle, threads: int):
        outdir = self.classify_dir / h.name

        if outdir.exists():
            print(f"[classify] classification {self.graph_dir.name} -> {h.name} exists at {outdir}, skipping.")
            return
        
        
        outdir.mkdir(parents=True, exist_ok=True)


class Classifier:
    def __init__(self, db: StaticDatabase, threads: int = 4):
        self.db = db
        self.threads = threads

    def run(self, sample_manifest: SampleManifest, output_dir: Path):
        output_dir.mkdir(parents=True, exist_ok=True)

        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            futures = []

            for sample in sample_manifest.samples:
                futures.append(
                    executor.submit(self._process_sample, sample, output_dir)
                )

            for f in futures:
                f.result()

    def _process_sample(self, sample, output_dir):
        result = self.classify_sample(sample)

        out_file = output_dir / f"{sample.name}.json"
        with open(out_file, "w") as f:
            f.write(result.to_json())

    def classify_sample(self, sample):
        reads = (sample.forward, sample.reverse)

        reads = self._run_filter(reads)

        assignments = {}

        for cid in self.db.clusters:
            cluster = self.db.clusters.get(cid)
            assignments[cid] = self._classify_cluster(cluster, reads)

        return ClassificationResult(assignments)

    def _run_filter(self, reads):
        return reads

    def _classify_cluster(self, cluster, reads):
        raise NotImplementedError