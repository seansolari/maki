
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


class ClassificationResult:
    def __init__(self, assignments):
        self.assignments = assignments

    def to_json(self):
        import json
        return json.dumps(self.assignments, indent=2)
      

class Classifier:
    def __init__(self, db, threads=4, confidence=0.1, min_hits=5):
        self.db = db
        self.threads = threads
        self.confidence = confidence
        self.min_hits = min_hits

    def run(self, sample_manifest, output_dir):
        output_dir = Path(output_dir)
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

        for cid, cluster in self.db.clusters.items():
            assignments[cid] = self._classify_cluster(cluster, reads)

        return ClassificationResult(assignments)

    def _run_filter(self, reads):
        return reads  # 🔌

    def _classify_cluster(self, cluster, reads):
        raise NotImplementedError  # 🔌