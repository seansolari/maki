from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
from maki.models.cluster import GenomeCluster
from maki.models.taxonomy import TaxonomyDB

class MetaGenomicDatabase:
    def __init__(self, root: str, kmer_size: int, mode="fixed"):
        self.root = Path(root)
        self.kmer_size = kmer_size
        self.mode = mode

        self.clusters = {}
        self.taxonomy = TaxonomyDB()
        self.metadata_file = self.root / "metadata.json"
        self.cluster_root = self.root / "clusters"

    @classmethod
    def load(cls, root: str):
        root = Path(root)
        meta = json.loads((root / "metadata.json").read_text())

        db = cls(root, meta["kmer_size"], meta["mode"])

        for cid, info in meta["clusters"].items():
            db.clusters[cid] = GenomeCluster(
                cid,
                info["tax_id"],
                info["rank"],
                db.kmer_size,
                db.cluster_root,
                mode=db.mode
            )

        db._rank = meta["rank"]
        db._manifest_hash = meta["manifest_hash"]

        return db
      
    # ====================
    # BUILD (Parallel)
    # ====================
    def build(self, manifest, rank, threads, force):
        groups = manifest.group_by_rank(self.taxonomy, rank)

        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = []

            for tax_id, records in groups.items():
                cid = f"{rank}_{tax_id}"

                cluster = GenomeCluster(
                    cid, tax_id, rank,
                    self.kmer_size,
                    self.cluster_root,
                    mode=self.mode
                )

                self.clusters[cid] = cluster
                futures.append(executor.submit(cluster.build, records, threads))

            for f in futures:
                f.result()

        self._save_metadata(manifest, rank)

    # ====================
    # TRUE INCREMENTAL UPDATE
    # ====================
    def update(self, new_manifest, threads, strategy, dry_run):
        if self.mode != "updatable":
            raise RuntimeError(
                "This database was built in 'fixed' mode and cannot be updated.\n"
                "Rebuild using:\n\n"
                "  metadb build --mode updatable ...\n"
            )

        new_hash = new_manifest.compute_hash()

        if new_hash == self._manifest_hash:
            print("No changes detected.")
            return

        new_groups = new_manifest.group_by_rank(self.taxonomy, self._rank)

        affected_clusters = []

        for tax_id, records in new_groups.items():
            cid = f"{self._rank}_{tax_id}"

            cluster = self.clusters.get(cid)

            if cluster is None:
                affected_clusters.append((cid, tax_id, records, "new"))
            else:
                if cluster.load_accessions() != {r.accession for r in records}:
                    affected_clusters.append((cid, tax_id, records, "modified"))

        if dry_run:
            print("Clusters to rebuild:")
            for cid, _, _, status in affected_clusters:
                print(f" - {cid} ({status})")
            return

        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = []

            for cid, tax_id, records, _ in affected_clusters:
                cluster = self.clusters.get(cid)

                if cluster is None:
                    cluster = GenomeCluster(
                        cid, tax_id, self._rank,
                        self.kmer_size,
                        self.cluster_root,
                        mode=self.mode
                    )
                    self.clusters[cid] = cluster

                    futures.append(executor.submit(cluster.build, records, threads))
                else:
                    futures.append(executor.submit(cluster.update, records, threads))

            for f in futures:
                f.result()

        self._save_metadata(new_manifest, self._rank)

    def _save_metadata(self, manifest, rank):
        meta = {
            "kmer_size": self.kmer_size,
            "rank": rank,
            "mode": self.mode,
            "manifest_hash": manifest.compute_hash(),
            "clusters": {
                cid: {
                    "tax_id": c.tax_id,
                    "rank": c.rank
                } for cid, c in self.clusters.items()
            }
        }

        with open(self.metadata_file, "w") as f:
            json.dump(meta, f, indent=2)
