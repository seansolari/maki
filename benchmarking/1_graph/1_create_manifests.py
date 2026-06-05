#!/usr/bin/env python3

from argparse import ArgumentParser
from dataclasses import dataclass
from itertools import chain, product
import logging
from pathlib import Path
import random
from typing import Dict, List, Set

from maki.atb.downloader import DiskAtbManifest
from maki.models.database.cluster import SequenceSourceDir


GENOME_COUNTS = [100, 500, 1000]
KMER_SIZES = [11, 21, 31, 41]
SUFFIX_SIZES = [4, 5, 6, 7, 8]
THREAD_COUNTS = [32, 64]

logger = logging.getLogger(__name__)


# Data members
# ------------


@dataclass
class AccessionList:
    name: str
    accessions: List[str]


# Auxilliary
# ----------


def group_by_taxonomy(mf: DiskAtbManifest) -> Dict[str, List[str]]:
    res = {}
    for rec in mf.records():
        res.setdefault(rec.taxonomy, []).append(rec.accession)
    return res
    
    
def sample_accessions(keys: List[str], groups: Dict[str, List[str]], n: int) -> List[str]:
    failed = 0
    selected: Set[str] = set()
    while len(selected) < n:
        for key in keys:
            nacn = random.choice(groups[key])
            if nacn not in selected:
                selected.add(nacn)
                if len(selected) == n:
                    break
            else:
                failed += 1
                if failed == 1000:
                    raise RuntimeError("Reached maximum attempted inserts")
    return sorted(selected)
  

def gather_accessions(datasets: List[SequenceSourceDir]) -> Set[str]:
    accessions: Set[str] = set()
    
    for dataset in datasets:
        accessions |= set(dataset.accessions())
        
    return accessions
    

# Main
# ----


def take_taxonomy_n(accns: List[str], taxonomy: str, n: int) -> AccessionList:
    assert len(accns) >= n, f"Not enough genomes belonging to group: {taxonomy}"
    
    logger.info("Sampling %d genomes from taxonomic group %s (which contains %d genomes).", n, taxonomy, len(accns))
    return AccessionList(
        name=f"{taxonomy.replace(" ", "")}_{n}",
        accessions=random.sample(accns, n)
    )


def generate_taxonomy_datasets(mf: DiskAtbManifest, taxonomy: str, seed: int):
    accns = sorted({r.accession for r in mf.records() if r.taxonomy == taxonomy})
    
    random.seed(seed)
    for count in GENOME_COUNTS:
        yield take_taxonomy_n(accns, taxonomy, count)


def generate_random_datasets(groups: Dict[str, List[str]], seed: int):
    random.seed(seed)
    
    keys = list(groups.keys())
    random.shuffle(keys)
    accns = sample_accessions(keys, groups, GENOME_COUNTS[-1])
    
    for count in GENOME_COUNTS:
        logger.info("Sampling %d genomes according to broad taxonomic distribution (seed=%d)", count, seed)
        yield AccessionList(name=f"seed{seed}_{count}", accessions=accns[:count])


def main(manifest_path: Path, threads: int, outdir: Path):
    logger.info("Reading manifest from %s.", manifest_path)
    mf = DiskAtbManifest.from_csv(manifest_path, threads, outdir)
    
    groups = group_by_taxonomy(mf)

    logger.info("Creating benchmarking datasets")
    datasets = list(chain(
        generate_taxonomy_datasets(mf, "Klebsiella pneumoniae", 1),
        generate_taxonomy_datasets(mf, "Salmonella enterica", 2),
        generate_taxonomy_datasets(mf, "Mycobacterium tuberculosis", 3),
        generate_random_datasets(groups, 1),
        generate_random_datasets(groups, 2),
        generate_random_datasets(groups, 3)
    ))

    # download data
    sources = [SequenceSourceDir(outdir / dset.name) for dset in datasets]
    
    accns = {accn for dset in datasets for accn in dset.accessions}
    missing = accns - gather_accessions(sources)
    
    if missing:
        logger.info("Downloading data for %d accessions (out of %d total)", len(missing), len(accns))
        
        with mf.retrieve_data(missing) as data:
            for dset, src in zip(datasets, sources):
                changed = False
                
                # Write sequences
                for accn in dset.accessions:
                    if accn in missing:
                        src.insert_genome(data[accn])
                        changed = True
                
                if changed:
                    logger.info("Writing manifest to %s", src.source_dir / "manifest.txt")
                    with (src.source_dir / "manifest.txt").open("w") as f:
                        src.write_manifest(f)

    # write parameter tsv
    with (outdir / "parameter-combinations.txt").open("w") as f:
        for k, s, t, src in product(KMER_SIZES, SUFFIX_SIZES, THREAD_COUNTS, sources):
            f.write(f"{k}\t{s}\t{t}\t{src.source_dir}\n")

    logger.info("Data and parameteres written to %s", outdir / "parameter-combinations.txt")


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("--manifest", type=str)
    parser.add_argument("--threads", type=int)
    parser.add_argument("--out", type=str)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    
    logging.basicConfig(
        filename=out / "create-manifests.log",
        format="%(asctime)s [%(name)s]:%(levelname)s %(message)s",
        encoding='utf-8',
        level=logging.DEBUG,
        force=True
    )
    
    main(Path(args.manifest), args.threads, out)
    