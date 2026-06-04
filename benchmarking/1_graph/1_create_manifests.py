#!/usr/bin/env python3

from argparse import ArgumentParser
from dataclasses import dataclass
from itertools import chain, product
from pathlib import Path
import random
from typing import Dict, List, Set

from maki.atb.downloader import MakiAtbManifest
from maki.models.database.cluster import SequenceSourceDir


GENOME_COUNTS = [100, 500, 1000]
KMER_SIZES = [11, 21, 31, 41]
SUFFIX_SIZES = [4, 5, 6, 7, 8]
THREAD_COUNTS = [32, 64]


# Data members
# ------------


@dataclass
class AccessionList:
    name: str
    accessions: List[str]


# Auxilliary
# ----------


def group_by_taxonomy(mf: MakiAtbManifest) -> Dict[str, List[str]]:
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
    

# Main
# ----


def take_taxonomy_n(accns: List[str], taxonomy: str, n: int) -> AccessionList:
    assert len(accns) >= n, f"Not enough genomes belonging to group: {taxonomy}"
    return AccessionList(
        name=f"{taxonomy}_{n}",
        accessions=random.sample(accns, n)
    )


def generate_taxonomy_datasets(mf: MakiAtbManifest, taxonomy: str, seed: int):
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
        yield AccessionList(name=f"seed{seed}_{count}", accessions=accns[:count])


def main(manifest_path: Path, threads: int, outdir: Path):
    mf = MakiAtbManifest.from_csv(manifest_path, threads, outdir)
    groups = group_by_taxonomy(mf)

    # create datasets
    datasets = list(chain(
        generate_taxonomy_datasets(mf, "Klebsiella pneumoniae", 1),
        generate_taxonomy_datasets(mf, "Salmonella enterica", 2),
        generate_taxonomy_datasets(mf, "Mycobacterium tuberculosis", 3),
        generate_random_datasets(groups, 1),
        generate_random_datasets(groups, 2),
        generate_random_datasets(groups, 3)
    ))

    # download data
    accns = {
        accn
        for dset in datasets
        for accn in dset.accessions
    }
    
    # create source directories
    data_directories: List[str] = []
    
    with mf.retrieve_data(accns) as data:
        for dset in datasets:
            src = SequenceSourceDir(outdir / dset.name)
            
            # Write sequences
            for accn in dset.accessions:
                src.insert_genome(data[accn])
                
            # Write manifest
            with (src.source_dir / "manifest.txt").open("w") as f:
                src.write_manifest(f)
                
            data_directories.append(str(src.source_dir))

    # write parameter tsv
    with (outdir / "parameter-combinations.txt").open("w") as f:
        for dir, k, s, t in product(data_directories, KMER_SIZES, SUFFIX_SIZES, THREAD_COUNTS):
            f.write(f"{k}\t{s}\t{t}\t{dir}\n")

    print(f"Data and parameteres written to {outdir}.")


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
    
    main(Path(args.manifest), args.threads, out)
    