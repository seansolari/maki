#!/usr/bin/env python3

from argparse import ArgumentParser
from pathlib import Path
import random
from typing import Dict, List, Set

from maki.atb.async_downloader import MakiAtbManifest
from maki.models.database.cluster import SequenceSourceDir
import maki.core as mx


GENOME_COUNTS = [100, 500, 1000]
KMER_SIZES = [11, 21, 31, 41]
SUFFIX_SIZES = [4, 5, 6, 7, 8]
THREAD_COUNTS = [32, 64]


def take_taxonomy(mf: MakiAtbManifest, taxonomy: str, seed: int = 1) -> List[List[str]]:
    accns = sorted({r.accession for r in mf.records() if r.taxonomy == taxonomy})
    assert len(accns) >= GENOME_COUNTS[-1], f"Not enough genomes belonging to group: {taxonomy}"
    
    random.seed(seed)
    return [
        random.sample(accns, count)
        for count in GENOME_COUNTS
    ]
    
  
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


def take_random(mf: MakiAtbManifest, reps: int, seed: int = 1) -> List[List[List[str]]]:
    groups = group_by_taxonomy(mf)

    result: List[List[List[str]]] = []
    for rep in range(reps):
        random.seed(seed + rep)
        keys = list(groups.keys())

        rep_result: List[List[str]] = []
        for count in GENOME_COUNTS:
            random.shuffle(keys)
            rep_result.append(sample_accessions(keys, groups, count))
            
        result.append(rep_result)
        
    return result


def main(manifest_path: Path, threads: int, outdir: Path):
    manifest = MakiAtbManifest.from_csv(manifest_path, threads, outdir)

    # create low-diversity genome samples
    kp_accns = take_taxonomy(manifest, "Klebsiella pneumoniae", 1)
    se_accns = take_taxonomy(manifest, "Salmonella enterica", 2)
    mt_accns = take_taxonomy(manifest, "Mycobacterium tuberculosis", 3)
    
    # random sampling
    rnd_accns = take_random(manifest, 3, 4)
    
    # download data
    all_accns: Set[str] = set()
    for grp in kp_accns:
        all_accns |= set(grp)
    for grp in se_accns:
        all_accns |= set(grp)
    for grp in mt_accns:
        all_accns |= set(grp)
    for batch in rnd_accns:
        for grp in batch:
            all_accns |= set(grp)
            
    with manifest.retrieve_data(all_accns) as data:
        for grp in kp_accns:
            src = SequenceSourceDir(outdir / f"kp_{len(grp)}", mx.FileType.GFF3)
            src.insert(data.select(grp))
            # write manifest
            ...
            
    # write parameter YAML
    ...
            

def parse_args():
    parser = ArgumentParser()
    parser.add_argument("--manifest", type=str)
    parser.add_argument("--threads", type=int)
    parser.add_argument("--out", type=str)
    return parser.parse_args()
    

if __name__ == "__main__":
    args = parse_args()
    main(Path(args.manifest), args.threads, Path(args.out))
    