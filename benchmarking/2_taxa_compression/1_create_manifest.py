#!/usr/bin/env python3

from argparse import ArgumentParser
import gzip
import logging
from pathlib import Path
import random
from typing import List

from maki.atb.downloader import DiskAtbManifest
from maki.models.database.cluster import SequenceSourceDir


logger = logging.getLogger(__name__)


# Auxilliary
# ----------


def sample_accessions(mf : DiskAtbManifest, n: int, seed: int) -> List[str]:
    accns = sorted({rec.accession for rec in mf.records()})
    
    if len(accns) < n:
        raise RuntimeError(f"Requested {n} genomes when only {len(accns)} records exist in manifest.")
    
    random.seed(seed)
    return random.sample(accns, n)
    


# Data members
# ------------


def main(manifest_path: Path, n: int, seed: int, threads: int, outdir: Path):
    logger.info("Reading manifest from %s.", manifest_path)
    mf = DiskAtbManifest.from_csv(manifest_path, threads, outdir)
    
    # download data
    src = SequenceSourceDir(outdir / f"n{n}_seed{seed}")
    
    accns = sample_accessions(mf, n, seed)
    missing = set(accns) - set(src.accessions())
    
    if missing:
        with mf.retrieve_data(missing) as data:
            src.pinsert([data[accn] for accn in missing], threads)
    
    # update manifest
    logger.info("Writing manifest to %s", src.source_dir / "manifest.txt")
    with (src.source_dir / "manifest.txt").open("w") as f:
        src.write_manifest(f)

    # write complete manifest
    with gzip.open(outdir / f"n{n}_seed{seed}.mani.csv.gz", "wt") as f:
        f.write(f"accession,taxonomy,fasta,gff\n")
        for rec in mf.select(accns).records():
            f.write(f"{rec.accession},{rec.taxonomy},,{src.source_dir.absolute() / f"{rec.accession}.gff.gz"}\n")
       
        
def parse_args():
    parser = ArgumentParser()
    parser.add_argument("--manifest", type=str)
    parser.add_argument("--n", type=int)
    parser.add_argument("--seed", type=int)
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
    
    main(Path(args.manifest), args.n, args.seed, args.threads, out)
        