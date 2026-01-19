#!/usr/bin/env python3
from __future__ import annotations
from argparse import ArgumentParser
import gzip
import multiprocessing
import os
import re
from typing import Mapping, Set, Tuple


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("--manifest", dest="manifest", type=str, required=True)
    parser.add_argument("--prefix", dest="prefix", type=str, required=True)
    return parser.parse_args()

def read_lines(fpath: str):
    with open(fpath, 'r') as f:
        return [line.strip() for line in f]
    
def summarise_bakta_inferences(fpath: str) -> Tuple[Set[str], Set[str]]:
    uref90, uref50 = set(), set()

    keys: Mapping[str, int] = {}
    with gzip.open(fpath, 'rb') as f:
        for line in f:
            if line.startswith(b'#Sequence Id'):
                for i, k in enumerate(line.decode("utf-8").strip().split("\t")):
                    keys[k] = i
            elif not line.startswith(b"#"):
                parts = line.decode("utf-8").split("\t")
                if parts[keys["Type"]] == "cds":
                    for uniref_xref in re.findall(r"UniRef\d+_[^,\s]+", parts[keys["DbXrefs"]]):
                        uniref_aa, uniref_cluster = uniref_xref.split("_")
                        if uniref_aa == "UniRef90":
                            uref90.add(uniref_cluster)
                        elif uniref_aa == "UniRef50":
                            uref50.add(uniref_cluster)

    return uref90, uref50

if __name__ == "__main__":
    args = parse_args()
    files = read_lines(args.manifest)
    with multiprocessing.Pool() as pool:
        data = pool.map(summarise_bakta_inferences, files)

    uref90: Mapping[str, int] = {}
    uref50: Mapping[str, int] = {}
    for g90, g50 in data:
        for g in g90:
            try:
                uref90[g] += 1
            except KeyError:
                uref90[g] = 1
        for g in g50:
            try:
                uref50[g] += 1
            except KeyError:
                uref50[g] = 1
        
        g90.clear()
        g50.clear()

    with gzip.open(args.prefix + "_UniRef90.txt.gz", 'wb') as f:
        for k, v in uref90.items():
            line = "%s,%d\n" % (k, v)
            f.write(line.encode("utf-8"))
    
    with gzip.open(args.prefix + "_UniRef50.txt.gz", 'wb') as f:
        for k, v in uref50.items():
            line = "%s,%d\n" % (k, v)
            f.write(line.encode("utf-8"))
