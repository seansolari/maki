#!/usr/bin/env python3
from argparse import ArgumentParser
import gzip
import multiprocessing
import os
from typing import Dict, List, Mapping, Set


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("--manifest", dest="manifest", type=str, required=True)
    parser.add_argument("--out", dest="out", type=str, required=True)
    return parser.parse_args()

def read_lines(fpath: str):
    with open(fpath, 'r') as f:
        return [line.strip() for line in f]
    
def parse_bakta_txt_file(fpath: str):
    data: Mapping[str, str] = {}
    data["genome"] = os.path.basename(os.path.dirname(fpath)).rstrip("_genomic")
    with gzip.open(fpath, 'rb') as f:
        for line in f:
            try:
                k, v = line.decode("utf-8").strip().split(":")
                if v:
                    data[k] = v.strip()
            except ValueError:
                continue
    return data

def common_keys(data: List[Dict[str, str]]) -> Set[str]:
    _keys = set()
    for dct in data:
        if len(_keys) == 0:
            _keys = set(dct.keys())
        else:
            _keys = _keys & dct.keys()
    return _keys

if __name__ == "__main__":
    args = parse_args()
    files = read_lines(args.manifest)
    with multiprocessing.Pool() as pool:
        data = pool.map(parse_bakta_txt_file, files)
    keys = sorted(common_keys(data) - {"DOI", "Database", "Software", "URL"})
    with open(args.out, 'w') as f:
        f.write("%s\n" % ",".join(keys))
        for dct in data:
            f.write("%s\n" % (",".join(dct[k] for k in keys)))
