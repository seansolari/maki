#!/usr/bin/env python3
from __future__ import annotations
from argparse import ArgumentParser

O = {
    'A': 0,
    'C': 1,
    'G': 2,
    'T': 3
}

def parse_args():
    parser = ArgumentParser()
    parser.add_argument("--seq", dest="seq", action="append", help="Sequence", required=True)
    parser.add_argument("-k", type=int, required=True, help="k-mer size")
    parser.add_argument("-s", type=int, required=True, help="Suffix size")
    return parser.parse_args()

def seqToInt(seq: str) -> int:
    val = 0
    for i, s in enumerate(seq):
        val += O[s] * (4**i)
    return val

def countSuffixes(seq: str, k: int, s: int):
    table = [0 for _ in range(4**s)]
    for i in range(k, len(seq)+1):
        table[seqToInt(seq[i-s:i])] += 1
    return table

def sumTables(tables):
    tot = [0 for _ in tables[0]]
    for table in tables:
        for i in range(len(table)):
            tot[i] += table[i]
    return tot

if __name__ == "__main__":
    args = parse_args()
    k, s = int(args.k), int(args.s)
    tables = [
        countSuffixes(seq, k, s)
        for seq in args.seq
    ]
    total = sumTables(tables)
    for table in tables:
        print(table)
    print("TOTAL:")
    print(total)   
