#!/usr/bin/env python3
from argparse import ArgumentParser

ORD = {
    'N': 0,
    'A': 1,
    'C': 2,
    'G': 3,
    'T': 4
}

COMPLEMENT = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A'
}

def parse_args():
    parser = ArgumentParser()
    parser.add_argument("-s", dest="seq", action="append", help="Sequence", required=True)
    parser.add_argument("-k", type=int, required=True, help="k-mer size")
    return parser.parse_args()

def reverse_complement(seq: str) -> str:
    return "".join(COMPLEMENT[c] for c in reversed(seq))

def get_kmers(seq: str, k: int):
    return [
        (tSeq[i:i+k], tSeq[i+k])
        for tSeq in [(k*'N') + seq + 'N', (k*'N') + reverse_complement(seq) + 'N']
        for i in range(0, len(tSeq)-k)
    ]

def kmer_to_key(kmer: str, edge: str, col: int):
    seq_key = "".join(
        str(ORD[cr])
        for seq in (kmer[::-1], edge)
        for cr in seq
    )
    return seq_key + str(col)

if __name__ == "__main__":
    args = parse_args()
    
    kmers = sorted(
        [
            (kmer, edge, i)
            for i, inSeq in enumerate(args.seq)
            for kmer, edge in get_kmers(inSeq, args.k)
        ],
        key = lambda d: kmer_to_key(*d))

    print("%d kmer-edge pairs" % len(kmers))

    for val, edge, col in kmers:
        print(val, edge, ORD[edge], col)
