#!/usr/bin/env python3
from argparse import ArgumentParser
from typing import List, TypeVar

ENCODER = {
    'A': "00",
    'C': "01",
    'G': "10",
    'T': "11"
}

EDGE_ENCODE = {
    'A': "000",
    'C': "001",
    'G': "010",
    'T': "011",
    'N': "111"
}

COMPLEMENT = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A'
}

def parse_args():
    parser = ArgumentParser()
    parser.add_argument("-k", type=int, required=True, help="k-mer size")
    parser.add_argument("-s", dest="seq", type=str, required=True, help="Sequence string")
    parser.add_argument("--suffix", dest="suffix", type=str, default=None, help="Filter by suffix")
    return parser.parse_args()

def yield_kmers(seq: str, k: int):
    for i in range(0, len(seq) - k):
        yield seq[i:i+k], seq[i+k]

T = TypeVar('T')

def yield_bytes(seq: List[T]):
    numBytes = (len(seq) + 3) // 4
    for i in range(numBytes):
        yield reversed(seq[4*i:4*(i+1)])

def print_kmers(input_sequence: str, input_k: int):
    for kmer, edge in yield_kmers(input_sequence + 'N', input_k):
        enc = [ENCODER[c] for c in kmer]

        print(
            " ".join("".join(b) for b in yield_bytes(enc)),
            "->", EDGE_ENCODE[edge],
            "(", kmer, "->", edge, ")")
        
def print_kmers_with_suffix(input_sequence: str, input_k: int, suffix: str):
    for kmer, edge in yield_kmers(input_sequence + 'N', input_k):
        if not kmer.endswith(suffix):
            continue

        enc = [ENCODER[c] for c in kmer]

        print(
            " ".join("".join(b) for b in yield_bytes(enc)),
            "->", EDGE_ENCODE[edge],
            "(", kmer, "->", edge, ")")

def reverse_complement(seq: str) -> str:
    return "".join(COMPLEMENT[c] for c in reversed(seq))

if __name__ == "__main__":
    args = parse_args()
    
    print_function = print_kmers if args.suffix is None else lambda seq, k: print_kmers_with_suffix(seq, k, args.suffix)

    print("[main] using k=%d" % args.k)

    print("[FORWARD]")
    print_function(args.seq, args.k)

    print("[REVERSE COMPLEMENT]")
    print_function(reverse_complement(args.seq), args.k)
