#!/usr/bin/env python3
from argparse import ArgumentParser
from typing import List, TypeVar

ENCODER = {
    '$': "00",
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
    return parser.parse_args()

def yield_terminals(seq: str, k: int):
    for i in range(k):
        yield (k-i) * "$" + seq[:i], seq[i], i

T = TypeVar('T')

def yield_bytes(seq: List[T]):
    numBytes = (len(seq) + 3) // 4
    for i in range(numBytes):
        yield reversed(seq[4*i:4*(i+1)])

def reverse_complement(seq: str) -> str:
    return "".join(COMPLEMENT[c] for c in reversed(seq))

def print_terminals(sequence: str, k: int):
    for terminal, edge, tsize in yield_terminals(sequence, k):
        enc_edge = EDGE_ENCODE[edge]
        tsize_bits = '{0:08b}'.format(tsize)

        print(
            ' '.join(''.join(ENCODER[c] for c in b) for b in yield_bytes(terminal)), '-', tsize_bits, '-', enc_edge,
            "\t|\t",
            terminal, '-', tsize, '-', edge)

if __name__ == "__main__":
    args = parse_args()

    print("[main] using k=%d" % args.k)

    print("[FORWARD]")
    print_terminals(args.seq, args.k)

    print("[REVERSE]")
    print_terminals(reverse_complement(args.seq), args.k)
