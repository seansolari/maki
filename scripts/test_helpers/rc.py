#!/usr/bin/env python3
import sys

COMPLEMENT = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A'
}

def reverse_complement(seq: str) -> str:
  return "".join(COMPLEMENT[c] for c in seq[::-1])

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("[error] expects one argument!")
        exit(1)

    seq = sys.argv[1]

    if not all(c in COMPLEMENT for c in seq):
        print("[error] must be DNA sequence: %s" % seq)
        
    rc = reverse_complement(seq)
    print(rc)
