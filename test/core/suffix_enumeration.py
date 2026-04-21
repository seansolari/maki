#!/usr/bin/env python3

from argparse import ArgumentParser


DNA = {0: "A", 1: "C", 2: "G", 3: "T"}

"""Convert binary representation of DNA sequence into string"""
def to_string(x: int, size: int) -> str:
  s = ""
  for i in range(size):
    bp = (x >> (2 * i)) & 0b11
    s += DNA[bp]
  return s

"""Number of DNA suffixes up to and including a given size"""
def num_suffixes(size: int) -> int:
  return ((1 << (2 * (size + 1))) - 1) // 3

def colex_from_rank(r: int, s: int):
  x, l = 0, s
  while r > 0:
    assert(l > 0)
    ls = num_suffixes(l-1)
    f = (r-1) // ls
    x = (x << 2) | f
    r -= (f * ls) + 1
    l -= 1
  return to_string(x, s-l)

"""Convert rank representation of DNA suffix into string"""
def from_rank(r: int) -> str:
  size = 0
  while num_suffixes(size) < r:
    size += 1
  if size == 0:
    return ""
  else:
    return to_string(r - num_suffixes(size - 1) - 1, size)

def parse_args():
  parser = ArgumentParser()
  parser.add_argument("size", type=int, help="List suffixes up to a given size")
  return parser.parse_args()

def main(s: int):
  for r in range(num_suffixes(s)):
    print(r, colex_from_rank(r, s), sep = "\t")

if __name__ == "__main__":
  args = parse_args()
  main(args.size)
