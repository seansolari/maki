#!/usr/bin/env python3

from argparse import ArgumentParser
from typing import Callable, List, Sequence, Tuple, TypeVar

from .sdbg import construct_sdbg, succinct_graph_container
from .rc import reverse_complement

T = TypeVar("T")

def remove_consecutive_duplicates(lst):
    if not lst:
        return

    write = 1

    for i in range(1, len(lst)):
        if lst[i] != lst[write - 1]:
            lst[write] = lst[i]
            write += 1

    del lst[write:]

def interleaving_indices(
    a: Sequence[T],
    b: Sequence[T],
    lessthaneq: Callable[[T, T], bool],
) -> List[int]:
    """
    Determine an interleaving between two sequences.

    Returns a list of 0s and 1s indicating which input sequence
    each successive element should be taken from.

    Args:
        a: First input sequence
        b: Second input sequence
        lessthaneq: Function taking (a_elem, b_elem) and returning
               True if a_elem should come before b_elem

    Returns:
        A list of integers (0 or 1) indicating the interleaving.
    """
    i = j = 0
    result = []

    while i < len(a) and j < len(b):
        if lessthaneq(a[i], b[j]):
            result.append(0)
            i += 1
        else:
            result.append(1)
            j += 1

    # Append remaining elements
    result.extend([0] * (len(a) - i))
    result.extend([1] * (len(b) - j))

    return result
  
def interleaving_overlap(
    a: Sequence[T],
    b: Sequence[T],
    interleaving: List[int],
    lessthaneq: Callable[[T, T], bool],
) -> List[int]:
  i = j = 0
  if interleaving[0] == 0:
    current = a[0]
    i += 1
  else:
    current = b[0]
    j += 1
      
  results = [1]
  for z in interleaving[1:]:
    if z == 0:
      if lessthaneq(current, a[i]) and lessthaneq(a[i], current):
        results.append(0)
      else:
        results.append(1)
      current = a[i]
      i += 1
    else:
      if lessthaneq(current, b[j]) and lessthaneq(b[j], current):
        results.append(0)
      else:
        results.append(1)
      current = b[j]
      j += 1
      
  return results

edge = Tuple[str, str]

class kmer_lessthaneq:
  _ALPHA = {"N": 0, "A": 1, "C": 2, "G": 3, "T": 4}
  
  def __init__(self, h: int) -> None:
     self.h = h
          
  def __call__(self, a: str, b: str) -> bool:
    return [self._ALPHA[c] for c in a[:-(self.h+1):-1]] <= [self._ALPHA[c] for c in b[:-(self.h+1):-1]]

def merge(g1: succinct_graph_container, g2: succinct_graph_container, h: int):
  g1nodes = [g1.edge(i)[0] for i in range(g1.num_edges)]
  g2nodes = [g2.edge(i)[0] for i in range(g2.num_edges)]

  # Make unique
  remove_consecutive_duplicates(g1nodes)
  remove_consecutive_duplicates(g2nodes)
  
  # Do interleaving
  Z = interleaving_indices(g1nodes, g2nodes, kmer_lessthaneq(h))
  B = interleaving_overlap(g1nodes, g2nodes, Z, kmer_lessthaneq(h))
  
  # Print results
  x = y = 0
  for i, b in zip(Z, B):
    if i == 0:
      print(i, g1nodes[x], b)
      x += 1
    else:
      print(i, g2nodes[y], b)
      y += 1
  
def prepare_sequences(k: int, *seqs: str):
  for seq in seqs:
    for seg in seq.split("N"):
      # Validate arguments
      if len(seg) >= k:
        yield seg
        yield reverse_complement(seg)
  
def main(seq1: str, seq2: str, k: int, m: int):
  assert m <= k, "m greater than k!"
  
  # Create graphs
  g1 = construct_sdbg(list(prepare_sequences(k, seq1)), k)
  g2 = construct_sdbg(list(prepare_sequences(k, seq2)), k)
  
  # Merge
  merge(g1, g2, m)
  
def parse_args():
  parser = ArgumentParser()
  parser.add_argument("seq1", type=str)
  parser.add_argument("seq2", type=str)
  parser.add_argument("-k", type=int, required=True)
  parser.add_argument("-m", type=int, required=True)
  return parser.parse_args()
  
if __name__ == "__main__":
  args = parse_args()
  main(args.seq1, args.seq2, args.k, args.m)

