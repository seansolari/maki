#!/usr/bin/env python3

from argparse import ArgumentParser
import sys
import timeit
from typing import Callable, List
import numpy as np
from parallelbar import progress_starmap


ALPHABET = np.array([0, 1, 2, 3], dtype=np.uint8)

class Kmers:
    def __init__(self, kmers: np.ndarray) -> None:
        self.kmers = kmers

    def print(self):
        print(self.kmers)

    def size(self) -> int:
        return self.kmers.shape[0]
    
    def kmer_size(self) -> int:
        return self.kmers.shape[1]
    
class RandomKmers(Kmers):
    def __init__(self, num_kmers: int, k: int, seed: int) -> None:
        rng = np.random.default_rng(seed=seed)
        super().__init__(np.unique(rng.choice(ALPHABET, (num_kmers, k)), axis=0))
    
def union_using_unique(a: Kmers, b: Kmers) -> int:
    ab = np.vstack((a.kmers, b.kmers))
    uq = np.unique(ab, axis=0)
    return uq.shape[0]

# https://stackoverflow.com/a/42937116
def intersect_using_search(a: Kmers, b: Kmers) -> int:
    dims = np.repeat(4, a.kmers.shape[1])

    a1D = np.ravel_multi_index(a.kmers.T, dims)
    b1D = np.ravel_multi_index(b.kmers.T, dims)

    fidx = np.searchsorted(a1D, b1D)
    fidx[fidx == a1D.size] = 0
    return np.sum(a1D[fidx] == b1D)

def union_using_search(a: Kmers, b: Kmers):
    ix = intersect_using_search(a, b) if a.size() < b.size() else intersect_using_search(b, a)
    return a.size() + b.size() - ix

class JaccardResult:
    def __init__(self, a: Kmers, b: Kmers) -> None:
        assert a.kmer_size() == b.kmer_size()

        U = union_using_search(a, b)

        self.A = a.size()
        self.B = b.size()
        self.AB = U
        self.J = (self.A + self.B - U) / U
        self.k = a.kmer_size()

    def print(self):
        print("%d\t%d\t%d\t%f" % (self.A, self.B, self.AB, self.J))

class RandomSequence:
    def __init__(self, seq: np.ndarray) -> None:
        self.seq = seq

    def unique_kmers(self, k: int):
        return Kmers(np.unique(np.lib.stride_tricks.sliding_window_view(self.seq, k), axis=0))

class RandomSequenceConfiguration:
    def __init__(self, length: int, seed: int) -> None:
        self.length = length
        self.seed = seed
        self.rng = np.random.default_rng(seed=seed)

    def generate(self) -> RandomSequence:
        return RandomSequence(self.rng.choice(ALPHABET, size=self.length))

def lower_triangular(n: int):
    for i in range(1, n):
        for j in range(i):
            yield i, j

def random_sequence_experiment(l1: int, s1: int, l2: int, s2: int, kmer_sizes: List[int]) -> List[JaccardResult]:
    seq1 = RandomSequenceConfiguration(l1, s1).generate()
    seq2 = RandomSequenceConfiguration(l2, s2).generate()
    return [
        JaccardResult(seq1.unique_kmers(k), seq2.unique_kmers(k))
        for k in kmer_sizes
    ]

def random_kmers_experiment(l1: int, s1: int, l2: int, s2: int, kmer_sizes: List[int]) -> List[JaccardResult]:
    return [
        JaccardResult(RandomKmers(l1 - k + 1, k, s1), RandomKmers(l2 - k + 1, k, s2))
        for k in kmer_sizes
    ]

ExperimentFunction = Callable[[int, int, int, int, List[int]], List[JaccardResult]]

def run_benchmark():
    # generate kmer data
    sequence_length = 100000
    k = 25
    a = RandomSequenceConfiguration(sequence_length, 1).generate().unique_kmers(k)
    b = RandomSequenceConfiguration(sequence_length, 2).generate().unique_kmers(k)
    
    assert union_using_unique(a, b) == union_using_search(a, b)

    # METHOD: unique
    uq_time = timeit.timeit(lambda: union_using_unique(a, b), number=10)

    # METHOD: search
    sr_time = timeit.timeit(lambda: union_using_search(a, b), number=10)

    print("===RESULTS===")
    print("METHOD: Unique... %f" % uq_time)
    print("METHOD: Binary Search... %f" % sr_time)

def run_experiments(fn: ExperimentFunction, sequence_lengths: List[int], sequence_replicates: int, kmer_sizes: List[int], seed_start: int = 0):
    random_sequence_parameters = [
        (l, seed_start + (i * sequence_replicates) + j)
        for i, l in enumerate(sequence_lengths)
        for j in range(sequence_replicates)
    ]
    nseqs = len(random_sequence_parameters)
    params = (
        (*random_sequence_parameters[i], *random_sequence_parameters[j], kmer_sizes)
        for i, j in lower_triangular(len(random_sequence_parameters))
    )
    jaccards = [
        j
        for proc_output in progress_starmap(fn, params, total=nseqs * (nseqs - 1) / 2)
        for j in proc_output
    ]
    return jaccards

def main(sequence_lengths: List[int], sequence_replicates: int, kmer_sizes: List[int], outfile: str):
    seq_results = run_experiments(random_sequence_experiment, sequence_lengths, sequence_replicates, kmer_sizes, seed_start=0)
    kmer_results = run_experiments(random_kmers_experiment, sequence_lengths, sequence_replicates, kmer_sizes, seed_start=len(seq_results))
    with open(outfile, 'w') as f:
        for j in seq_results:
            f.write(f"{j.A:d}\t{j.B:d}\t{j.AB:d}\t{j.J:e}\t{j.k:d}\tseq\n")
        for j in kmer_results:
            f.write(f"{j.A:d}\t{j.B:d}\t{j.AB:d}\t{j.J:e}\t{j.k:d}\tkmer\n")

def parse_args():
    parser = ArgumentParser()
    parser.add_argument("-l", action="append", dest="l", type=int, required=True, help="Sequence length")
    parser.add_argument("--reps", dest="reps", type=int, required=True, help="Replicates per sequence length")
    parser.add_argument("-k", action="append", dest="k", type=int, required=True, help="K-mer size")
    parser.add_argument("-o", "--out", dest="out", type=str, required=True, help="Output file")
    return parser.parse_args()

if __name__ == "__main__":
    if sys.argv[1] == "benchmark":
        run_benchmark()
    else:
        args = parse_args()
        main(args.l, args.reps, args.k, args.out)
