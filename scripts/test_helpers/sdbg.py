#!/usr/bin/env python3
from __future__ import annotations
import random
from typing import Iterable, List, Sequence

_NI = ord("N")
_AI = ord("A")
_CI = ord("C")
_GI = ord("G")
_TI = ord("T")

def vadd_modn(arr: List[int], l: int, r: int, delta: int, n: int):
    for i in range(l, r):
        result = arr[i] + delta
        if result >= n or result < 0:
            raise TypeError("invalid value!")
        arr[i] = result

def zeros(n: int):
    return [0] * n

def make_dna_table(n: int = 0, a: int = 0, c: int = 0, g: int = 0, t: int = 0):
    return [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, a, 0, c, 0, 0, 0, g, 0, 0, 0, 0, 0, 0, n, 0,
        0, 0, 0, 0, t, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    ]

ORDER = make_dna_table(0, 1, 2, 3, 4)

def isort(arr: str, order: List[int], o: int, l: int, r: int) -> int:
    """perform insertion sort on unique values in `arr[order[l..r]]`,
    with the resulting ordered pointers written to `order[o..]`, for `o <= l`.
    Returns the number of unique values written.
    """
    # base case
    order[o] = order[l]
    written = 1
    # iterative case
    for tmp in order[l+1:r]:
        for i in range(o,o+written):
            if arr[order[i]] == arr[tmp]:
                break
            elif ORDER[ord(arr[order[i]])] > ORDER[ord(arr[tmp])]:
                tmp, order[i] = order[i], tmp
        else:
            order[o+written] = tmp
            written += 1
    return written

def resize(n: int, *arrs: List[int]):
    for arr in arrs:
        yield arr[:n]

# radix sort

def msb_radix_sort(seq: str, order: List[int], l: int, r: int, rem: int, block: List[int], bval: int):
    # sort range
    ni = i = l
    ti = r
    while i < ti:
        if ord(seq[order[i]]) == _NI:
            order[i], order[ni] = order[ni], order[i]
            i += 1
            ni += 1
        elif ord(seq[order[i]]) == _TI:
            ti -= 1
            order[i], order[ti] = order[ti], order[i]
        else:
            i += 1
    ai = i = ni
    gi = ti
    while i < gi:
        if ord(seq[order[i]]) == _AI:
            order[i], order[ai] = order[ai], order[i]
            i += 1
            ai += 1
        elif ord(seq[order[i]]) == _GI:
            gi -= 1
            order[i], order[gi] = order[gi], order[i]
        else:
            i += 1
    # compare prefixes
    prv = None
    for i in range(l, r):
        x = seq[order[i]]
        if (prv is None or x != prv) and block[i] == 0:
            block[i] = bval
        prv = x
    # recursive step
    if rem > 0:
        vadd_modn(order, l, r, -1, len(seq))
        msb_radix_sort(seq, order,  l, ni, rem-1, block, bval+1)
        msb_radix_sort(seq, order, ni, ai, rem-1, block, bval+1)
        msb_radix_sort(seq, order, ai, gi, rem-1, block, bval+1)
        msb_radix_sort(seq, order, gi, ti, rem-1, block, bval+1)
        msb_radix_sort(seq, order, ti,  r, rem-1, block, bval+1)

# k-mer extraction

def radix_extract_kmers(seq: str, k: int, order: List[int]):
    """extract k-mers using a in-place most-significant-bit Radix sort and return
    B array as defined in Egidi et al.
    """
    n = len(order)
    b = zeros(n)
    # recursive case - sort k-mers
    msb_radix_sort(seq, order, 0, n, k-1, b, 1)
    vadd_modn(order, 0, n, k, len(seq))
    return order, b

# nucleotide counter

class nt_tracker:
    def __init__(self) -> None:
        self._table = make_dna_table(False, False, False, False, False)
    
    def set(self, nt: str):
        self._table[ord(nt)] = True

    def is_set(self, nt: str):
        return self._table[ord(nt)]
    
    def reset(self):
        self._table = make_dna_table(False, False, False, False, False)

class nt_counter:
    def __init__(self) -> None:
        self._table = make_dna_table()

    def add(self, nt: str):
        self._table[ord(nt)] += 1

    def to_list(self):
        # accumulate values
        NX = self._table[_NI]
        AX = self._table[_AI]
        CX = self._table[_CI]
        GX = self._table[_GI]
        return [NX, NX+AX, NX+AX+CX, NX+AX+CX+GX]

# graph construction

class succinct_graph_container:
    def __init__(self, k: int, seq: str, order: Sequence[int], blocks: Sequence[int], wminus: Sequence[int], wplus: Sequence[int], F: Sequence[int]) -> None:
        self._k = k
        self._seq = seq
        self._order = order
        self._blocks = blocks
        self._w = [seq[i] for i in order]
        self._wminus = wminus
        self._wplus = wplus
        assert len(F) == 4
        self._F = F

    def __iter__(self):
        return zip(self._order, self._blocks, self._w, self._wminus, self._wplus)
    
    @property
    def seq(self):
        return self._seq
    
    @property
    def kmer_size(self) -> int:
        return self._k

    @property
    def num_nodes(self) -> int:
        """count the number of nodes in the de Bruijn graph (excluding
        the initial `$..$` node).
        """
        return sum(self._wminus)

    @property
    def num_edges(self) -> int:
        return len(self._w)
    
    @property
    def W(self):
        return self._w

    @property
    def wminus(self):
        return self._wminus
    
    @property
    def wplus(self):
        return self._wplus

    @property
    def F(self):
        return self._F
    
    def count_node_last_characters(self) -> List[int]:
        tbl = make_dna_table()
        for e, wm in zip(self.W, self.wminus):
            if wm == 1:
                tbl[ord(e)] += 1
        return [tbl[_AI], tbl[_CI], tbl[_GI], tbl[_TI]]
    
    def find_first_node(self) -> int:
        for i, wp in enumerate(self._wplus):
            if wp == 1:
                return i
        else:
            raise ValueError("no nodes found!")
          
    def edge(self, x: int):
        i = self._order[x]
        suff_i = max(0, i-self._k)
        suff = self.seq[suff_i:i]
        if len(suff) < self._k:
              df = self._k - len(suff)
              suff = self.seq[-df:] + suff
        return suff, self.seq[i]

def finalise_graph(seq: str, k: int, order: List[int], blocks: List[int]):
    NE = len(order)
    F = nt_counter()
    assert NE == len(order)
    assert NE == len(blocks)
    # initialise wplus array data
    wplus = zeros(NE)
    # initialise wminus array data
    wminus = zeros(NE)
    nts = nt_tracker()
    # begin edge sorting
    l = p = 0
    while l < NE:
        # find block of edges corresponding to the same graph node
        r = l + 1
        while r < NE and blocks[r] == 0:
            r += 1
        # sort edges and keep unique (k+1)-mers
        delta = isort(seq, order, p, l, r)
        # set last outgoing edge
        wplus[p+delta-1] = 1
        # check for first incoming edge per node
        if blocks[l] < k:  # by construction, blocks[l] > 0
            nts.reset()
        for e in range(p, p+delta):
            # count edges for F array
            F.add(seq[order[e]-1])
            # check last occurrence for wpos
            if not nts.is_set(seq[order[e]]):
                wminus[e] = 1
                nts.set(seq[order[e]])
        # next block
        l = r
        p += delta
    # resize to unique values
    order, wminus, wplus = resize(p, order, wminus, wplus)
    return succinct_graph_container(k, seq, order, blocks, wminus, wplus, F.to_list())

def construct_sdbg(sequences: Iterable[str], k: int):
    # construct a concatenated sequence from the inputs
    buffer, buffer_len = "", 0
    order = []
    for sequence in sequences:
        # add valid start indices
        seq_len = len(sequence)
        order.extend(range(buffer_len+k-1,buffer_len+k+seq_len))
        # append sequence to buffer
        buffer += (k*"N") + sequence
        buffer_len += k + seq_len
    buffer += "N"
    # extract k-mers using radix sort
    order, block = radix_extract_kmers(buffer, k, order)
    # sort on edges and construct graph arrays
    graph = finalise_graph(buffer, k, order, block)
    return graph

# testing

def print_kmers(st: str, k: int, p: Sequence[int], b: Sequence[int]):
    for r, h in zip(p, b):
        suff_len = k if r > k else r
        suff = st[r-suff_len:r]
        pre = " "*(k-suff_len)
        print(pre + suff + " %s\t%d\t%d" % (st[r], r, h))

def print_graph(g: succinct_graph_container):
    print("<<< BEGIN GRAPH >>>")
    print("Node\tEdge\tW(-)\tW(+)")
    for i, b, w, wmin, wpls in g:
        suff_i = max(0, i-g._k)
        suff = g.seq[suff_i:i]
        if len(suff) < g._k:
            df = g._k - len(suff)
            suff = g.seq[-df:] + suff
        print(suff + "\t%s\t%d\t%d" % (g.seq[i], wmin, wpls))
    print("\nF Array: ", ", ".join(str(v) for v in g._F))
    print("<<< END GRAPH >>>")

def random_graph_test(l: int = 50, n: int = 3, k: int = 5):
    import random
    g = construct_sdbg(["".join(random.choices("ACGT", k=l)) for _ in range(n)], k)
    print("STRING: ", g.seq)
    print_graph(g)

def known_graph_test():
    k = 3
    g = construct_sdbg(["TACACT", "TACTCG", "GACTCA"], k)
    print("STRING: ", g.seq)
    print_graph(g)

# succinct graph interface

class succinct_graph:
    def __init__(self, container: succinct_graph_container) -> None:
        self._k = container._k
        self._w = container._w
        self._wminus = container._wminus
        self._wplus = container._wplus
        self._F = container._F

    @property
    def kmer_size(self) -> int:
        return self._k

    @property
    def num_nodes(self) -> int:
        """count the number of nodes in the de Bruijn graph (excluding
        the initial `$..$` node).
        """
        return sum(self._wminus)

    @property
    def num_edges(self) -> int:
        return len(self._w)
    
    @property
    def W(self):
        return self._w

    @property
    def wminus(self):
        return self._wminus
    
    @property
    def wplus(self):
        return self._wplus

    @property
    def F(self):
        return self._F
    
    def count_node_last_characters(self) -> List[int]:
        tbl = make_dna_table()
        for e, wm in zip(self.W, self.wminus):
            if wm == 1:
                tbl[ord(e)] += 1
        return [tbl[_AI], tbl[_CI], tbl[_GI], tbl[_GI]]

def to_graph(sequences: Iterable[str], k: int):
    return succinct_graph(construct_sdbg(sequences, k))

# testing

def main():
    k = 3
    reference = construct_sdbg(["TACACT", "TACTCG", "AGTGTA", "CGAGTA"], k)
    print_graph(reference)
    query = construct_sdbg(["GACTCA", "TGAGTC"], k)
    print_graph(query)
    other = construct_sdbg(["ACGGTCGG", "CCGACCGT"], k)
    print_graph(other)

# main

if __name__ == "__main__":
    random.seed(1)
    print("(succinct graph testing module)")
    main()
