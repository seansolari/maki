
"""Genome generation and mutation helpers.

This module provides small, dependency-light helpers to:
- generate random DNA genomes of a specified length; and
- mutate genomes using an efficient Poisson-based scheme.

Prefer NumPy if available (for speed); fall back to the Python stdlib.
All randomness must be controlled by an RNG object supplied by callers.
"""
from __future__ import annotations
from pathlib import Path
from typing import Protocol, Sequence

class RNG(Protocol):
    """Protocol for random generators used in tests."""
    # NumPy's Generator implements choice and poisson; stdlib Random does not.
    def choice(self, a, size=None, replace=True): ...

def _has_numpy_rng(rng) -> bool:
    """Return True if `rng` looks like a NumPy Generator (has poisson/choice)."""
    return hasattr(rng, "poisson") and hasattr(rng, "choice")

def make_random_genome(length: int, rng: RNG) -> str:
    """Create a random DNA genome of `length` with uniform bases A/C/G/T."""
    bases: Sequence[str] = ("A", "C", "G", "T")
    if _has_numpy_rng(rng):
        # Vectorized sampling
        seq_arr = rng.choice(bases, size=length, replace=True)
        return "".join(seq_arr.tolist())
    else:
        # Stdlib fallback
        import random
        # If rng is stdlib Random, use it; else fallback to random
        r = rng if isinstance(rng, random.Random) else random
        return "".join(r.choice(bases) for _ in range(length))

def _mutate_base(b: str, rng) -> str:
    """Mutate base `b` to one of the other three nucleotides uniformly."""
    if b not in ("A", "C", "G", "T"):
        # Robustness for unexpected symbols; treat as A
        b = "A"
    others = [x for x in ("A", "C", "G", "T") if x != b]
    if _has_numpy_rng(rng):
        return rng.choice(others)
    else:
        import random
        r = rng if isinstance(rng, random.Random) else random
        return r.choice(others)

def mutate_genome(template: str, mu_rate: float, rng: RNG) -> str:
    """Mutate `template` using Poisson-efficient scheme.

    Steps:
      1) Sample k ~ Poisson(mu_rate * L) (NumPy) or Binomial fallback (stdlib).
      2) Choose k unique positions uniformly.
      3) Mutate each chosen base to one of the other three nucleotides.
    """
    L = len(template)
    if L == 0 or mu_rate <= 0.0:
        return template

    if _has_numpy_rng(rng):
        # NumPy path
        k = int(rng.poisson(mu_rate * L))
        if k <= 0:
            return template
        k = min(k, L)
        positions = rng.choice(L, size=k, replace=False).tolist()
    else:
        # Stdlib fallback — sample k via Binomial(L, mu_rate)
        import random
        r = rng if isinstance(rng, random.Random) else random
        # sum of Bernoulli trials; O(L) but fine for test sizes
        k = sum(1 for _ in range(L) if r.random() < mu_rate)
        if k <= 0:
            return template
        positions = r.sample(range(L), k)

    seq_list = list(template)
    for pos in positions:
        seq_list[pos] = _mutate_base(seq_list[pos], rng)
    return "".join(seq_list)

def write_fasta(path: Path, name: str, seq: str) -> None:
    """Write a single-record FASTA file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        # No line wrap for simplicity
        fh.write(seq + "\n")
