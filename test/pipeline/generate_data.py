
#!/usr/bin/env python3
"""
Small command-line application to generate synthetic datasets that match the
pytest helpers and Nextflow pipelines.

It interfaces directly with the helpers you added under `tests/helpers/` and
produces files such as FASTA genomes, GFF3+FASTA annotated genomes, manifests,
phylogeny trees, and cluster maps.

Usage examples:
  # Jaccard control (identical genomes)
  python maki_data_cli.py jaccard-control --L 10000 --N 3 --outdir ./data/jaccard_ctl

  # Jaccard mutation (per-base rate R)
  python maki_data_cli.py jaccard-mutation --L 10000 --N 3 --R 0.05 --outdir ./data/jaccard_mut

  # Random Gene Distance control (annotated identical genomes)
  python maki_data_cli.py randist-control --L 10000 --N 3 --A 300 --outdir ./data/randist_ctl

  # Random Gene Distance mutation (per-gene rates in [0.05,0.6])
  python maki_data_cli.py randist-mutation --L 10000 --N 3 --A 300 --outdir ./data/randist_mut

  # K-mer prevalence control/random datasets (FASTA only)
  python maki_data_cli.py prevalence-control --L 10000 --N 3 --outdir ./data/prev_ctl
  python maki_data_cli.py prevalence-random --L 10000 --N 3 --R 0.05 --outdir ./data/prev_rand

Environment:
- Tries to import helpers from `tests/helpers/`. If not found, it modifies
  `sys.path` to locate a local `helpers/` directory.
- Random seed can be controlled via `--seed` (default: 42), using NumPy if
  available, otherwise stdlib `random.Random`.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import sys

# Resolve helpers (prefer tests/helpers)
try:
    from tests.helpers import genomes, gff, phylogeny, clusters, io_utils
except Exception:
    sys.path.insert(0, str(Path(__file__).parent / 'tests'))
    sys.path.insert(0, str(Path(__file__).parent))
    from helpers import genomes, gff, phylogeny, clusters, io_utils  # type: ignore


def make_rng(seed: int):
    """Create a reproducible RNG (NumPy if available, else stdlib)."""
    try:
        import numpy as np
        return np.random.default_rng(seed)
    except Exception:
        import random
        return random.Random(seed)


def ensure_dir(d: Path) -> Path:
    d.mkdir(parents=True, exist_ok=True)
    return d


def write_manifest_and_tree(outdir: Path, fasta_paths: list[Path], genome_names: list[str]) -> tuple[Path, Path]:
    manifest = outdir / 'genomes.txt'
    io_utils.write_manifest(manifest, fasta_paths)
    newick = phylogeny.make_fork_newick(genome_names)
    tree = outdir / 'tree.nwk'
    phylogeny.write_newick(tree, newick)
    return manifest, tree


def cmd_jaccard_control(L: int, N: int, outdir: Path, seed: int):
    rng = make_rng(seed)
    ensure_dir(outdir)
    genomes_dir = ensure_dir(outdir / 'genomes')

    template = genomes.make_random_genome(L, rng)
    fasta_paths, genome_names = [], []
    for i in range(N):
        name = f'genome_{i+1:03d}'
        p = genomes_dir / f'{name}.fa'
        genomes.write_fasta(p, name, template)
        fasta_paths.append(p)
        genome_names.append(name)

    manifest, tree = write_manifest_and_tree(outdir, fasta_paths, genome_names)
    print(f'Wrote manifest: {manifest}')
    print(f'Wrote phylogeny: {tree}')
    print(f'Genomes dir: {genomes_dir}')


def cmd_jaccard_mutation(L: int, N: int, R: float, outdir: Path, seed: int):
    rng = make_rng(seed)
    ensure_dir(outdir)
    genomes_dir = ensure_dir(outdir / 'genomes')

    template = genomes.make_random_genome(L, rng)
    fasta_paths, genome_names = [], []
    for i in range(N):
        name = f'genome_{i+1:03d}'
        mutated = genomes.mutate_genome(template, R, rng)
        p = genomes_dir / f'{name}.fa'
        genomes.write_fasta(p, name, mutated)
        fasta_paths.append(p)
        genome_names.append(name)

    manifest, tree = write_manifest_and_tree(outdir, fasta_paths, genome_names)
    print(f'Wrote manifest: {manifest}')
    print(f'Wrote phylogeny: {tree}')
    print(f'Genomes dir: {genomes_dir}')


def _build_annotations(L: int, A: float, rng):
    return gff.build_nonoverlapping_annotations(L, A, rng)


def cmd_randist_control(L: int, N: int, A: float, outdir: Path, seed: int):
    rng = make_rng(seed)
    ensure_dir(outdir)
    ann_dir = ensure_dir(outdir / 'annotated')

    template = genomes.make_random_genome(L, rng)
    annotations = _build_annotations(L, A, rng)

    gff_paths, genome_names = [], []
    accession_map: dict[str, str] = {}
    for i in range(N):
        name = f'genome_{i+1:03d}'
        accession_map[name] = name
        p = ann_dir / f'{name}.gff3'
        gff.write_gff3_with_fasta(p, accession=name, seq=template, annotations=annotations)
        gff_paths.append(p)
        genome_names.append(name)

    # manifest (paths to GFF3) and fork tree
    manifest, tree = write_manifest_and_tree(outdir, gff_paths, genome_names)

    # clusters: representative = first genome
    gene_ids = [ann.gene_id for ann in annotations]
    cluster_file = outdir / 'clusters.tsv'
    clusters.write_cluster_map(cluster_file, representative_genome=genome_names[0], accession_map=accession_map, gene_ids=gene_ids)

    print(f'Wrote manifest: {manifest}')
    print(f'Wrote phylogeny: {tree}')
    print(f'Wrote clusters: {cluster_file}')
    print(f'Annotated dir: {ann_dir}')


def cmd_randist_mutation(L: int, N: int, A: float, outdir: Path, seed: int):
    rng = make_rng(seed)
    ensure_dir(outdir)
    ann_dir = ensure_dir(outdir / 'annotated')
    rlow, rhigh = 0.01, 0.1

    template = genomes.make_random_genome(L, rng)
    annotations = _build_annotations(L, A, rng)
    gene_ids = [ann.gene_id for ann in annotations]

    # Per-gene mutation rate in [rlow, rhigh]
    mu_map: dict[str, float] = {}
    # NumPy path
    if hasattr(rng, 'uniform'):
        for gid in gene_ids:
            mu_map[gid] = float(rng.uniform(rlow, rhigh))
    else:
        import random
        rr = rng if hasattr(rng, 'random') else random.Random(seed)
        for gid in gene_ids:
            mu_map[gid] = rlow + (rhigh - rlow) * rr.random()

    gff_paths, genome_names = [], []
    accession_map: dict[str, str] = {}
    for i in range(N):
        name = f'genome_{i+1:03d}'
        accession_map[name] = name
        mutated_seq = gff.mutate_annotated_template_per_segment(template, annotations, mu_map, rng)
        attrs = {gid: {"rate": f"{mu_map[gid]:.6f}"} for gid in gene_ids}
        p = ann_dir / f'{name}.gff3'
        gff.write_gff3_with_fasta(p, accession=name, seq=mutated_seq, annotations=annotations, extra_attrs_for_gene=attrs)
        gff_paths.append(p)
        genome_names.append(name)

    manifest, tree = write_manifest_and_tree(outdir, gff_paths, genome_names)
    cluster_file = outdir / 'clusters.tsv'
    clusters.write_cluster_map(cluster_file, representative_genome=genome_names[0], accession_map=accession_map, gene_ids=gene_ids)

    print(f'Wrote manifest: {manifest}')
    print(f'Wrote phylogeny: {tree}')
    print(f'Wrote clusters: {cluster_file}')
    print(f'Annotated dir: {ann_dir}')


def cmd_prevalence_control(L: int, N: int, outdir: Path, seed: int):
    # Same data form as jaccard-control (FASTA + manifest + tree)
    return cmd_jaccard_control(L=L, N=N, outdir=outdir, seed=seed)


def cmd_prevalence_random(L: int, N: int, R: float, outdir: Path, seed: int):
    # Same data form as jaccard-mutation (FASTA + manifest + tree)
    return cmd_jaccard_mutation(L=L, N=N, R=R, outdir=outdir, seed=seed)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Synthetic dataset generator for maki Nextflow pipelines")
    sub = p.add_subparsers(dest='cmd', required=True)

    # Shared defaults
    DEFAULT_L = 10000
    DEFAULT_N = 3
    DEFAULT_A = 300.0
    DEFAULT_R = 0.05

    # jaccard-control
    jc = sub.add_parser('jaccard-control', help='FASTA identical genomes for Jaccard pipeline')
    jc.add_argument('--L', type=int, default=DEFAULT_L)
    jc.add_argument('--N', type=int, default=DEFAULT_N)
    jc.add_argument('--outdir', type=Path, required=True)
    jc.add_argument('--seed', type=int, default=42)

    # jaccard-mutation
    jm = sub.add_parser('jaccard-mutation', help='FASTA mutated genomes for Jaccard pipeline')
    jm.add_argument('--L', type=int, default=DEFAULT_L)
    jm.add_argument('--N', type=int, default=DEFAULT_N)
    jm.add_argument('--R', type=float, default=DEFAULT_R)
    jm.add_argument('--outdir', type=Path, required=True)
    jm.add_argument('--seed', type=int, default=42)

    # randist-control
    rc = sub.add_parser('randist-control', help='GFF3 identical annotated genomes for Random Gene Distance pipeline')
    rc.add_argument('--L', type=int, default=DEFAULT_L)
    rc.add_argument('--N', type=int, default=DEFAULT_N)
    rc.add_argument('--A', type=float, default=DEFAULT_A)
    rc.add_argument('--outdir', type=Path, required=True)
    rc.add_argument('--seed', type=int, default=42)

    # randist-mutation
    rm = sub.add_parser('randist-mutation', help='GFF3 mutated annotated genomes for Random Gene Distance pipeline')
    rm.add_argument('--L', type=int, default=DEFAULT_L)
    rm.add_argument('--N', type=int, default=DEFAULT_N)
    rm.add_argument('--A', type=float, default=DEFAULT_A)
    rm.add_argument('--outdir', type=Path, required=True)
    rm.add_argument('--seed', type=int, default=42)

    # prevalence-control
    pc = sub.add_parser('prevalence-control', help='FASTA identical genomes for K-mer Prevalence pipeline')
    pc.add_argument('--L', type=int, default=DEFAULT_L)
    pc.add_argument('--N', type=int, default=DEFAULT_N)
    pc.add_argument('--outdir', type=Path, required=True)
    pc.add_argument('--seed', type=int, default=42)

    # prevalence-random
    pr = sub.add_parser('prevalence-random', help='FASTA mutated genomes for K-mer Prevalence pipeline')
    pr.add_argument('--L', type=int, default=DEFAULT_L)
    pr.add_argument('--N', type=int, default=DEFAULT_N)
    pr.add_argument('--R', type=float, default=DEFAULT_R)
    pr.add_argument('--outdir', type=Path, required=True)
    pr.add_argument('--seed', type=int, default=42)

    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    cmd = args.cmd
    if cmd == 'jaccard-control':
        cmd_jaccard_control(L=args.L, N=args.N, outdir=args.outdir, seed=args.seed)
    elif cmd == 'jaccard-mutation':
        cmd_jaccard_mutation(L=args.L, N=args.N, R=args.R, outdir=args.outdir, seed=args.seed)
    elif cmd == 'randist-control':
        cmd_randist_control(L=args.L, N=args.N, A=args.A, outdir=args.outdir, seed=args.seed)
    elif cmd == 'randist-mutation':
        cmd_randist_mutation(L=args.L, N=args.N, A=args.A, outdir=args.outdir, seed=args.seed)
    elif cmd == 'prevalence-control':
        cmd_prevalence_control(L=args.L, N=args.N, outdir=args.outdir, seed=args.seed)
    elif cmd == 'prevalence-random':
        cmd_prevalence_random(L=args.L, N=args.N, R=args.R, outdir=args.outdir, seed=args.seed)
    else:
        raise SystemExit(f'Unknown command: {cmd}')

if __name__ == '__main__':
    main()
