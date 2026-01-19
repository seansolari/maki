
"""
Test 4: Random Gene Distance mutation process.

Generate N annotated genomes where each gene is mutated with its own rate
uniformly sampled in [0.05, 0.6]. Assert that for rows with same gene sampled,
R 95% CI contains the true gene mutation rate.
"""
import os
import pytest

from helpers import genomes, gff, phylogeny, clusters, io_utils, pipeline_runner

DEFAULT_L = int(os.getenv("MAKI_TEST_L", "10000"))
DEFAULT_N = int(os.getenv("MAKI_TEST_N", "3"))
DEFAULT_A = float(os.getenv("MAKI_TEST_A", "300"))
DEFAULT_KMER = int(os.getenv("MAKI_TEST_KMER", "31"))

@pytest.mark.parametrize("L,N,A", [(DEFAULT_L, DEFAULT_N, DEFAULT_A)])
def test_randist_mutation(tmp_path, rng, L, N, A):
    # Template and annotations
    template = genomes.make_random_genome(L, rng)
    annotations = gff.build_nonoverlapping_annotations(L, A, rng)
    gene_ids = [ann.gene_id for ann in annotations]

    # Sample per-gene mutation rates in [0.05, 0.6]
    has_numpy = hasattr(rng, "uniform")
    import random
    rr = rng if hasattr(rng, "random") else random.Random(42)
    mu_map = {}
    for gid in gene_ids:
        if has_numpy:
            mu = float(rng.uniform(0.05, 0.6))
        else:
            mu = 0.05 + 0.55 * rr.random()
        mu_map[gid] = mu

    # Build N mutated genomes (per segment) with attributes noting the rate
    gff_paths = []
    genome_names = []
    accession_map = {}
    for i in range(N):
        name = f"genome_{i+1:03d}"
        accession_map[name] = name
        mutated_seq = gff.mutate_annotated_template_per_segment(template, annotations, mu_map, rng)
        attrs = {gid: {"rate": f"{mu_map[gid]:.6f}"} for gid in gene_ids}
        gff_path = tmp_path / f"{name}.gff3"
        gff.write_gff3_with_fasta(gff_path, accession=name, seq=mutated_seq, annotations=annotations, extra_attrs_for_gene=attrs)
        gff_paths.append(gff_path)
        genome_names.append(name)

    manifest = tmp_path / "genomes.txt"
    io_utils.write_manifest(manifest, gff_paths)

    tree = tmp_path / "tree.nwk"
    phylogeny.write_newick(tree, phylogeny.make_fork_newick(genome_names))

    cluster_file = tmp_path / "clusters.tsv"
    clusters.write_cluster_map(cluster_file, representative_genome=genome_names[0], accession_map=accession_map, gene_ids=gene_ids)

    params = {
        "genomes": str(manifest),
        "phylogeny": str(tree),
        "cluster_file": str(cluster_file),
        "kmer_size": DEFAULT_KMER,
        "outdir": "results",
    }
    run = pipeline_runner.run_pipeline(
        entry="random_gene_distances", params=params, workdir=tmp_path, profile=os.getenv("MAKI_TEST_PROFILE", "conda")
    )
    assert run.returncode == 0, f"Pipeline failed:\nSTDERR:\n{run.stderr}\nSTDOUT:\n{run.stdout}"

    pipeline_runner.expect_output([run.expected_outputs["feature_table"]], base=tmp_path / "results")
    published = list((tmp_path / "results").rglob(run.expected_outputs["feature_table"].name))
    assert published, "Feature table output not found under results/**"
    rows = io_utils.read_tsv_gz(published[0])
    assert len(rows) == N**2, f"Expected {N**2} rows, found {len(rows)}"

    total_same = 0
    covered = 0
    for row in rows:
        rec = io_utils.parse_feature_row(row)
        if rec.cluster_1 != rec.cluster_2:
            continue
        # cluster_id format: '{ACCESSION}-gene-{GENEID}' (REP accession)
        try:
            gene_id = rec.cluster_1.split('-gene-')[1]
        except Exception:
            continue
        mu_true = mu_map.get(gene_id)
        if mu_true is None:
            continue
        total_same += 1
        if rec.r_lo <= mu_true <= rec.r_hi:
            covered += 1
    assert total_same > 0, "No rows with same gene sampled (cluster_id != 'NA')"
    assert covered / total_same >= 0.95, f"CI coverage {covered}/{total_same} < 95%"
