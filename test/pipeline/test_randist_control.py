
"""
Test 3: Random Gene Distance control test.

Generate N identical annotated genomes (GFF3 + FASTA) with non-overlapping segments.
Run the Random Gene Distance pipeline and assert:
  - rows with cluster_id != 'NA' have jaccard == 1.0
  - output contains a mix of cluster_id values (ideally some 'NA').
"""
import os
import pytest

from helpers import genomes, gff, phylogeny, clusters, io_utils, pipeline_runner

DEFAULT_L = int(os.getenv("MAKI_TEST_L", "10000"))
DEFAULT_N = int(os.getenv("MAKI_TEST_N", "3"))
DEFAULT_A = float(os.getenv("MAKI_TEST_A", "300"))
DEFAULT_KMER = int(os.getenv("MAKI_TEST_KMER", "31"))

@pytest.mark.parametrize("L,N,A", [(DEFAULT_L, DEFAULT_N, DEFAULT_A)])
def test_randist_control(tmp_path, rng, L, N, A):
    # Template genome and annotations
    template = genomes.make_random_genome(L, rng)
    annotations = gff.build_nonoverlapping_annotations(L, A, rng)

    # N copies with unique accessions
    gff_paths = []
    genome_names = []
    accession_map = {}
    for i in range(N):
        name = f"genome_{i+1:03d}"
        accession_map[name] = name  # use genome name as accession
        gff_path = tmp_path / f"{name}.gff3"
        # NOTE: corrected filename (remove stray brace)
        gff_path = tmp_path / f"{name}.gff3"
        gff.write_gff3_with_fasta(gff_path, accession=name, seq=template, annotations=annotations)
        gff_paths.append(gff_path)
        genome_names.append(name)

    # Manifest and fork phylogeny
    manifest = tmp_path / "genomes.txt"
    io_utils.write_manifest(manifest, gff_paths)

    tree = tmp_path / "tree.nwk"
    phylogeny.write_newick(tree, phylogeny.make_fork_newick(genome_names))

    # Cluster file: representative = first genome
    gene_ids = [ann.gene_id for ann in annotations]
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

    seen_na = False
    seen_match = False
    for row in rows:
        rec = io_utils.parse_feature_row(row)
        if rec.cluster_1 != rec.cluster_2:
            seen_na = True
        elif rec.kmer_matches > 0:
            seen_match = True
            assert abs(rec.jaccard - 1.0) < 1e-9, f"Expected 1.0 jaccard for same gene, got {rec.jaccard}"
    # Prefer both categories; if randomness yields only NA, assert at least one row exists
    assert len(rows) > 0, "No rows in feature-table"
    assert seen_match, "No rows with same gene sampled (cluster_id != 'NA')"
