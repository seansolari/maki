
"""
Test 1: Jaccard control test.

Generate N identical genomes (length L) and a fork phylogeny.
Run the Jaccard pipeline and assert jaccard==1.0 for all pairs.
"""
import os
import pytest

from helpers import genomes, phylogeny, io_utils, pipeline_runner

DEFAULT_L = int(os.getenv("MAKI_TEST_L", "10000"))
DEFAULT_N = int(os.getenv("MAKI_TEST_N", "3"))
DEFAULT_KMER = int(os.getenv("MAKI_TEST_KMER", "31"))

@pytest.mark.parametrize("L,N", [(DEFAULT_L, DEFAULT_N)])
def test_jaccard_control(tmp_path, rng, L, N):
    # 1) Generate template genome and N identical copies
    template = genomes.make_random_genome(L, rng)
    fasta_paths = []
    genome_names = []
    for i in range(N):
        name = f"genome_{i+1:03d}"
        fasta = tmp_path / f"{name}.fa"
        genomes.write_fasta(fasta, name, template)
        fasta_paths.append(fasta)
        genome_names.append(name)

    # 2) Manifest and phylogeny
    manifest = tmp_path / "genomes.txt"
    io_utils.write_manifest(manifest, fasta_paths)

    newick = phylogeny.make_fork_newick(genome_names)
    tree = tmp_path / "tree.nwk"
    phylogeny.write_newick(tree, newick)

    # 3) Run pipeline
    params = {
        "genomes": str(manifest),
        "phylogeny": str(tree),
        "kmer_size": DEFAULT_KMER,
        "outdir": "results",
    }
    run = pipeline_runner.run_pipeline(
        entry="jaccard", params=params, workdir=tmp_path, profile=os.getenv("MAKI_TEST_PROFILE", "conda")
    )
    assert run.returncode == 0, f"Pipeline failed:\nSTDERR:\n{run.stderr}\nSTDOUT:\n{run.stdout}"

    # 4) Check outputs exist and parse
    pipeline_runner.expect_output([run.expected_outputs["jaccard"]], base=tmp_path / "results")
    published = list((tmp_path / "results").rglob(run.expected_outputs["jaccard"].name))
    assert published, "Jaccard output not found under results/**"
    rows = io_utils.read_tsv_gz(published[0])

    # 5) Assertions: jaccard==1.0 for all pairs
    for row in rows:
        rec = io_utils.parse_jaccard_row(row)
        assert abs(rec.jaccard - 1.0) < 1e-9, f"Expected 1.0, got {rec.jaccard}"
        assert rec.lca, "LCA should be non-empty"
