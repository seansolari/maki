
"""
Test 5: K-mer Prevalence control test.

Generate N identical genomes and run the kmer_prevalence pipeline.
Assert signature_prevalence == node_size for all rows.
"""
import os
import pytest

from helpers import genomes, phylogeny, io_utils, pipeline_runner

DEFAULT_L = int(os.getenv("MAKI_TEST_L", "10000"))
DEFAULT_N = int(os.getenv("MAKI_TEST_N", "3"))
DEFAULT_KMER = int(os.getenv("MAKI_TEST_KMER", "31"))

@pytest.mark.parametrize("L,N", [(DEFAULT_L, DEFAULT_N)])
def test_prevalence_control(tmp_path, rng, L, N):
    template = genomes.make_random_genome(L, rng)
    fasta_paths = []
    genome_names = []
    for i in range(N):
        name = f"genome_{i+1:03d}"
        fasta = tmp_path / f"{name}.fa"
        # NOTE: corrected filename (remove stray brace)
        fasta = tmp_path / f"{name}.fa"
        genomes.write_fasta(fasta, name, template)
        fasta_paths.append(fasta)
        genome_names.append(name)

    manifest = tmp_path / "genomes.txt"
    io_utils.write_manifest(manifest, fasta_paths)

    tree = tmp_path / "tree.nwk"
    phylogeny.write_newick(tree, phylogeny.make_fork_newick(genome_names))

    params = {
        "genomes": str(manifest),
        "phylogeny": str(tree),
        "kmer_size": DEFAULT_KMER,
        "outdir": "results",
    }
    run = pipeline_runner.run_pipeline(
        entry="kmer_prevalence", params=params, workdir=tmp_path, profile=os.getenv("MAKI_TEST_PROFILE", "conda")
    )
    assert run.returncode == 0, f"Pipeline failed:\nSTDERR:\n{run.stderr}\nSTDOUT:\n{run.stdout}"

    pipeline_runner.expect_output([run.expected_outputs["cooc_summary"]], base=tmp_path / "results")
    published = list((tmp_path / "results").rglob(run.expected_outputs["cooc_summary"].name))
    assert published, "Co-occ summary output not found under results/**"
    rows = io_utils.read_tsv_gz(published[0])

    for row in rows:
        signature_prev, node_size, _, _ = io_utils.parse_prevalence_row(row)
        assert abs(signature_prev - float(node_size)) < 1e-9, (
            f"signature_prevalence {signature_prev} != node_size {node_size}"
        )
