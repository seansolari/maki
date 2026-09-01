from __future__ import annotations
import logging
from pathlib import Path
from typing import Annotated

import typer


logger = logging.getLogger(__name__)


classify_app = typer.Typer()


def _samples_from_cli(input_data: Path):
    from maki.models.analyse.files import (
        discover_paired_fastqs,
        group_fastqs,
        parse_paired_fastq_list
    )
    
    if not input_data.exists():
        logger.error("Input does not exist: %s", input_data)
        return 1
    elif input_data.is_file() and any(ext in input_data.suffixes for ext in (".fastq", ".fq")):
        logger.error("Input file is FASTQ, must be either directory or manifest")
        return 1
    
    files = discover_paired_fastqs(input_data) if input_data.is_dir() else parse_paired_fastq_list(input_data)
    
    if not files:
        logger.warning("No files discovered from %s", input_data)
        return 0
    
    samples = group_fastqs(files)
    return samples


@classify_app.command(help="""
Initialise output folder.
""")
def init(
    output_root: Annotated[Path, typer.Option("-o", "--output", help="Classification output path")],
    kmer_size: Annotated[int, typer.Option("-k", "--kmer-size", help="K-mer size for analysis")],
    scale: Annotated[int, typer.Option(help="Sketching scale parameter"),] = 1000,
    seed: Annotated[int, typer.Option(help="Seed for sketching"),] = 42,
    permissive: Annotated[bool, typer.Option("--permissive", "-p", help="Do not fail if one already exists")] = False
):
    from maki.sketch.core import SketchParameters
    from maki.models.analyse.hook import OutputHook, ClassifyParameters
    
    OutputHook.init(
        output_root,
        ClassifyParameters(),
        SketchParameters(kmer_size, scale, seed),
        exist_ok=permissive
    )
    
    logger.info("Output folder initialised at %s.", output_root)
    
    return 0


@classify_app.command(help="""
Sketch input metagenomes for pre-filtering.
""")
def sketch(
    input_data: Annotated[Path, typer.Option("-i", "--input", help="Sample manifest file or base directory."),],
    output_root: Annotated[Path, typer.Option("-o", "--output", help="Classification output path")],
    workers: Annotated[int, typer.Option("-w", "--workers", help="Number of Sourmash sketching workers to run in parallel.")],
    overwrite: Annotated[bool, typer.Option("--overwrite", help="Overwrite old results if they have a hash mismatch")] = False
):
    samples = _samples_from_cli(input_data)
    
    if isinstance(samples, int):
        return samples
    
    # sketch samples
    from maki.models.analyse.hook import OutputHook
    
    hook = OutputHook(output_root=output_root)
    
    sketched_samples = hook.sketches.sketch_many(samples, clear_old=overwrite, workers=workers)
    logger.info("Sketched %d metagenomes into output root %s.", len(sketched_samples), output_root)


@classify_app.command(help="""
Screen current metagenome sketches against database.
""")
def screen(
    output_root: Annotated[Path, typer.Option("-o", "--output", help="Classification output path")],
    index_path: Annotated[Path, typer.Option("-X", "--sourmash-index", help="Inverted index path")],
    workers: Annotated[int, typer.Option("-w", "--workers", help="Number of Sourmash search threads.")],
):
    from maki.models.analyse.hook import OutputHook
    
    hook = OutputHook(output_root=output_root)


def main():
    SystemExit(classify_app())


if __name__ == "__main__":
    main()
