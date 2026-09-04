
from __future__ import annotations
from concurrent.futures import ProcessPoolExecutor
import gzip
from itertools import repeat
import logging
import os
from pathlib import Path
import re
import shutil
from typing import Iterable

from maki.models.database.manifest import GenomeData
from maki.utils.io import open_maybe_gzip, read_fasta
import maki.core as mx

from .annotations import AnnotationRecord


logger = logging.getLogger(__name__)

ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
DBXREF_RE = re.compile(r"(?:^|;)Dbxref=([^;]+)")


def _record_worker(args: tuple[Path, GenomeData]):
    source_dir, record = args
    return _try_write_record(source_dir, record)


def _try_write_record(source_dir: Path, record: GenomeData) -> list[AnnotationRecord]:
    ext = "gff" if record.gff else "fna"
    dest = source_dir / f"{record.accession}.{ext}.gz"
    
    if dest.exists():
        logger.warning("skipping writing %s as it already exists", record.accession)
        return []
    
    else:
        tmp_dest = dest.with_suffix(".tmp")
        
        try:
            annots = _write_record__dispatch(record, tmp_dest)
            os.replace(tmp_dest, dest)
            
        finally:
            tmp_dest.unlink(missing_ok=True)
            
        return annots


def _write_record__dispatch(record: GenomeData, file: Path) -> list[AnnotationRecord]:
    if not record.fasta and not record.gff:
        raise RuntimeError(
            f"No data suppled for accession {record.accession}"
        )
    
    elif record.fasta and not record.gff:
        return _write_fasta_as_gff(record, file)
        
    else:
        missing_annots, annots = _write_gff_to(record, file)
        
        if missing_annots:
            logger.warning(
                "GFF file %s contains %d annotations without ID attribute. These "
                "regions will not be indexed.",
                os.path.basename(record.gff), # type: ignore
                missing_annots
            )
        
        return annots


def _write_fasta_as_gff(
    record: GenomeData,
    file: Path,
    per_genome: bool = True
) -> list[AnnotationRecord]:
    assert record.fasta
    
    logger.warning(
        "Writing fasta %s as dummy GFF, with IDs per %s.",
        record.fasta,
        "genome" if per_genome else "contig"
    )
    
    annots: list[AnnotationRecord] = []
    
    contigs = list(read_fasta(record.fasta))
    
    with gzip.open(file, "wt") as fh:
        fh.write("##gff-version 3\n")

        for contig_id, seq in contigs:
            length = len(seq)

            accession = record.accession if per_genome else contig_id

            fh.write(
                f"##sequence-region {contig_id} 1 {length}\n"
            )

            fh.write(
                f"{contig_id}\tdummy\tregion"
                f"\t1\t{length}\t.\t+\t.\tID=base\n"
            )

            fh.write(
                f"{contig_id}\tdummy\tregion"
                f"\t1\t{length}\t.\t-\t.\tID=base\n"
            )
            
            annots.append(
                AnnotationRecord(accession, f"{contig_id}-base", None, None)
            )

        fh.write("##FASTA\n")

        for contig_id, seq in contigs:
            fh.write(f">{contig_id}\n")

            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")
    
    return annots


def _gff_get_accn(line: str):
    return line[:line.index("\t")]


def _gff_get_attrs(line: str):
    attrs = line[(line.rindex("\t")+1):]
    
    feature_id: str | None = None
    m = ID_RE.search(attrs)
    if m:
        feature_id = m.group(1)

    dbxrefs: list[tuple[str, str]] = []
    m = DBXREF_RE.search(attrs)
    if m:
        dbxrefs = [
            tuple(kvp.split(":", 1))
            for kvp in m.group(1).split(",")
        ] # type: ignore
        
    return feature_id, dbxrefs


def _parse_gff_attributes(line: str):
    try:
        accn = _gff_get_accn(line)
        annot_id, xrefs = _gff_get_attrs(line)
    except ValueError:
        print(line)
        raise
    
    seed = None if annot_id is None else f"{accn}-{annot_id}"
    return seed, xrefs


def _write_gff_to(
    record: GenomeData,
    file: Path
) -> tuple[int, list[AnnotationRecord]]:
    assert record.gff
    
    annots: list[AnnotationRecord] = []
    missed_annots = 0
    
    with gzip.open(file, "wt") as f:
        
        # Write GFF component
        with open_maybe_gzip(Path(record.gff), "rt") as gff:
            for line in gff:
                if "\tbakta\tregion\t" in line:
                    f.write(f"#{line}")
                    
                elif line.startswith("##FASTA"):
                    if record.fasta:
                        logger.warning(
                            "Record supplies %s sequence in both "
                            "GFF and FNA, preferring GFF.",
                            record.accession
                        )
                    
                    f.write(line)
                    f.writelines(gff)
                    
                    return missed_annots, annots
                
                else:
                    if not line.startswith("#"):
                        seed_name, xrefs = _parse_gff_attributes(line.strip())
                        
                        if seed_name:
                            for xname, xlabel in xrefs:
                                annots.append(
                                    AnnotationRecord(record.accession, seed_name, xname, xlabel)
                                )
                        
                        else:
                            missed_annots += 1
                    
                    f.write(line)
        
        # Write remaining FASTA component
        if not record.fasta:
            raise RuntimeError(f"No FASTA sequence supplied for {record.accession}")
        else:
            f.write("##FASTA\n")
            
            with open_maybe_gzip(Path(record.fasta), "rt") as fna:
                f.writelines(fna)
    
    return missed_annots, annots


class SequenceSourceDir:
    def __init__(self, source_dir: Path) -> None:
        self.source_dir = source_dir
        self.source_dir.mkdir(parents=True, exist_ok=True)
    
    def insert_genome(self, record: GenomeData) -> list[AnnotationRecord]:
        return _try_write_record(self.source_dir, record)
    
    def insert(self, data: Iterable[GenomeData]) -> list[AnnotationRecord]:
        return [
            annot
            for record in data
            for annot in self.insert_genome(record)
        ]
    
    def pinsert(self, data: Iterable[GenomeData], concurrency: int) -> list[AnnotationRecord]:
        tasks = zip(repeat(self.source_dir), data)
        
        with ProcessPoolExecutor(max_workers=concurrency) as executor:
            new_annots = [
                annot
                for annot_list in executor.map(_record_worker, tasks)
                for annot in annot_list
            ]
        
        return new_annots
        
    def sources(self):
        for p in self.source_dir.iterdir():
            if p.name.endswith(".fna.gz") or p.name.endswith(".gff.gz"):
                yield p
        
    def accessions(self):
        for p in self.sources():
            yield p.name.rsplit(".", 2)[0]
    
    def write_manifest(self, fh):
        for p in self.sources():
            fh.write(f"{p}\n")
        fh.flush()


class PreIndex(SequenceSourceDir):
    def __init__(self, root: Path) -> None:
        self.root = root
        self.root.parent.mkdir(parents=True, exist_ok=True)
        
        self.sqlite_file = self.root.with_suffix(".sql")
        
        SequenceSourceDir.__init__(
            self,
            self.root.parent / f"{self.root.name}-source"
        )
    
    def clear_sources(self):
        shutil.rmtree(self.source_dir)
        
    def write_digest(self, digest: str):
        if not self.root.exists():
            raise RuntimeError(f"Index does not exist at {self.root}")
        
        (self.root / "digest").write_text(digest)
        
    def _temp_build_file_locs(self):
        tmp_build_dir = self.root.parent / f".{self.root.name}-temp"
        tmp_manifest_file = self.root.parent /  f".{self.root.name}-temp-manifest.txt"
        return tmp_build_dir, tmp_manifest_file
    
    def cleanup_temp(self):
        tmp_build_dir, tmp_manifest_file = self._temp_build_file_locs()
        
        if tmp_build_dir.exists():
            shutil.rmtree(tmp_build_dir)
        
        tmp_manifest_file.unlink(missing_ok=True)
    
    def build(self, k: int, s: int, threads: int):
        tmp_build_dir, tmp_manifest_file = self._temp_build_file_locs()
        tmp_build_dir.mkdir(parents=True)
        
        # Write build manifest
        with tmp_manifest_file.open("wt") as fh:
            self.write_manifest(fh)
        
        # Build
        try:
            manifest = mx.read_manifest(str(tmp_manifest_file), mx.FileType.GFF3)
            opts = mx.build_opts(k, s, tmp_build_dir, threads)
            
            mx.construct_cdbg(manifest, opts)
            
            # Finalise
            os.replace(tmp_build_dir, self.root)

        except Exception:
            self.cleanup_temp()
            
            raise
