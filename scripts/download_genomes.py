#!/usr/bin/env python3

"""Download genomes from a manifest.

To create the environment:

```{console}
$ mamba create -n ncbi-datasets-cli conda-forge::ncbi-datasets-cli
```
"""

from argparse import ArgumentParser
from collections import defaultdict
import gzip
import io
from itertools import repeat
from multiprocessing import Pool
import os
import re
import shutil
import subprocess
from tempfile import NamedTemporaryFile
from typing import List, Mapping, Tuple

ACC_COL = 0
SRC_COL = 1
TYPE_COL = 2
START_COL = 3
END_COL = 4
SCORE_COL = 5
STRAND_COL = 6
PHASE_COL = 7
ATTR_COL = 8

GffRecordType = Tuple[str, str, str, int, int, str, str, str, str]

def parse_args():
    parser = ArgumentParser()
    parser.add_argument("-f", "--inputfile", type=str, required=True, help="File containing genomes to download")
    parser.add_argument("-o", "--outdir", type=str, default=None, help="Where to save downloaded genomes")
    parser.add_argument("-p", "--procs", dest="num_procs", type=int, default=5, help="Number of filtering processes")
    parser.add_argument("--api-key", dest="api_key", type=str, default=None, help="NCBI api key")
    parser.add_argument("--skip-download", dest="skip_download", action="store_true", default=False, help="Format only")
    return parser.parse_args()

def parse_manifest(file: str) -> Tuple[List[str], List[List[str]]]:
    with open(file, 'r') as f:
        header = f.readline().strip().split(",")
        data = [line.strip().split(",") for line in f]
        return header, data

def list_dirs(base: str) -> List[str]:
    return [f.path for f in os.scandir(base) if f.is_dir()]

def yield_chunks(seq: str, chunk_size: int):
    i = 0
    while i < len(seq):
        yield seq[i:i+chunk_size]
        i += chunk_size

def strip_compressed_extensions(file: str) -> str:
    if file.endswith(".gz"):
        return file[:-3]
    else:
        return file

def is_fasta(file: str) -> bool:
    return file.endswith(".fa") or file.endswith(".fna") or file.endswith(".faa") or file.endswith(".fasta")

def get_accessions_from_archive(md5file: str):
    print("[get_accessions_from_archive] reading accessions from %s" % md5file)
    with open(md5file, 'r') as f:
        for line in f:
            md5, fname = line.strip().split()
            genome_name = re.findall("ncbi_dataset/data/([^/]+)/.*", fname)
            if genome_name:
                yield genome_name[0]
            
class ArchiveDoesNotExist(Exception):
    pass

def download_by_accession(outdir: str, *accessions: str, api_key = None) -> None:
    with NamedTemporaryFile(mode='w', suffix=".txt") as f:
        f.writelines("\n".join(accessions))
        f.flush()

        output_file = os.path.join(outdir, "ncbi_dataset.zip")
        output_archive = os.path.join(outdir, "gff3")
        ncbi_base = os.path.join(output_archive, "ncbi_dataset", "data")

        try:
            if os.path.exists(output_archive):
                # files have been downloaded previously, check all accessions are accounted for

                archive_md5 = os.path.join(output_archive, "md5sum.txt")
                
                query_set = set(accessions)
                ref_set = set(get_accessions_from_archive(archive_md5))
                print("[download_by_accession] found %d md5sum values for reference genomes" % len(ref_set))

                remaining_queries = query_set - ref_set
                if remaining_queries:
                    raise ArchiveDoesNotExist()
            else:
                raise ArchiveDoesNotExist()
        except ArchiveDoesNotExist:
            # download genomes

            download_cmd = [
                "datasets", "download", "genome", "accession",
                "--annotated",
                "--dehydrated",
                "--inputfile", f.name,
                "--include", "genome,gff3",
                "--filename", output_file
            ]
            if api_key is not None:
                download_cmd.extend(["--api-key", api_key])

            print("[download_by_accession] downloading dehydrated data")
            subprocess.check_call(download_cmd)

            # unzip archive

            unzip_cmd = ["unzip", output_file, "-d", output_archive]
            print("[download_by_accession] unzipping archive")
            subprocess.check_call(unzip_cmd)

    # rehydrate dataset

    hydrate_cmd = ["datasets", "rehydrate", "--directory", output_archive]
    print("[download_by_accession] rehydrating archive")
    subprocess.check_call(hydrate_cmd)

    return output_archive, ncbi_base

def find_assembly_files(base: str) -> Tuple[str, str]:
    fna_files = [
        f.path
        for f in os.scandir(base)
        if f.is_file() and is_fasta(strip_compressed_extensions(f.path))
    ]
    if len(fna_files) != 1:
        print("[error] assembly data folder %s contains multiple fna files, can't establish genome" % base)
        exit(1)

    gff_file = os.path.join(base, "genomic.gff")
    if not os.path.exists(gff_file):
        print("[error] could not find annotation file, expected at %s" % gff_file)
        exit(1)
    return fna_files[0], gff_file

class Contig:
    def __init__(self, header: str, sequence: str):
        self.header = header
        self.sequence = sequence

    def size(self) -> int:
        return len(self.sequence)

def load_genome(genome_file: str) -> Mapping[str, Contig]:
    data = {}

    accn = None
    header = None
    buffer = io.BytesIO()

    try:
        with open(genome_file, 'rb') as f:
            for line in f:
                line = line.strip()

                if line.startswith(b">"):
                    if accn is not None:
                        data[accn] = Contig(header, buffer.getvalue().decode("utf-8"))
                        buffer.seek(0)
                        buffer.truncate()
                    accn = re.findall(b"^>(\\S+)", line)[0].decode("utf-8")
                    header = line[1:].decode("utf-8")
                else:
                    buffer.write(line)
            else:
                data[accn] = Contig(header, buffer.getvalue().decode("utf-8"))
    finally:
        buffer.close()

    return data

def load_annotations(gff_file: str) -> List[GffRecordType]:
    with open(gff_file, 'r') as f:
        data = [line.strip().split("\t") for line in f if not line.startswith("#")]
    
    for r in data:
        r[START_COL] = int(r[START_COL])
        r[END_COL] = int(r[END_COL])
    
    return data

def annot_length(rec: GffRecordType) -> int:
    return rec[END_COL] - rec[START_COL] + 1

def check_whole_regions(annots: List[GffRecordType], fna: Mapping[str, Contig]) -> List[int]:
    return [
        i
        for i, rec in enumerate(annots)
        if annot_length(rec) == fna[rec[ACC_COL]].size()
    ]

def group_annotations(annots: List[GffRecordType]) -> Mapping[Tuple[int, int, str], List[int]]:
    groups = {}

    for i, r in enumerate(annots):
        k = (r[START_COL], r[END_COL], r[STRAND_COL])

        try:
            groups[k].append(i)
        except KeyError:
            groups[k] = [i]
    
    return groups

def select_most_features(gff: List[GffRecordType], *inds: int) -> int:
    return max(inds, key = lambda i: len(gff[i][ATTR_COL]))

def apply_gff_selection(gff: List[GffRecordType], gff_groups: Mapping[Tuple[int, int, str], List[int]]) -> List[int]:
    selected_annots = []
    
    for inds in gff_groups.values():
        # get annotation types

        rec_types = {}
        for i in inds:
            try:
                rec_types[gff[i][TYPE_COL]].append(i)
            except KeyError:
                rec_types[gff[i][TYPE_COL]] = [i]

        # selection priority
        
        if "CDS" in rec_types:
            selected_annots.append(select_most_features(gff, *rec_types["CDS"]))
        elif "gene" in rec_types:
            selected_annots.append(select_most_features(gff, *rec_types["gene"]))
        elif "pseudogene" in rec_types:
            selected_annots.append(select_most_features(gff, *rec_types["pseudogene"]))
        else:
            selected_annots.append(select_most_features(gff, *(j for subinds in rec_types.values() for j in subinds)))

    return selected_annots

def make_ids_unique(gff: List[GffRecordType], inds: List[int]) -> None:
    id_counts = defaultdict(int)

    for idx in inds:
        try:
            annot_name = re.findall("ID=([^;]+)[;\n]", gff[idx][ATTR_COL])[0]
        except IndexError as e:
            raise Exception("Caught index error at line 256, caused by gff entry: %s" % ", ".join(str(v) for v in gff[idx]))
        id_counts[annot_name] += 1
        gff[idx][ATTR_COL] = gff[idx][ATTR_COL].replace("ID=%s" % annot_name, "ID=%s-%d" % (annot_name.replace(" ", "_"), id_counts[annot_name]), 1)

class CoverageStats:
    num_annots = 0
    covered_length = 0
    coding_length = 0

class CoverageStatSummary:
    def __init__(self, genome: str, accn: str, size: int, annots: int, cov: float, coding: float, annot_size: float):
        self.genome = genome
        self.accn = accn
        self.size = size
        self.annots = annots
        self.cov = cov
        self.coding = coding
        self.annot_size = annot_size

def calculate_coverage_stats(genome_name: str, fna: Mapping[str, Contig], gff: List[GffRecordType], gff_inds: List[int]):
    # count gff length

    gff_cov = defaultdict(CoverageStats)

    for i in gff_inds:
        gff_cov[gff[i][ACC_COL]].num_annots += 1
        gff_cov[gff[i][ACC_COL]].covered_length += annot_length(gff[i])
        if gff[i][TYPE_COL] == "CDS":
            gff_cov[gff[i][ACC_COL]].coding_length += annot_length(gff[i])

    # calculate coverage

    for accn, contig in fna.items():
        stats = gff_cov[accn]
        yield CoverageStatSummary(genome_name, accn, contig.size(), stats.num_annots, stats.covered_length / contig.size(), stats.coding_length / contig.size(), 0.0 if (stats.num_annots == 0) else stats.covered_length / stats.num_annots)

def filter_annotations(fna_file: str, gff_file: str, outdir: str) -> List[CoverageStatSummary]:
    genome_name = re.findall("/?([^/]+)/genomic\\.gff$", gff_file)[0]

    fna = load_genome(fna_file)
    gff = load_annotations(gff_file)

    inds_to_remove = check_whole_regions(gff, fna)
    inds_to_remove.sort()

    print("[main] removing %d annotation(s) from %s as they cover the whole contig" % (len(inds_to_remove), genome_name))
    for i in inds_to_remove[::-1]:
        gff.pop(i)

    gff_groups = group_annotations(gff)
    inds = apply_gff_selection(gff, gff_groups)
    inds.sort()
    
    print("[main] selected %d/%d annotations for %s" % (len(inds), len(gff), genome_name))

    # make IDs unique
    
    make_ids_unique(gff, inds)

    # calculate coverage stats

    stats = list(calculate_coverage_stats(genome_name, fna, gff, inds))

    output_file = os.path.join(outdir, genome_name + ".gff.gz")
    print("[main] writing output to %s" % output_file)

    with gzip.open(output_file, 'wb') as f:
        for gff_i in inds:
            f.write(b"%s\n" % (b"\t".join(str(v).encode("utf-8") for v in gff[gff_i])))
        f.write(b"##FASTA\n")
        for accn in sorted(fna.keys()):
            f.write(b">%s\n" % fna[accn].header.encode("utf-8"))
            for seq_chunk in yield_chunks(fna[accn].sequence, 60):
                f.write(b"%s\n"% seq_chunk.encode("utf-8"))

    return stats

if __name__ == "__main__":
    if not shutil.which("datasets"):
        print("[error] cannot find `NCBI datasets` tool - download from https://github.com/ncbi/datasets")
        exit(1)

    args = parse_args()

    # parse manifest
    header, manifest = parse_manifest(args.inputfile)

    try:
        aidx = header.index("Accession")
    except ValueError:
        print("[error] manifest does not contain `Accession` column")
        exit(1)

    queries = {r[aidx] for r in manifest}
    print("[main] parsed %d queries accessions" % len(queries))

    # run download command

    outdir = os.path.dirname(args.inputfile) if args.outdir is None else args.outdir
    
    if not args.skip_download:
        archive_dir, ncbi_base = download_by_accession(outdir, *queries, api_key=args.api_key)
    else:
        archive_dir = os.path.join(outdir, "gff3")
        ncbi_base = os.path.join(archive_dir, "ncbi_dataset", "data")

        if not os.path.exists(archive_dir):
            print("[error] opted to skip download, but archive %s does not exist" % archive_dir)
            exit(1)
        elif not os.path.exists(ncbi_base):
            print("[error] opted to skip download, but dataset %s does not exist" % ncbi_base)
            exit(1)
        else:
            print("[main] skipping data download")

    # filter accessions

    filtered_out = os.path.join(archive_dir, "filtered")
    filtered_stats_file = os.path.join(archive_dir, "filtered.stats.txt")
    if not os.path.exists(filtered_out):
        os.mkdir(filtered_out)

    assembly_folders = list_dirs(ncbi_base)

    # create work queue

    tasks = [find_assembly_files(base) for base in assembly_folders]

    with Pool(processes=args.num_procs) as ppl:
        stats = [
            l
            for result in ppl.starmap(filter_annotations, zip(*zip(*tasks), repeat(filtered_out)))
            for l in result
        ]

    with open(filtered_stats_file, 'w') as f:
        f.write("genome\taccession\tcontig_length\tnum_annots\tannot_coverage\tcoding_coverage\taverage_annot_size\n")
        for s in stats:
            f.write("%s\t%s\t%d\t%d\t%.4f\t%.4f\t%.1f\n" % (s.genome, s.accn, s.size, s.annots, s.cov, s.coding, s.annot_size))
    print("[main] wrote stats to %s" % filtered_stats_file)
